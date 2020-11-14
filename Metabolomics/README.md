# Analysis of untargeted metabolomics from intestinal environment of rainbow trout (_Oncorhynchus mykiss_)

__Author:__ Jacob Agerbo Rasmussen

__Contact:__ <genomicsisawesome@gmail.com>

__Date:__ 09.nov.2020

## Data availability

* Processed data for redoing analysis is available in underlying repository. 

* The metabolomics datasets generated and analysed during the current study is available in the MSV000084364 repository, ftp://massive.ucsd.edu/MSV000084364/.

* Raw spectres for metabolomics were processed through GNPS
**  Guidelines for process of raw spectres can be found at https://ccms-ucsd.github.io/GNPSDocumentation/gnpsanalysisoverview/
**  Please refer to main text for details of GNPS processing

__________________________________________________________________________________________________________
###  load dependencies
```{r}
library(ggplot2)
library(readxl)
library(phyloseq)
library(ape)
library(dplyr)
library(MetaboDiff)
library(RColorBrewer)
library(cowplot)
library(EnhancedVolcano)
library(ggpubr)
library(rstatix)
library(beeswarm)
library(cluster)
```
### Load data
```{r}
ft <- read_excel("featuretable.xlsx")
tax <- read_excel("Chemical_Class.xlsx")
md <- read_excel("metadata.xlsx")
MetDNA <- read.csv("HF_pseudomonas_MRN.annotation.result.csv")

# Create colour palette
My_pal <- c("#3889A0","#F2AD00","#9C964A","#D95F02","#1B9E77","#046C9A","#0B775E","#35274A","#F2300F","#666666") 
group_pal <- c("#e3aa74", "#ed828c", "#7bb6bd")
```

### Data filtering
To minimise false positives we included procedural blanks in the data generation. Metabolites which are present in the blanks should be removed (at least to my believe).
To include a little bit of plasticity in filtering of blanks i include only metabolites, which are lower than 5,000,000 in blanks.
```{r remove metabolites present in blanks}
## Remove features present in blanks
blk_ids <- which(colnames(ft) %in% as.character(md$SampleID[grep("BLANK",md$SampleType )]))
blanks <- ft[,blk_ids]
blanks <- blanks[apply(blanks[,-1], 1, function(x) !all(x < 1*10^5)),]
blank_ids <- rownames(blanks)
ft.filtered <- ft[,-c(blk_ids)]
ft.filtered = ft.filtered[-which(rownames(ft.filtered) %in% blank_ids), ]
# remove QCs
QC_ids <- which(colnames(ft.filtered) %in% as.character(md$SampleID[grep("QC",md$SampleType )]))
ft.filtered_all <- ft.filtered[,-c(QC_ids)]
ft <- ft.filtered_all
# set ID as rownames 
ft <- as.matrix(ft)
rownames(ft) <- ft[,1]
ft <- ft[,-c(1)]
```
To minimise highly zero inflated dataset, I remove metabolites which are present in less than 50 of the 60 samples.
```{r remove metabolites present in blanks}
# Remove metabolites not present in 50 of 60 of the samples
ft_curated <- ft[apply(ft, 1, FUN = function(x){sum(x == 0)}) < 10,]
# post blank removal
dim(ft)
# post inflation removal
ft <- ft_curated
```
Make sure metadata and annotation follows
```{r Curate tax adn md files}
## check if datasets are identically ordered
## Due to a bug, there is one metabolite no. "124151", which is NA in tax
ft <- ft[rownames(ft) != "125141",]
## Filter tax and MetDNA annotation files of remove metabolites
tax <- tax[match(rownames(ft),tax$ID),]
MetDNA <- MetDNA[match(rownames(ft),MetDNA$name),]

#combine tax and MetDna
Annotation <- as.data.frame(tax)

cmp <- MetDNA$compound.name
Annotation$Annotation <- paste(Annotation$Compound_Name,cmp)
## remove "NA NA" and make it to NA, so pasted text looks better
Annotation$Annotation <- gsub("NA NA", NA, as.character(Annotation$Annotation))
Annotation$Annotation <- gsub("NA ", "", as.character(Annotation$Annotation))
Annotation$Annotation <- gsub(" NA", "", as.character(Annotation$Annotation))

identical(tax$ID,as.numeric(rownames(ft)))
identical(MetDNA$name,as.integer(rownames(ft)))
identical(Annotation$ID,as.numeric(rownames(ft)))
## Cool seems to be identical

#Filter out blanks and QC from metadata file
md <- md[match(colnames(ft),md$SampleID),]
rownames(md) <- md$SampleID
identical(rownames(md),colnames(ft))
rownames(md) <- colnames(ft) # to minimise difficulties with phyloseq
md <- as.matrix(md) # to minimise difficulties with phyloseq
md <- as.data.frame(md) # to minimise difficulties with phyloseq - i dont know what is going on?! Probably a bug with class in rownames
```
I tend to clean up my working space, because i'm one of those who get lost in my own r codes....
```{r Clean up - some peepz just hate the mess...}
rm("blanks", "ft.filtered","ft.filtered_all", 
   "blank_ids", "blk_ids", "QC_ids","ft_curated", "cmp")
```
I use phyloseq to order my data (both 16S and metabolomics), but anything could be used. :)
```{r generate phyloseq elements, because phyloseq is genious for this type of abundance data}
# physeq element
# since Phyloseq likes taxa to be in matrices, i will make data to matrix
Ant <- as.matrix(Annotation)
rownames(Ant) <- rownames(ft)

# combine procedural replicates 
#agg = aggregate(t(ft),
#                by = list(md$Individuals),
#                FUN = sum)
#rownames(agg) <- agg$Group.1
#agg <- agg[,-c(1)]
#ft <- t(agg)

## Make metadata compatible to metabolite data
#md <- md[md$RUN == "B",]
#colnames(ft) <- md$SampleID

#make phyloseq object for ordination
physeq <- phyloseq(otu_table(ft, taxa_are_rows = TRUE), tax_table(Ant),sample_data(md))
# sum normalise physeq object for ranked analysis
physeq.norm = transform_sample_counts(physeq, function(x) x/sum(x))
#make dataset only for metabolites versus mycoplasma
Myco_x_meta = subset_samples(physeq, Keep=="Yes")
```
We start with a simple rank-abundance barplot, using the cumulative fractional abundance of each metabolite in the dataset.

This really visualise the sparsity of the dataset. Simply a lot of metabolites, which are very low abundant.
```{r barplot stuff,message=FALSE}
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
barplot(sort(taxa_sums(physeq.norm), TRUE)[1:740]/nsamples(physeq.norm), las=2, col = "#1B9E77")
```
![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin/cumulative_plot.pdf)

### Metabolites Ordination
I use PCoA to ordinate variantion of metabolites across all metabolites, amino acids, terpenoids, bile acids.

I am using jaccard distances to be conservative about the intensity measurements and i would like to see how composition is determined by presence/absence of metabolites. 
```{r Class specific Ordination}
GP.ord <- ordinate(Myco_x_meta, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta, GP.ord, type="samples", color="Feed", axes = c(1,2)) 
md_ord <- md[md$Keep == "Yes",]
p_ALL =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)
  
Myco_x_meta_Class <- subset_taxa(Myco_x_meta, CF_subclass=="Amino acids, peptides, and analogues")
GP.ord <- ordinate(Myco_x_meta_Class, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta_Class, GP.ord, type="samples", color="Feed", axes = c(1,2)) 
md_ord <- md[md$Keep == "Yes",]
p_AA =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)

Myco_x_meta_Class <- subset_taxa(Myco_x_meta, CF_subclass=="Monoterpenoids" | CF_subclass=="Diterpenoids" |CF_subclass=="Triterpenoids" | CF_subclass=="Sesquiterpenoids")
GP.ord <- ordinate(Myco_x_meta_Class, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta_Class, GP.ord, type="samples", color="Feed", axes = c(1,3)) 
md_ord <- md[md$Keep == "Yes",]
p_terp =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)

Myco_x_meta_Class <- subset_taxa(Myco_x_meta, CF_subclass=="Bile acids, alcohols and derivatives")
GP.ord <- ordinate(Myco_x_meta_Class, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta_Class, GP.ord, type="samples", color="Feed", axes = c(1,2)) 
md_ord <- md[md$Keep == "Yes",]
p_Bile =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)

#pdf("jaccard_ordination_metabolites_x_mycoplasma_classes.pdf", height = 8, width = 12)
cowplot::plot_grid(p_ALL,
                   p_AA,
                   p_terp,                   
                   p_Bile,
                   labels = c('All Classes of metabolites','Amino acids, peptides, and analogues', 'Terpenoids',
                              'Bile acids, alcohols and derivatives'), ncol = 2)
#dev.off()
```
Clean up - some peepz just hate the mess...
```{r}
rm("Myco_x_meta_Class","md_ord", 
   "MetDNA", "tax","physeq.norm", "p_Bile", "p_FA","p_AA","p_terp","p2")
```
###    Differential intensity of metabolites across samples         
We use metabodiff to analyse the dataset, since it is includes knn imputation and SVN normalisation (Mock et al. 2015). 
See also: https://github.com/andreasmock/MetaboDiff. 
```{r - }
met <- create_mae(ft,Ant,md)

met = knn_impute(met,cutoff=0.4)
met <- normalize_met(met)
quality_plot(met,
             group_factor="Feed",
             label_colors=group_pal)
met = diff_test(met,
                group_factors = c("Feed"))             
```
![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin//QC_plot_metabodiff.pdf)

Despite that Metabodiff can du volcano plots, i like to use EnhacedVolcano Libs, which i find really cool :) 
```{r plot differential test, using Enhanced volcano}
# Extract data from test to use for selection of significantly different metabolites for heatmap
df <- met@metadata[["anova_Feed_Control_vs_Probiotics_vs_Synbiotics"]]
df$Metabolites <- as.integer(rownames(met@ExperimentList@listData[["raw"]]@assays@data@listData[[1]]))
tax_diff <-  Ant[order(match(df$Metabolites,rownames(Ant))),]
df <- cbind(df,tax_diff)

# Intercept is CTRL -----> PRO & SYN
p1 = EnhancedVolcano(df,
                     lab = c(ifelse(df$adj_pval<10e-3, 
                                    ifelse(df$CF_class != "no matches", 
                                           df$CF_class, ""), "")),
                     boxedLabels = TRUE,
                     col = c('grey0', 'grey0', 'red3', 'red3'),
                     x = 'fold_change',
                     y = 'adj_pval',
                     xlim = c(-10,10),
                     pCutoff = 0.05,
                     FCcutoff = 1.0,
                     labSize = 3,
                     pointSize = c(ifelse(df$adj_pval<10e-3, 5, 1)),
                     title = 'Differential Abundance of metabolites related to feed',
                     subtitle = 'Comparison between CTRL, PRO, and SYN',
                     xlab = bquote(~Log[2]~ 'Fold Change'),
                     .legend = c('NS','P & Log2 FC'),
                     legendLabels = c('NS',expression(p-value~and~log[2]~FC)),
                     legendPosition = 'bottom')

pdf("Differential Abundance of metabolites MetaboDiff template.pdf", height = 7.5, width = 15)
cowplot::plot_grid(p1, ncol = 1)
dev.off()

#write.csv(df,file = "Diff_test_Metabodiff.csv")
#write.csv(df_sig,file = "Sig_Diff_test_Metabodiff.csv")
```

