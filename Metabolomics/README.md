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

## load dependencies
```r
library(readxl)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(phyloseq)
library(MetaboDiff)
library(ape)
library(RColorBrewer)
library(cowplot)
library(EnhancedVolcano)
```

```{r load data}
# Load data
ft <- read.csv("MS1_curated.csv")
rt_mz <- ft[,1:3]
rownames(ft) <- ft[,1]
ft <- ft[,4:63]

tax <- read_excel("Chemical_Class.xlsx")
md <- read.csv("metadata.csv")
MetDNA <- read.csv("HF_pseudomonas_MRN.annotation.result.csv")

My_pal <- c("#3889A0","#F2AD00","#9C964A","#D95F02","#1B9E77","#046C9A","#0B775E","#35274A","#F2300F","#666666") 
group_pal <- c("#e3aa74", "#ed828c", "#7bb6bd")
```

```{r Curate tax files}
tax <- tax[match(rownames(ft),tax$ID),]
MetDNA <- MetDNA[match(rownames(ft),MetDNA$name),]
cmp <- MetDNA$compound.name
Ant <- cbind(tax,cmp)
Ant <- as.matrix(Ant)
rownames(Ant) <- Ant[,1]
rownames(md) <- colnames(ft)
Ant <- as.data.frame(Ant)
Ant$Annotation <- paste(Ant$Compound_Name,Ant$cmp)
#is.na(Ant$ID)

# metobolite 253 seems be NA in the dataset
Ant <- Ant[-c(253),]

names <- Ant$ID
Ant <- as.matrix(Ant)
rownames(Ant) <- names
```
## Ordination Analysis
I use phyloseq for ordination analysis. To me it's very handy, when handy abundance/intensity based data together with classification and metadata.

I start with generating a phyloseq object and generate a PCoA with jaccard based distances, since this is a conservative measurement 
to minimise biases of relative intensity.

I included abundance of _Mycoplasma_ in metadata to visualise how feeding type affects both _Mycoplasma_ abundance and the metabolic landscape. 
```{r generate phyloseq elements, because phyloseq is genious for this type of abundance data}
# physeq element
physeq <- phyloseq(otu_table(ft, taxa_are_rows = TRUE), tax_table(Ant),sample_data(md))
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq = merge_phyloseq(physeq, random_tree)
physeq.sort = transform_sample_counts(physeq, function(x) x/sum(x))
Myco_x_meta = subset_samples(physeq.sort, Keep=="Yes")
```
```{r PCOA plot - sum normalised}
#r NMDS/PCOA plot - sum normalised
MM.ord <- ordinate(Myco_x_meta, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta, MM.ord, type="samples", color="Feed") 
#pdf("jaccard_ordination_metabolites_x_mycoplasma.pdf", height = 8, width = 12)
md_ord <- md[md$Keep == "Yes",]
p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)
#dev.off()
```
![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin/PCoA_Myco_Meta.png)

## Class specific Ordination
Amino acids, lipids bile acids, and stereoids are highly relevant for high performing predators, like _O. mykiss_. 
Therefore, we look in to the variation of these classes of metabolites.
```{r}
Myco_x_meta_Class <- subset_taxa(Myco_x_meta, CF_subclass=="Amino acids, peptides, and analogues")
GP.ord <- ordinate(Myco_x_meta_Class, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta_Class, GP.ord, type="samples", color="Feed") 
md_ord <- md[md$Keep == "Yes",]
p_AA =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)

Myco_x_meta_Class <- subset_taxa(Myco_x_meta, CF_superclass=="Lipids and lipid-like molecules")
GP.ord <- ordinate(Myco_x_meta_Class, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta_Class, GP.ord, type="samples", color="Feed") 
md_ord <- md[md$Keep == "Yes",]
p_Lipids =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)

Myco_x_meta_Class <- subset_taxa(Myco_x_meta, CF_subclass=="Bile acids, alcohols and derivatives")
GP.ord <- ordinate(Myco_x_meta_Class, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta_Class, GP.ord, type="samples", color="Feed") 
md_ord <- md[md$Keep == "Yes",]
p_Bile =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)

Myco_x_meta_Class <- subset_taxa(Myco_x_meta, CF_class =="Steroids and steroid derivatives")
GP.ord <- ordinate(Myco_x_meta_Class, "PCoA", "jaccard")
p2 = plot_ordination(Myco_x_meta_Class, GP.ord, type="samples", color="Feed") 
md_ord <- md[md$Keep == "Yes",]
p_Steroids =p2 + geom_point(size=7.5, alpha = as.numeric(md_ord$Mycoplasma)) + 
  theme_minimal() + scale_fill_manual(values = group_pal) +
  scale_color_manual(values = group_pal)

# Remove hashtag below to write out output
#pdf("jaccard_ordination_metabolites_x_mycoplasma_classes.pdf", height = 8, width = 12)
cowplot::plot_grid(p_AA,
                   p_Lipids,
                   p_Bile,
                   p_Steroids,
                   labels = c('Amino acids, peptides, and analogues', 'Lipids and lipid-like molecules',
                              'Bile acids, alcohols and derivatives','Steroids and steroid derivatives'), ncol = 2)
#dev.off()
```
![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin/PCoA_Myco_Meta_Class.png)

## Differential intensity of metabolites across feeding types

This is a home knitted version of a differential test for relative intensity of metabolites, using non paramental kruskal wallis, comparing mean of metabolites between feeding groups

```{r - Differential Abundance testing - Categorical - w/ Kruskal Wallis}
library("beeswarm")
ft.diff <- t(ft)
md.diff <- md

### Choose category for differential abbundance testing
cats <- droplevels(as.factor(md.diff$Feed_Additives))
cats
ftt <- t(apply(ft.diff, 1, function(x) (x)/sum((x))))
ftt[1:5,1:5]

### Define function
krusk <- function (x) {
  out <- tryCatch(unlist(kruskal.test(scale(ftt[,x])[,1] ~ as.factor(cats))[c("p.value")]), error = function(e) return(NA))
  return(out)
}
corres <- t(sapply(1:ncol(ftt), krusk ))
corres <- cbind(unlist(corres[1,]),p.adjust(unlist(corres[1,]), method = "fdr"))
rownames(corres) <- colnames(ftt)
corres <- cbind(corres,rownames(corres))
colnames(corres) <- c('p.value','p.value.corrected','feature id')

#write.table(corres, 'DifferentialAbundance.txt', sep = '\t', quote = F, row.names = F)
```
Only significant adjusted p-values were chosen and boxplots were created for annotated metabolites
```
length(which(as.numeric(corres[,2]) < 0.05))
corressig <- corres[which(as.numeric(corres[,2]) < 0.05),]

ftplot <- ftt

par(mfrow = c(2,3))
  for (i in 1:length(rownames(corressig))){
    selid <- rownames(corressig)[i]
    comp <- ftplot[,which(colnames(ftplot) == selid)]
    
    df <- cbind(as.character(md.diff$Feed),comp)
    colnames(df) <- c('SampleGroup','Summed Precursor Ion Intensities/MS2 Ion')
    df <- as.data.frame(df)
    df$SampleGroup <- as.character(df$SampleGroup)
    df$`Summed Precursor Ion Intensities/MS2 Ion` <- as.numeric(as.character(df$`Summed Precursor Ion Intensities/MS2 Ion`))
  }
#dev.off()

sig_tax <- Ant[order(match(rownames(corressig), Ant[,1])),]
st <- sig_tax[,c(1,8,9)]

data_sig <- cbind(corressig,st)

data <- cbind(corressig,st)
data <- as.data.frame(data)
data$Annotation <- paste(data$Compound_Name,data$cmp)
data <- data[data$Annotation != "NA NA",]
krusk_data <- data

plot_list = list()
for (i in 1:length(rownames(data))){
    selid <- rownames(corressig)[i]
    comp <- ftplot[,which(colnames(ftplot) == selid)]
    
    df <- cbind(as.character(md.diff$Feed),comp)
    colnames(df) <- c('SampleGroup','Summed Precursor Ion Intensities/MS2 Ion')
    df <- as.data.frame(df)
    df$SampleGroup <- as.character(df$SampleGroup)
    df$`Summed Precursor Ion Intensities/MS2 Ion` <- as.numeric(as.character(df$`Summed Precursor Ion Intensities/MS2 Ion`))
    
    #boxplot(`Summed Precursor Ion Intensities/MS2 Ion` ~ SampleGroup, data = df, outline = FALSE, main = data$Annotation[i], col = group_pal)
    #beeswarm(`Summed Precursor Ion Intensities/MS2 Ion` ~ SampleGroup, data = df, col = "black", pch = 16, add = TRUE)
    # Add p-values onto the box plots
    stat.test <- df %>%
      group_by("SampleGroup") %>%
      wilcox_test(`Summed Precursor Ion Intensities/MS2 Ion` ~ SampleGroup) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj")
    stat.test <- stat.test %>%
      add_x_position(x = "SampleGroup", dodge = 0.8) %>%
      add_y_position()
    
    # Create a box plot
    p = ggboxplot(
      df, x = "SampleGroup", y = "Summed Precursor Ion Intensities/MS2 Ion", 
      color = "black",
      fill = "SampleGroup", palette = group_pal,
      outlier.shape = 8,
      size = 0.5,
      title = data$Annotation[i])  + 
      stat_pvalue_manual(
        stat.test,  label = "{p.adj} {p.adj.signif}", tip.length = 0.02,
        step.increase = 0.075, 
        bracket.nudge.y = max(df$`Summed Precursor Ion Intensities/MS2 Ion`))
    
    ## Save  output in lists
    plot_list[[i]] = p
}

#Remove hashtag below to write out output
#pdf("Boxplot_of_significant_compounds_w_sig_notification.pdf")
plot_list[5]
#dev.off()
```

![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin/Acetyl_Ornithine_boxplot.png)


```{r Generate summary table}
## Generate summary of wilcoxon statistics of metabolites between feeding groups - same stats, which are shown in boxplot figure
stat_list_kw = list()
for (i in 1:length(rownames(data))){
    selid <- rownames(corressig)[i]
    comp <- ftplot[,which(colnames(ftplot) == selid)]
    
    df <- cbind(as.character(md.diff$Feed),comp)
    colnames(df) <- c('SampleGroup','Summed Precursor Ion Intensities/MS2 Ion')
    df <- as.data.frame(df)
    df$SampleGroup <- as.character(df$SampleGroup)
    df$`Summed Precursor Ion Intensities/MS2 Ion` <- as.numeric(as.character(df$`Summed Precursor Ion Intensities/MS2 Ion`))
    
    #boxplot(`Summed Precursor Ion Intensities/MS2 Ion` ~ SampleGroup, data = df, outline = FALSE, main = data$Annotation[i], col = group_pal)
    #beeswarm(`Summed Precursor Ion Intensities/MS2 Ion` ~ SampleGroup, data = df, col = "black", pch = 16, add = TRUE)
    # Add p-values onto the box plots
    stat.test <- df %>%
      group_by("SampleGroup") %>%
      wilcox_test(`Summed Precursor Ion Intensities/MS2 Ion` ~ SampleGroup) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj")
    stat.test <- stat.test %>%
      add_x_position(x = "SampleGroup", dodge = 0.8) %>%
      add_y_position()
    
    stat_list_kw[[i]] = stat.test
}

## Make statistic summary of wilcoxon pairwise test of metabolites between groups
stat_list_kw  <-  as.data.frame(matrix(unlist(stat_list_kw), nrow=length(unlist(stat_list_kw[1]))))
stat_list_kw <- t(stat_list_kw)
rownames(stat_list_kw) <- rownames(data)
# remove redundant columns
stat_list_kw <- stat_list_kw[,-c(1,2,3,4,5,7,9,11,13,14,15,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45)]
# Add understanbable column names
names <- c("Y","gr.1","gr.2_1","gr.2_2","n1","n2_1","n2_2","statistic.1","statistic.2_1","statistic.2_2","p.1","p.2_1","p.2_2","p.adj.1","p.adj.2_1","p.adj.2_2","signif.1","signif.2_1","signif.2_2")
colnames(stat_list_kw) <- names


#Remove hashtag below to write out output
#write.csv(stat_list_kw, "Wilcoxon_Summary.csv")
```
## Metabodiff Analysis
Use Metabodiff to apply normalisation and knn imputation. Further also apply differential test between usage of feed additives and feeding groups. 
```{r - Metabodiff normalisation}
met <- create_mae(ft,tax,md)
md$Feed_Additives
met = knn_impute(met,cutoff=0.4)
met <- normalize_met(met)
quality_plot(met,
             group_factor="Feed_Additives",
             label_colors=group_pal)
```
```{r - Metabodiff diff test, include=FALSE}
met = diff_test(met,
                group_factors = c("Feed","Feed_Additives"))
str(metadata(met), max.level=2)
```
```{r plot differential test, using Enhanced volcano}
#df <- met@metadata[["anova_Feed_CTRL_vs_Probiotics_vs_Synbiotics"]]
df <- met@metadata[["ttest_Feed_Additives_Yes_vs_No"]]

df$Metabolites <- as.integer(rownames(met@ExperimentList@listData[["raw"]]@assays@data@listData[[1]]))

tax_diff <-  tax[order(match(df$Metabolites,rownames(tax))),]
df <- cbind(df,tax)

df_sig <- df[df$adj_pval<0.5,]

p1 = EnhancedVolcano(df,
    lab = df$Compound_Name,
    boxedLabels = FALSE,
    x = 'fold_change',
    y = 'adj_pval',
    xlim = c(-5,5),
    pCutoff = 10e-2,
    FCcutoff = 1.0,
    labSize = 3,
    pointSize = c(ifelse(df$adj_pval<10e-3, 5, 1)),
    title = 'Differential Abundance of metabolites related to feed additives',
    legendPosition = 'bottom')
p1

#pdf("Differential Abundance of metabolites MetaboDiff template.pdf", height = 7.5, width = 15)
#cowplot::plot_grid(p1, ncol = 1)
#dev.off()

#write.csv(df,file = "Diff_test_Metabodiff.csv")
#write.csv(df_sig,file = "Sig_Diff_test_Metabodiff.csv")
```

![alt_text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin/MetaboDiff_Volcanoplot.png)

```{r Create scatterplots with correlations within each group}
ft_selected <- ft[match(df_sig$Metabolites,rownames(ft)),]
ftt_selected <- t(apply(ft_selected, 1, function(x) (x)/sum((x))))
ftt_selected <- as.data.frame(ftt_selected)
ftt_selected <- t(ftt_selected)
colnames(ftt_selected) <- paste("X",colnames(ftt_selected), sep = ) #Add X to metabolite names, since its and integer name

test <- cbind(ftt_selected,md)
plot_list_weigth = list()
plot_list_myco = list()
# Scatter plot colored by groups ("Species")
for (i in 1:length(colnames(ftt_selected))){
selid <- colnames(ftt_selected)[i]
title  <- as.data.frame(df_sig[df_sig$Metabolites == rownames(ft_selected)[i],])
sp <- ggscatter(test, x = selid, y = "Weight",
                color = "Feed", palette = group_pal,
                title = title$CF_Dparent ,
                size = 3, alpha = 0.6, add = "reg.line", conf.int = TRUE) +
  border() +
  stat_cor(aes(color = `Feed`),method = "pearson")
sp
xplot <- ggdensity(test, selid, fill = "Feed",
                   palette = group_pal)
yplot <- ggdensity(test, "Weight", fill = "Feed", 
                   palette = group_pal)+ rotate()
# Cleaning the plots
sp <- sp +  theme(legend.position="bottom")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
p = plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))
plot_list_weigth[[i]] = p
}

# Scatter plot colored by groups ("Species")
for (i in 1:length(colnames(ftt_selected))){
selid <- colnames(ftt_selected)[i]
title  <- as.data.frame(df_sig[df_sig$Metabolites == rownames(ft_selected)[i],])
sp <- ggscatter(test, x = selid, y = "Mycoplasma",
                color = "Feed", palette = group_pal,
                title = title$CF_Dparent ,
                size = 3, alpha = 0.6, add = "reg.line", conf.int = TRUE) +
  border() +
  stat_cor(aes(color = `Feed`),method = "pearson")
sp
xplot <- ggdensity(test, selid, fill = "Feed",
                   palette = group_pal)
yplot <- ggdensity(test, "Mycoplasma", fill = "Feed", 
                   palette = group_pal)+ rotate()
# Cleaning the plots
sp <- sp +  theme(legend.position="bottom")
yplot <- yplot + clean_theme() + rremove("legend")
xplot <- xplot + clean_theme() + rremove("legend")
# Arranging the plot using cowplot
p = plot_grid(xplot, NULL, sp, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2, 1), rel_heights = c(1, 2))
plot_list_myco[[i]] = p
}
```
![alt_text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin/Cor_Metabolites_Weight.png)
![alt_text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/Metabolomics/data/bin/Cor_Metabolites_and_Myco.png)

```{r - plot scatter plots}
# Remove hashtag below to write out output
#pdf("Sig_Metabolites_and_weight.pdf")
plot_list_weigth[10]
#dev.off()
#pdf("Sig_Metabolites_and_Mycoplasma.pdf")
plot_list_myco[10]
#dev.off()
```
Calculate Pearson correlation as plotted, between metabolites and weight for each feed group & calculate Pearson correlation as plotted, between metabolites and Mycoplasma for each feed group
```{r Pearson correlation of metabolite intensity and mycoplasma and weight, include=FALSE}
stat_list = list()
for (i in 1:length(colnames(ftt_selected))){
selid <- colnames(ftt_selected)[i]
lm_data <- test[,c(selid,"Weight","Feed")]
lm_stat.test_CTRL <- lm_data[lm_data$Feed == "CTRL",]
lm_stat.test_PRO <- lm_data[lm_data$Feed == "Probiotics",]
lm_stat.test_SYN <- lm_data[lm_data$Feed == "Synbiotics",]
CTRL <- cor.test(lm_stat.test_CTRL$Weight, lm_stat.test_CTRL[,selid], 
                    method = "pearson")
CTRL.p.value <- CTRL$p.value
CTRL.cor <- CTRL$estimate
PRO <- cor.test(lm_stat.test_PRO$Weight, lm_stat.test_PRO[,selid], 
                    method = "pearson")
PRO.p.value <- PRO$p.value
PRO.cor <- PRO$estimate
SYN <- cor.test(lm_stat.test_SYN$Weight, lm_stat.test_SYN[,selid], 
                    method = "pearson")
SYN.p.value <- SYN$p.value
SYN.cor <- SYN$estimate
stat <- cbind(selid,CTRL.p.value,CTRL.cor,PRO.p.value,PRO.cor,SYN.p.value,SYN.cor)
stat_list[[i]] = stat
}
stat_list  <-  as.data.frame(matrix(unlist(stat_list), nrow=length(unlist(stat_list[1]))))
stat_list.weight <- t(stat_list)
colnames(stat_list.weight) <- c("ID","CTRL.p.value.weight","CTRL.cor.weight","PRO.p.value.weight","PRO.cor.weight","SYN.p.value.weight","SYN.cor.weight")

####
stat_list = list()
for (i in 1:length(colnames(ftt_selected))){
selid <- colnames(ftt_selected)[i]
lm_data <- test[,c(selid,"Mycoplasma","Feed")]
lm_stat.test_CTRL <- lm_data[lm_data$Feed == "CTRL",]
lm_stat.test_PRO <- lm_data[lm_data$Feed == "Probiotics",]
lm_stat.test_SYN <- lm_data[lm_data$Feed == "Synbiotics",]
CTRL <- cor.test(lm_stat.test_CTRL$Mycoplasma, lm_stat.test_CTRL[,selid], 
                    method = "pearson")
CTRL.p.value <- CTRL$p.value
CTRL.cor <- CTRL$estimate
PRO <- cor.test(lm_stat.test_PRO$Mycoplasma, lm_stat.test_PRO[,selid], 
                    method = "pearson")
PRO.p.value <- PRO$p.value
PRO.cor <- PRO$estimate
SYN <- cor.test(lm_stat.test_SYN$Mycoplasma, lm_stat.test_SYN[,selid], 
                    method = "pearson")
SYN.p.value <- SYN$p.value
SYN.cor <- SYN$estimate
stat <- cbind(selid,CTRL.p.value,CTRL.cor,PRO.p.value,PRO.cor,SYN.p.value,SYN.cor)
stat_list[[i]] = stat
}
stat_list  <-  as.data.frame(matrix(unlist(stat_list), nrow=length(unlist(stat_list[1]))))
stat_list.myco <- t(stat_list)
colnames(stat_list.myco) <- c("ID","CTRL.p.value.Myco","CTRL.cor.Myco","PRO.p.value.Myco","PRO.cor.Myco","SYN.p.value.Myco","SYN.cor.Myco")

stat_list <- cbind(stat_list.weight,stat_list.myco)
stat_list <- as.data.frame(stat_list)
stat_list$ID <- gsub("X ", "", as.character(stat_list$ID))
stat_list <- stat_list[,-c(8)]
data <- stat_list[match(stat_list$ID,df_sig$ID),]
data <- cbind(df_sig,stat_list)
data <- data[,-c(13)]
```

Write out statistics from Metabodiif between feed addtives or not, further also Pearson correlations between significant different abundant metabolites and Host weight and mycoplasma presence.
```{r, include=FALSE}
# Remove hashtag below to write out output
#write.csv(data, "Metabolite_statistics.csv")
```
