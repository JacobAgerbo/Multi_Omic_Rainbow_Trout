# Analysis of 16S data from intestinal environment of rainbow trout (_Oncorhynchus mykiss_)

__Author:__ Jacob Agerbo Rasmussen

__Contact:__ <genomicsisawesome@gmail.com>

__Date:__ 09.nov.2020

## Data availability

* Data for redoing analysis is available in underlying repository. 

* Raw 16S .fq data is available at ENA project repository: XXXXX, when article is accepted.

* Bioinformatic pipeline can be found at https://github.com/JacobAgerbo/Data2Result/tree/master/Metabarcoding



## Setup dependencies
```r
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(iNEXT); packageVersion("iNEXT")
library(ape); packageVersion("ape")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(cowplot); packageVersion("cowplot")
library(plyr); packageVersion("plyr")
library(sjmisc); packageVersion("sjmisc")
library(data.table); packageVersion("data.table")
library(metacoder);packageVersion("metacoder")
library(hilldiv);packageVersion("hilldiv")
library(car);packageVersion("car")
library(lme4);packageVersion("lme4")
library(RColorBrewer);packageVersion("RColorBrewer")
```
I use Phyloseq to organise our data for profiling the V3-V4 region of 16S rRNA gene, thank you @Joey711!

```r
#  Import and arrange data for ASVs

ASVs <- read.csv("Curated_Table.txt", sep = ",")
tax <- read.csv("Curated_Tax.csv", sep = ",")
md <- read.csv("metadata.csv", sep = ";")
rownames(tax) <- rownames(ASVs)
rownames(md) <- md$Sample_ID

# Generate phyloseq object
physeq <- phyloseq(otu_table(ASVs,taxa_are_rows=TRUE),
               tax_table(as.matrix(tax)),
               sample_data(md))
               
# Concatenate ASVs within same genus, since we only look a genus level. 
# Proper taxanomy of any bacteria found in the environment will be analysed better with metagenomics
physeq = tax_glom(physeq, "Genus")
## Create tree for physeq to use unifrac distances for ordination
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
physeq = merge_phyloseq(physeq, random_tree)         

## Abundance of ASVs are normalised
physeq.norm = transform_sample_counts(physeq, function(x) x/sum(x))
## Rare occuring ASVs were remove to minimise inflation
physeq.sort = prune_taxa(taxa_sums(physeq.norm) > 0.02*length(colnames(ASVs)), physeq.norm)
```

## Barplots of bacterial composition across distal and mid gut of _O. mykiss_
We survey the genera of the top 10 most abundant phyla across samples, using relative abundance barplots

```{r barplot stuff,message=FALSE}
physeq_samples <- subset_samples(physeq.sort, Sample_Type=="Sample")

My_pal <- c("#E6AB02","#D95F02","#7570B3","#1B9E77","#666666")

p_sample = plot_bar(physeq_samples, "Sample_ID", fill = "Family", facet_grid = ~Gut_Section) + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values=My_pal) +
  scale_color_manual(values=My_pal) +
  labs(x = "Individuals") +
  labs(y = "Relative Abundance")

plot_grid(p_sample, labels = 'AUTO', nrow = 1)

```
![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/16S/data/bin/16S_barplot.png)

## Bacterial composition analysis, using PCoA with unifrac distances

We check the variance between sample types, using PCoA and Unifrac distances. 
Unifrac distances were chosen to include the phylogenetic distance. Several other measurements were also tested, included Bray-Curtis and Jaccard. 
Here we see a clear pattern, that CTRL are clustering alone, away from both SYN and PRO feed additives (Fig. A). Furthermore, we see that sample types don't seem to affect composition of bacteria.

```{r Ordination,message=FALSE}
GP.ord <- ordinate(physeq_samples, "PCoA", "unifrac", weighted=TRUE)
p1 = plot_ordination(physeq_samples, GP.ord, type="samples",
                     color="Feed_Type",
                     shape = "Gut_Section",
                     title="Microbiome Composition of Rainbow Trout Reared on Different Feeding Types")
p1 = p1 + facet_wrap(~Feed_Type, 1) + theme_bw()

physeq_ord <- subset_samples(physeq.sort, Sample_Type !="Blank")

sample_types.ord <- ordinate(physeq_ord, "PCoA", "unifrac", weighted=TRUE)
p2 = plot_ordination(physeq_ord, sample_types.ord, type="samples",
                     color="Feed_Type",
                     shape = "Sample_Type",
                     title="Microbiome Composition of Rainbow Trout Reared on Different Sample Types")

plot_grid(p1,p2, labels = 'AUTO', nrow = 2)
```
![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/16S/data/bin/16S_PCoA.png)

## Differential abundance of bacteria

We do a differential abundance test between feeding types and include their taxonomic relation. We use metacoder to get a higher resolution of taxonomy, sinec the phylogenetic tree is more informative than actual barplots.

```{r Create metacoder dataseets, include=FALSE}
### Create meacoder environment
TopNOTUs = names(sort(taxa_sums(physeq.norm), TRUE)[1:50])
physeq_100 = prune_taxa(TopNOTUs, physeq.norm)
#substract only samples
metacoder <- subset_samples(physeq_100, Gut_Section=="Distal Gut content")
metacoder <- parse_phyloseq(metacoder)
metacoder$data$tax_abund <- calc_taxon_abund(metacoder, data = "otu_table")
metacoder$data$tax_occ <- calc_n_samples(metacoder, "tax_abund", groups = "Feed_Type")
metacoder$data$diff_table <- compare_groups(metacoder, data = "tax_abund", cols = metacoder$data$sample_data$sample_id,
                                      groups = metacoder$data$sample_data$Feed_Type)
metacoder$data$diff_table$adjusted_p_value <- p.adjust(metacoder$data$diff_table$wilcox_p_value,
                                                 method = "fdr")
```
```{r Plot Heat Trees}
metacoder.plot <- heat_tree_matrix(metacoder,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-5, 5),
                 edge_color_interval = c(-5, 5),
                 key_size = 0.75,
                 overlap_avoidance = 5,
                 layout = "fruchterman-reingold", # The primary layout algorithm
                 initial_layout = "reingold-tilford", # The layout algorithm that initializes node locations
                 node_size_axis_label = "Number of ASVs",
                 node_color_axis_label = "Log2 ratio median proportions") 
metacoder.plot
```

![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/16S/data/bin/16S_Metacoder.png)




###### Supplementary

group_pal <- c("#e3aa74", "#ed828c", "#7bb6bd")
rich_samples <- subset_samples(physeq, Sample_Type=="Sample")
rich_samples <- prune_taxa(taxa_sums(rich_samples) > 0, rich_samples)

plot_richness(rich_samples, x="Feed_Type", measures=c("Chao1", "Shannon"))

hill_data <- rich_samples@otu_table@.Data
hill_data <- as.matrix(hill_data)

hill_0 <- hill_div(rich_samples@otu_table,0)
hill_1 <- hill_div(rich_samples@otu_table,1)
hill_2 <- hill_div(rich_samples@otu_table,2)

summary(hill_0)
mean(hill_0)
sd(hill_0)


hierarchy <- md[md$Sample_Type=="Sample",]
hierarchy <- hierarchy[,c(1,6)]

div_profile(hill_data)
div <- div_profile(hill_data, hierarchy=hierarchy, level="alpha")

div_profiles <- div_profile_plot(div, colour = group_pal)

#With post-hoc analyses
divtest <- div_test(hill_data,qvalue=2,hierarchy=hierarchy,posthoc=TRUE)
divtestdata <- divtest$data
divtestdata$Feed <- as.factor(divtestdata$Group)
divtestdata$Feed <- factor(divtestdata$Feed, levels = as.character(unique(divtestdata$Feed)))

stat.test <- divtestdata %>%
  group_by("Feed") %>%
  tukey_hsd(Value ~ Feed) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  add_x_position(x = "Feed", dodge = 0.8) %>%
  add_y_position()

    # Create a box plot
plot = ggboxplot(
      divtestdata, x = "Feed", y = "Value", 
      color = "black",
      fill = "Feed", palette = group_pal,
      outlier.shape = 8, order = c("Control", "Probiotics", "Synbiotics"),
      size = 0.5,
      title = "")  + 
      stat_pvalue_manual(
        stat.test,  label = "{p.adj.signif}", tip.length = 0.045,
        step.increase = 0.09,
        position = "identity", 
        y.position = 12) 
plot = plot + xlab("Feeding Type") + #changing labels
  ylab("Effective number of ASVs")


cowplot::plot_grid(div_profiles,plot, labels = 'AUTO')

![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/16S/data/bin/Diversity_analysis.jpg)
