# Analysis of 16S data from intestinal environment of rainbow trout (_Oncorhynchus mykiss_)

### Author: Jacob Agerbo Rasmussen
### Date: 09.nov.2020

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
I use Phyloseq to organise our data for profiling the V3-V4 region of 16S rRNA gene, thank you Joey711!
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
![alt text](https://github.com/JacobAgerbo/Multi_Omic_Rainbow_Trout/blob/main/16S/data/bin/16S_Barplot.pdf)
