# Intestinal multi omics unravels a differentiated metabolic landscape and loss of important gut bacteria Mycoplasma salmoninae in rainbow trout (_Oncorhynchus mykiss_)  using non marine pre- and probiotics 

_Please find more information of study Rasmussen et al. 2021_

## Abstract
Animal protein production is increasingly looking towards microbiome associated services as to further improve production sustainability. Here, we investigate the functional effects of   pro- and synbiotic feed additives on microbiome associated functions in relation to growth performance in the commercially important rainbow trout (Oncorhynchus mykiss). We screen multiple omics datasets from gut content samples, including 16S bacterial profiling, whole metagenomes, and untargeted metabolomics, to functionally investigate presence of bacteria and their molecular interactions with host metabolism. Our findings reveal __(I)__ fish reared on the control feed outperforms fish fed with either pro- or synbiotics additives in scores of gain of weight, feed conversion ratio, and a significant higher protein efficiency ratio, __(II)__ significant changes of the metagenome across fish fed with different feeding types, where rainbow trout reared on the pro- or synbiotics feeds had a reduced abundance of Mycoplasma salmoninae in the intestine,  __(III)__ inferring functionality of Mycoplasma salmoninae related to feed, using genome resolve metagenomics, __(IV)__ significant changes of the intestinal metabolomic landscape including an increase of lipids, lipid-like metabolites, and amino amides classified metabolites in fish reared on the control feed.

## 16S rRNA gene profiling of bacteria in rainbow trout gut
_Please find scripts for processing 16S data at: https://github.com/JacobAgerbo/Data2Result._

Pipeline is briefly described below:

* Demultiplexing and removal of adaptors and low quality reads were done with AdapterRemoval/v2.2.4, with a quality base of 30 and a minimum length of 50bp. 
* Microbial 16S data were further filtered, trimmed according to error rate, and ASV-clustered analysed through DADA2(Callahan et al., 2016). 
* Taxonomy was assigned through DADA2, using Silva/v138. 
* Post clustering algorithms were applied to minimise false positives, using LULU(Frøslev et al., 2017) and subsequently  contaminations were removed from samples, using decontam(Davis et al., 2018). 
* Composition analysis were carried out using phyloseq (McMurdie & Holmes, 2013) 
* Differential abundance analysis across feeding groups were carried out using metacoder, using wilcoxon rank sum test and FDR correction for multiple comparisons ( Foster et al., 2017).

## Genome resolve metagenomics
Metagenomic pipeline were performed each individual of feeding type inclduing: 

Feed| n |Seq Chemistry | Gut Section|ENA Archive|
--- | --- | --- | ---|---
__Control Feed__|2 |  PE150 | Mid Gut |PRJEB40990
__Control Feed__|2 |  PE150 | Distal Gut |PRJEB40990
__Probiotic Feed__|2 |  PE150 | Mid Gut |PRJEBXXXX
__Probiotic Feed__|2 |  PE150 | Distal Gut |PRJEB40990 & PRJEBXXXX
__Synbiotic Feed__|2 |  PE150 | Mid Gut |PRJEBXXXX
__Synbiotic Feed__|2 |  PE150 | Distal Gut |PRJEBXXXX

Host filtered metagenomic sequences can be found on respective ENA archives.

### Genome resolve metagenomics preprocessing

_Please find information of parameters and modules used to trim and filter data in Metagenomics/Metagenomic_pipeline.sh_

Pipeline is briefly described below:

* Quality Control, using FastQC
* Adaptor Removal, using AdapterRemoval
* Duplicate and singleton removal, using bbmap and re-pair.sh
* Filtering of phiX, human, and rainbow trout, using BWA
* Co-assembly to contigs of minimum 1000 bp, using MegaHIT
* Quality control of contigs, using Quast

### Binning and MAG curation

_Please find information of parameters and modules used in Anvio_pipeline.sh_

After assembly of host filtered reads were contigs processed in Anvio. 
The subsequent workflow is outlined at http://merenlab.org/2016/06/22/anvio-tutorial-v2/. Briefly; 
* anvi’o was used to profile the scaffolds using Prodigal/v2.6.338 with default parameters to identify genes and HMMER/v.3.339 to identify genes matching to archaeal, protista (based on http://merenlab.org/delmont-euk-scgs), and bacterial single-copy core gene collections. Also, ribosomal RNA based HMMs were identified (based on https://github.com/tseemann/barrnap). The HMMs were used to determine completeness of  metagenome assembled genomes (MAGs); 
* Kaiju 41 was used with NCBI’s non-redundant protein database ‘nr’ to infer the taxonomy of genes (as described in http://merenlab.org/2016/06/18/importing-taxonomy/); 
* we mapped short reads from the metagenomic set to the scaffolds using BWA/v0.7.1596 (minimum identity of 95%) and stored the recruited reads as BAM files using samtools 
* anvi'o profiled each BAM file to estimate the coverage and detection statistics of each scaffold, and combined mapping profiles into a merged profile database for each metagenomic set. Contigs were binned automatically, using CONCOCT, by constraining the number of clusters per metagenomic set to 10.
* Bin and MAGs where curated, following guideline from Veronika Kivenson: http://merenlab.org/2017/05/11/anvi-refine-by-veronika/
* Functional infeering were carried using KEGG and RAST annotations

## Metabolomic Analysis

_The metabolomics datasets generated and analysed during the current study is available in the MSV000084364 repository, ftp://massive.ucsd.edu/MSV000084364/.

Data for redoing statistical analysis can be found in metabolomics data folder._

### Statistical Analysis
Pipeline is briefly described below:

* Removal of false positives and minimise zero inflation of dataset
* Cumulative fractional abundance of each metabolite were assesed
* PCoA Ordination of metabolite composition across feeding types
* Differential Intensity test of Metabolites across feeding types
* Clustering and composition analysis of differential abundant metabolites across feeding types and _Mycoplasma_

### Network Analysis to improve metabolite annotation

Pipeline is briefly described below, please refer to main text for GNPS details.

A molecular network was created using the online workflow https://ccms-ucsd.github.io/GNPSDocumentation on the GNPS website http://gnps.ucsd.edu.

In order to enhance identification of unknown metabolites, unsupervised substructures were discovered by combining MS2LDA based on MS2 spectre of all samples(van der Hooft et al., 2016), in silico structures of MS2 spectre annotated using Network Annotation Propagation (NAP)(da Silva et al., 2018).

Peptidic natural products (PNPs) were identified, using DEREPLICATOR(Mohimani et al., 2017) with VARQUEST. 

Subsequently, chemical classifications were assessed, using ClassyFire(Djoumbou Feunang et al., 2016). 

All structural annotations from GNPS were combined, using MolNetEnhancer(Ernst et al., 2019). 

Further annotation of metabolites was carried out, using MetDNA(Shen et al., 2019). Networks were visualised using Cytoscape/v3.8.0. 
