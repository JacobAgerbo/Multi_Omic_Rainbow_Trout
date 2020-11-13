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

```{r load dependencies}
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
