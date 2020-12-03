setwd("/Users/jacob/Dropbox/Arbejde/PhD/HappyFish/Multi Omic Paper/Metagenomics/Functional stuff")

library("dplyr")
library("readxl")


Coverage <- read.table("Gene_Coverage.txt-GENE-COVERAGES.txt", sep = "\t", header = TRUE)
COG_cat <- read.table("COG_Cats.txt", sep = "\t", header = TRUE)
COG_Func <- read_excel("COG_Function.xlsx")
KEGG <- read_excel("KEGG_Function.xlsx")
Pfam <- read_excel("Pfam_Function.xlsx")

Coverage <- Coverage[order(Coverage$key),]
data <- merge(Coverage,COG_cat, by.x="key", by.y="gene_callers_id", all.x = TRUE)
data <- merge(data,COG_Func, by.x="key", by.y="gene_callers_id", all.x = TRUE)
colnames(data) <- c("key","CTRL_D01","CTRL_D02","CTRL_M01","CTRL_M02","PRO_D01","PRO_D02","PRO_M01","PRO_M02",
                    "SYN_D01","SYN_D02","SYN_M01","SYN_M02","source.COG_CATEGORY","accession.COG","COG_CATEGORY",
                    "e_value.COG_CAT","source.COG_FUNCTION","accession.COG_FUNCTION","COG_FUNCTION","e_value.COG_FUNCTION")
colnames(KEGG) <- c("gene_callers_id","source.KEGG","accession.KEGG","KEGG","e_value.KEGG")
data <- merge(data,KEGG, by.x="key", by.y="gene_callers_id", all.x = TRUE)
Pfam <- Pfam[order(Pfam$gene_callers_id),]
Pfam <- Pfam[!duplicated(Pfam$gene_callers_id),]
colnames(Pfam) <- c("gene_callers_id","source.Pfam","accession.Pfam","Pfam","e_value.Pfam")
Pfam <- na.omit(Pfam)
data <- merge(data,Pfam, by.x="key", by.y="gene_callers_id", all.x = TRUE)
write.csv(data,file = "Merged_Functions.csv")
