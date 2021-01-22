# Performance Analysis


Load dependencies
```
library(readxl)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(dplyr)
library(ggpubr)
library(rstatix)
```


Load data 
```
WGG <- read_excel("performance_data_R_Trial4.xlsx")
wgg <- stack(WGG)
wgg$Feed <- c(rep("PRO",5),rep("CTRL",5), rep("SYN",5))
wgg <- wgg %>% arrange(Feed)

FCR <- read_excel("performance_data_R_Trial4.xlsx",sheet = "FCR")
fcr <- stack(FCR)
fcr$Feed <- c(rep("PRO",5),rep("CTRL",5), rep("SYN",5))
fcr <- fcr %>% arrange(Feed)

FER <- read_excel("performance_data_R_Trial4.xlsx",sheet = "FER")
fer <- stack(FER)
fer$Feed <- c(rep("PRO",5),rep("CTRL",5), rep("SYN",5))
fer <- fer %>% arrange(Feed)

SGR <- read_excel("performance_data_R_Trial4.xlsx",sheet = "SGR")
sgr <- stack(SGR)
sgr$Feed <- c(rep("PRO",5),rep("CTRL",5), rep("SYN",5))
sgr <- sgr %>% arrange(Feed)

SFR <- read_excel("performance_data_R_Trial4.xlsx",sheet = "SFR")
sfr <- stack(SFR)
sfr$Feed <- c(rep("PRO",5),rep("CTRL",5), rep("SYN",5))
sfr <- sfr %>% arrange(Feed)

LER <- read_excel("performance_data_R_Trial4.xlsx",sheet = "LER")
ler <- stack(LER)
ler$Feed <- c(rep("PRO",5),rep("CTRL",5), rep("SYN",5))
ler <- ler %>% arrange(Feed)

PER <- read_excel("performance_data_R_Trial4.xlsx",sheet = "PER")
per <- stack(PER)
per$Feed <- c(rep("PRO",5),rep("CTRL",5), rep("SYN",5))
per <- per %>% arrange(Feed)

## Group colours!
group_pal <- c("#e3aa74", "#ed828c", "#7bb6bd")
```

Pairwise Tukey HSD and boxplots for aquaculture related abilities
```
stat.test <- wgg %>%
  group_by("Feed") %>%
  tukey_hsd(values ~ Feed) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

wgg.stats <- stat.test
wgg.stats$Data <- c(rep("wgg",3))

stat.test <- stat.test %>%
  add_x_position(x = "Feed", dodge = 0.8) %>%
  add_y_position()
    # Create a box plot

WGG.plot = ggboxplot(
      wgg, x = "Feed", y = "values", 
      color = "black",
      fill = "Feed", palette = group_pal,
      outlier.shape = 8,
      #order = c("epsilon", "beta", "zeta"),
      size = 0.5,
      title = "")  + 
      stat_pvalue_manual(
        stat.test,  label = "{p.adj} {p.adj.signif}", tip.length = 0.02,
        step.increase = 0.075) 
WGG = WGG.plot + xlab("Feeding Type") + #changing labels
  ylab("Weight Gain (%)")

stat.test <- fcr %>%
  group_by("Feed") %>%
  tukey_hsd(values ~ Feed) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

fcr.stats <- stat.test
fcr.stats$Data <- c(rep("fcr",3))

stat.test <- stat.test %>%
  add_x_position(x = "Feed", dodge = 0.8) %>%
  add_y_position()
    # Create a box plot
FCR.plot = ggboxplot(
      fcr, x = "Feed", y = "values", 
      color = "black",
      fill = "Feed", palette = group_pal,
      outlier.shape = 8,
      #order = c("epsilon", "beta", "zeta"),
      size = 0.5,
      title = "")  + 
      stat_pvalue_manual(
        stat.test,  label = "{p.adj} {p.adj.signif}", tip.length = 0.2,
        step.increase = 0.075) 
FCR = FCR.plot + xlab("Feeding Type") + #changing labels
  ylab("Feed Conversion Ratio")

stat.test <- fer %>%
  group_by("Feed") %>%
  tukey_hsd(values ~ Feed) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

fer.stats <- stat.test
fer.stats$Data <- c(rep("fer",3))

stat.test <- stat.test %>%
  add_x_position(x = "Feed", dodge = 0.8) %>%
  add_y_position()

FER.plot = ggboxplot(
      fer, x = "Feed", y = "values", 
      color = "black",
      fill = "Feed", palette = group_pal,
      outlier.shape = 8,
      #order = c("epsilon", "beta", "zeta"),
      size = 0.5,
      title = "")  + 
      stat_pvalue_manual(
        stat.test,  label = "{p.adj} {p.adj.signif}", tip.length = 0.2,
        step.increase = 0.075) 
FER = FER.plot + xlab("Feeding Type") + #changing labels
  ylab("Feed Efficiency Ratio (%)")

stat.test <- ler %>%
  group_by("Feed") %>%
  tukey_hsd(values ~ Feed) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

ler.stats <- stat.test
ler.stats$Data <- c(rep("ler",3))

stat.test <- stat.test %>%
  add_x_position(x = "Feed", dodge = 0.8) %>%
  add_y_position()
    # Create a box plot
LER.plot = ggboxplot(
      ler, x = "Feed", y = "values", 
      color = "black",
      fill = "Feed", palette = group_pal,
      outlier.shape = 8,
     # order = c("epsilon", "beta", "zeta"),
      size = 0.5,
      title = "")  + 
      stat_pvalue_manual(
        stat.test,  label = "{p.adj} {p.adj.signif}", tip.length = 0.2,
        step.increase = 0.075) 
LER = LER.plot + xlab("Feeding Type") + #changing labels
  ylab("Lipid Efficiency Ratio")

stat.test <- per %>%
  group_by("Feed") %>%
  tukey_hsd(values ~ Feed) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

per.stats <- stat.test
per.stats$Data <- c(rep("per",3))

stat.test <- stat.test %>%
  add_x_position(x = "Feed", dodge = 0.8) %>%
  add_y_position()
    # Create a box plot
PER.plot = ggboxplot(
      per, x = "Feed", y = "values", 
      color = "black",
      fill = "Feed", palette = group_pal,
      outlier.shape = 8,
     # order = c("epsilon", "beta", "zeta"),
      size = 0.5,
      title = "")  + 
      stat_pvalue_manual(
        stat.test,  label = "{p.adj} {p.adj.signif}", tip.length = 0.2,
        step.increase = 0.075) 
PER = PER.plot + xlab("Feeding Type") + #changing labels
  ylab("Protein Efficiency Ratio")
```

Merge stats and make plots
```
stats <- rbind(wgg.stats,fcr.stats,fer.stats,ler.stats,per.stats)
#write.csv(stats, "Feed_Stat_summary.csv")

#pdf("Performance_Boxplots.pdf", height = 8, width = 8)
plot_grid(WGG,FER, LER, PER)
#dev.off()
```

Correlation between aquacultural abilities, like Feed conveersion ratio, growth, protein efficiency ratio, and lipid efficiency ratio
```
corr_data <- as.data.frame(wgg$values)
corr_data$PER <- per$values
corr_data$LER <- ler$values
corr_data$FCR <- fer$values
colnames(corr_data) <- c("WGG", "PER", "LER", "FER")
cor(wgg$values, per$values, method = c("pearson", "kendall", "spearman"))
cor.test(wgg$values, per$values, method=c("pearson", "kendall", "spearman"))
color <- c("#e3aa74", "#ed828c", "#7bb6bd")
library("ggpubr")
WGG_PER <- ggscatter(corr_data, x = "PER", y = "WGG", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Protein Efficiency Ratio", ylab = "Gain of Weight")

WGG_LER <- ggscatter(corr_data, x = "LER", y = "WGG", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "spearman",
                     xlab = "Lipid Efficiency Ratio", ylab = "Gain of Weight",
                     pallette = color)

FER_PER <- ggscatter(corr_data, x = "PER", y = "FER", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "spearman",
                     xlab = "Protein Efficiency Ratio", ylab = "Feed Conversion Ratio")

FER_LER <- ggscatter(corr_data, x = "LER", y = "FER", 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "spearman",
                     xlab = "Lipid Efficiency Ratio", ylab = "Feed Conversion Ratio",
                     pallette = color)


#pdf("Performance_Correlations.pdf", height = 8, width = 8)
plot_grid(WGG_PER, FER_PER, WGG_LER,FER_LER, labels = "AUTO", nrow = 2)
#dev.off()
```

Create a radar plot to sum up data
```
library(fmsb)
RADAR <- read_excel("performance_data_R_Trial4.xlsx", 
                  sheet = "radar chart")
name <- RADAR$feed.code
RADAR <- RADAR[,-c(1)]
Max <- c(5.1,2,640,1.25,1,5,6)
Min <- c(4.2,1.5,600,1.1,0.5,4.8,5)
RADAR <- rbind(Max,Min, RADAR)
RADAR <- RADAR[,-c(4)]
row.names(RADAR) <- c("Max","Min","Beta","Epsilon","Zeta")
colnames(RADAR) <- c("Lipid Efficiency Ratio","Protein Efficiency Ratio","Gain of Weight", "Feed Efficiency Ratio", "Specific Growth Rate", "Speficic FCR")

colors_border=c(rgb(0.2,0.5,0.5,0.8), rgb(0.7,0.5,0.1,0.8),rgb(0.8,0.2,0.2,0.8) )
colors_in=c( rgb(0.2,0.5,0.5,0.2), rgb(0.7,0.5,0.1,0.2),rgb(0.8,0.2,0.2,0.2) )

#pdf(file = "Radar_Feed_performance.pdf")
radar.plot <- radarchart( RADAR  , axistype=2 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", cglwd=0.8,
            #custom labels
            vlcex=0.8)
radar.plot + legend(x=-2.2, y=.4, legend = name, bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.2, pt.cex=3)# Add a legend
#dev.off()
```
