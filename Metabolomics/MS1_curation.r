df <- read.csv("MS1_table.csv")
rt <- df[,1:3]
rownames(df) <- df$name
df <- df[4:63]
df <- df[apply(df, 1, FUN = function(x){sum(x == 0)}) < 30,]

rt <- rt[match(rownames(df), rt$name),]
rownames(rt) <- rt$name
data <- cbind(rt,df)
write.csv(data,"MS1_curated.csv", row.names = FALSE)
