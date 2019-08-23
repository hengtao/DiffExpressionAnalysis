library(edgeR)
library(airway)
args = =commandArgs(T) 
data <- args[1]
outfile <- args[2]
dispersion <- args[3]

G1_vs_G2 <- read.table(data, header = T, sep = "\t", row.names = 1)

group <- 1:2
y <- DGEList(counts=counts, group = group)
y <- DGEList(counts=G1_vs_G2, group = group)

keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]

y <- calcNormFactors(y)
y_bcv <- y
bcv <- dispersion
et <- exactTest(y_bcv, dispersion = bcv ^ 2)

write.table(et$table,paste(outfile, ".Diff.edgeR.txt", sep = ""), sep = "\t")