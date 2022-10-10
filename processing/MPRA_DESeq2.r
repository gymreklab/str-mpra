#!/usr/bin/Rscript --vanilla
.Library
.Library.site
.libPaths()

# load packages 
if(!require("DESeq2")){
    install.packages("DESeq2")
    library("DESeq2")
}

# get arguments 
args = commandArgs(trailingOnly=TRUE)
matrix <- args[1]
out_dir <- args[2]
group_num <- strtoi(args[3])

# data
data <- read.csv(matrix,
                 row.name="STR")

# column re-ordering
data_order <- vector("integer", group_num * 2)
for (x in 1:group_num){
    data_order[x] <- (2*x - 1)
    data_order[x + group_num] <- (2*x)
}

data <- data[, data_order]

cts <- as.matrix(data)
conditions <- c(rep("gDNA", group_num),
                rep("cDNA", group_num))
samples <- data.frame("run" = c(rep("gDNA", group_num),
                                rep("cDNA", group_num)),
                      "condition"=conditions)

dds <- DESeqDataSetFromMatrix(countData = cts, 
                              colData = samples, 
                              design = ~condition)

# perform median of ratios normalization
# generate size factors 
dds <- estimateSizeFactors(dds)

# retrieve normalized count matrix 
normalized_counts <- counts(dds, normalized=TRUE)

# perform DESeq2
dds <- DESeq(dds)
res <- results(dds)

# write output 
write.csv(as.data.frame(normalized_counts),
          file=paste(out_dir, "normalized_aggregate_count_matrix.csv", sep=""))
write.csv(as.data.frame(res),
          file=paste(out_dir, "deseq2_result.csv", sep=""))