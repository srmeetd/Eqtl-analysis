
#! /usr/bin/Rscript

library(readr)
library(tidyr)
library(dplyr)
library(tximportData)
library(tximport)
library(GenomicFeatures)

#####################################

args <- commandArgs(trailingOnly = TRUE)
TxDb <- makeTxDbFromGFF(paste0(args[1]))
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(TxDb, k, "GENEID", "TXNAME")
workingdir = "quantification.dir/"
file = list.files(path =workingdir, pattern=".sf$")
Files=unlist(file)
name = gsub("trimmed-|.sf", "",Files)
names(Files) <- name
txi <- tximport(paste0(workingdir,Files), type = "salmon", tx2gene = tx2gene,countsFromAbundance = "lengthScaledTPM")
Expression = txi$counts
colnames(Expression) = name
Expression = as.data.frame(Expression)
Expression$ID = rownames(Expression)
EXPRESSION = Expression[moveme(names(Expression), "ID first")]
rownames(EXPRESSION) = NULL
write.table(EXPRESSION, file = paste0(args[3],"_quntified.tsv"), sep="\t", row.names=TRUE, quote=FALSE)
