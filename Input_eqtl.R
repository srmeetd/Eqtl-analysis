#! /usr/bin/Rscript

library (data.table)

args <- commandArgs(trailingOnly = TRUE)
names = read.table (args[1], header = FALSE)
names1  = gsub ("trimmed-|.star","",names$V1)
SNP = fread(args[2],header =FALSE)
SNP1 =  SNP[!grepl("^.$",SNP$V3),]
SNP = SNP1
snp_location = subset(SNP, select = c("V3","V1","V2"))
colnames(snp_location) = c("snp","chr","pos")
write.table (snp_location,file = args[3],sep = "\t", row.names = FALSE, quote = FALSE)

genotype = SNP
genotype$V1 = NULL
genotype$V2 = NULL
colnames(genotype) = c("id",paste0(names1))
#(genotype) = make.names(genotype$id,unique = TRUE)
genotyp1 = as.data.frame(gsub ("^0/1$|^0/2$|^0/3$|^0/4$","1",as.matrix(genotype)))
genotyp1 = as.data.frame(gsub ("^1/1$|^2/2$|^3/3$|^4/4$","2",as.matrix(genotyp1)))
genotyp1 = as.data.frame(gsub("^0/0$","0",as.matrix(genotyp1)))
#genotyp1 = as.data.frame(gsub ("[./.]","3",as.matrix(genotyp1)))

colnames(genotyp1) = gsub('(.*)_\\w+', '\\1',colnames(genotyp1))
expression = read.table (args[5],sep = "\t",header = TRUE)
g = genotyp1[,colnames(genotyp1)%in%colnames(expression)]

###################################
#genotyp1$ID = rownames(genotyp1)
write.table (g, file = args[4],sep = "\t", quote = FALSE,row.names = FALSE)
