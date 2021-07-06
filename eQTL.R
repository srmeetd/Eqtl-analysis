#! /usr/bin/Rscript

library(MatrixEQTL) 
library (data.table)

args <- commandArgs(trailingOnly = TRUE)
SNP_data = fread (args[1], header  = FALSE)
snp_location = fread (args[2], header = TRUE)
colnames(snp_location) = c("snp","chr","pos")

SNP_file_name = args[1]
snps = SlicedData$new()
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "./."; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = length(rownames(SNP_data));     # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

#################################
expression_data = fread(args[3],header = FALSE)
expression_file_name = args[4]
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = length(rownames(expression_data));      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

####################################
cvrts  = fread(args[6],header = FALSE)
cvrt_link = args[6]
cvrt= SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = length(rownames(cvrts));      # read file in slices of 2,000 rows
cvrt$LoadFile(cvrt_link)

#####################################

output_file_name = tempfile()
pvOutputThreshold = 1e-2
errorCovariance = numeric()
useModel = modelLINEAR # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();
pvOutputThreshold_cis = 1e-2;
pvOutputThreshold_tra = 1e-2;
cisDist = 1e6

#######################################
snpspos =  snp_location
genepos  = fread(args[5],header = TRUE)
colnames(genepos) = c("geneid","chr","s1","s2")
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = FALSE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra)
unlink(output_file_name_cis)

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
a = as.data.frame (me$trans$eqtls)
head (a)
write.table (me$trans$eqtls, file = args[7], sep = "\t", quote = FALSE)
write.table (me$cis$eqtls, file = args[8], sep = "\t", quote = FALSE)


##########################################


