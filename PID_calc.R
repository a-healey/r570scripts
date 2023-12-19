#!/usr/bin/env Rscript
options(warn=-1)

################################################################################
# -- Load R's argument parser
if(!requireNamespace("argparser",quietly = T)){
  cat("Installing the argparser package\n")
  install.packages("argparser", repos='http://cran.us.r-project.org')
}
if(!requireNamespace("Biostrings",quietly = T)){
  cat("Installing the Biostings package\n")
  install.packages("Biostrings", repos='http://cran.us.r-project.org')
}
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("Biostrings"))
################################################################################
# -- Set up the command line arguments
p <- arg_parser(description = "Align and calculate PID score for protein pairs")
p <- add_argument(p, arg = "file1", help="protein1 to align")
argv <- parse_args(p)
#print(argv)
f1 <- parse_args(p)$file1


######################
seqs <- read.table(f1, quote="\"", comment.char="")

x1 <- AAString(seqs$V3)
x2 <- AAString(seqs$V5)

palign1 <- pairwiseAlignment(x1,x2)
#palign1
match<-nmatch(palign1)
mismatch<-nmismatch(palign1)
alignLen<-nchar(palign1)
PID<-pid(palign1, type = "PID2")

out<-data.frame(key= seqs$V1, seq1=seqs$V2, seq2=seqs$V4,Match=match,Mismatch=mismatch,AlignLength=alignLen, ID=PID)
kaksf <- paste0("pair_",seqs$V1,"_",seqs$V2,"_",seqs$V4,".protAlign")
write.table(out, file = kaksf, sep = "\t", row.names = F, quote = F,col.names = F)