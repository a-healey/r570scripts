library(data.table)

#the purpose of this script is to combine multiple minimap alignment files (paf format)
#and calculate the average depth over a specified window length

#this is useful because of how messy read alignments can get in some regions when 
#aligning reads in complex polyploids

#function to count reads in paf alignment files
#pafFile = minimap alignment file (paf format)
#faiFile = samtools faidx index file for fasta file of interest
#stepSize= window size to calculate average depth call
#minAlignmentSize= minimum length of alignment to consider
#minPercId= minimum percent identity between read and reference to consider
#minMapQ= minimum mapping quality score between read and reference to consider

#also, these are default parameters, but can be changed when calling the function, see below
count_readsPAF <- function(pafFiles,
                           faiFile,
                           stepSize = 1e4,
                           minAlignSize = 5e3,
                           minPercId = 90,
                           minMapQ = 50){
  # -- make the windows
  fai <- fread(faiFile, select = c(1, 2), col.names = c("chr", "end"))
  winds <- fai[,list(start = seq(from = 1, to = end, by = stepSize)), by = "chr"]
  winds[, end := start]
  setkey(winds, chr, start, end)
  # -- read in the paf files, parse to best alignments
  dl <- rbindlist(lapply(pafFiles, function(i){
    x <- fread(
      i, select = c(3:4, 6, 8:9, 10, 12),
      col.names = c("readStart", "readEnd",
                    "chr", "start", "end", "nMatches", "mapQ"),
      showProgress = F,
      key = c("chr", "start", "end"))
    x <- subset(x, mapQ >= minMapQ & mapQ != 255)
    x[,mapQ := NULL]
    x[,alignSize := readEnd - readStart]
    x[,percID := 100 * (nMatches/alignSize)]
    x <- subset(x, alignSize >= minAlignSize & percID >= minPercId)
    fo <- foverlaps(x, winds)[,list(n = .N), by = c("chr", "start", "end")]
    return(fo)
  }))
  cnts <- dl[,list(cov = sum(n)), by = c("chr","start", "end")]
  out <- merge(winds, cnts, by = c("chr","start", "end"), all.x = T)
  out[is.na(out)] <- 0
  setkey(out, chr, start, end)
  return(out)
}

faiFile <- "reference.fasta.fai"

pafFiles <- "alignment.paf.bz2"

#or vectorize the libs to run all at once from a directory
pafFiles <- list.files(path = "/path/to/directory/containing/paf/files/", pattern = ".paf", full.names = T)

#generate a count file, but use zero for minAlignSize and minPercId
countFile <- count_readsPAF(faiFile = faiFile, pafFiles = pafFiles, minAlignSize = 0, minPercId = 0)

#write to a file
fwrite(countFile,file='combinedPAF_counts.dat', quote =F)
