library(data.table)
library(Biostrings)
library(GENESPACE)

#this is a superful script that can finds search exact kmers in a genome

#read in the kmer fasta file and reverse complement
kmerFa <- readDNAStringSet("kmers.fasta")
kmerFar <- reverseComplement(kmerFa)
names(kmerFar) <- sprintf("%s_rev", names(kmerFar))
kmerv <- as.character(c(kmerFa, kmerFar))
kmerv <- kmerv[!duplicated(kmerv)]
di <- names(kmerv); names(di) <- kmerv

#functions to search for and keep only unique kmer matches
find_manyKmers <- function(dnass, kmers, revComp = T, verbose = T){
  # -- convert kmer strings into more accessible format
  kmerList <- lapply(1:length(kmers), function(i) DNAString(kmers[i]))
  # -- optionally add reverse complement seq
  if(revComp){
    kmerListr <- lapply(1:length(kmers), function(i)
      reverseComplement(DNAString(kmers[i])))
    kmerList <- c(kmerList, kmerListr)
  }
  # -- only keep unique kmers
  ndup <- which(!duplicated(sapply(kmerList, as.character)))
  kmerSS <- DNAStringSet(kmerList[ndup])
  # -- make the preprocessed pattern dictionary
  pdict0 <- PDict(kmerSS)
  pmn <- as.character(kmerSS)
  outa <- rbindlist(lapply(names(dnass), function(j){
    if(verbose)
      cat(sprintf("\tSequence %s: ", j))
    ss1 <- dnass[[j]]
    m <- matchPDict(pdict0, ss1)
    whm <- which(sapply(m, length) > 0)
    if(length(whm) > 0){
      out <- rbindlist(lapply(whm,function(i)
        data.table(
          chr = j,
          kmer = pmn[[i]],
          start = start(m[[i]]),
          end = end(m[[i]]))))
      if(verbose)
        cat(sprintf("found %s hits of %s unique kmers\n",
                    nrow(out), uniqueN(out$kmer)))
      return(out)
    }else{
      if(verbose)
        cat("found no exact matches\n")
    }
  }))
  return(outa)
}

#read in the reference genome that you would like to search
fa <- readDNAStringSet("Saccharum_officinarum_X_spontaneum_var_R570.mainGenome.fasta")

#find regions of exact match between kmers of interest and the reference genome
kmerLoc <- find_manyKmers(dnass = fa, kmers = kmerv)
kmerLoc[,kmerID := di[kmer]]

#print the locations of the kmers of interest
fwrite(kmerLoc, file = "kmerLocations.csv")