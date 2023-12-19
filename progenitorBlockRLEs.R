library(data.table)
library(Biostrings)
library(GENESPACE)
library(ggplot2)

#this is a superful script that can be generalized to quick search exact kmers in a genome
#in this example, we use to assign sequences to progenitor blocks, but this 
#can easily be addapted to quickly perform exact kmer matches across the genome

#script runs in 2 parts:
#first, exact kmer matches are searched in a reference genome

#second, the kmer matches are converted into run-length equivalents (RLE).  In this case, 
#it is used to assign blocks to one progenitor or the other (based on the kmers)


################Part1#############################
#read in the kmer fasta file and reverse complement
kmerFa <- readDNAStringSet("combKmer.fa")
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
fwrite(kmerLoc, file = "subsetKmerLocs.csv")

#######################

######################Part2##################################

#read in teh kmer location file
kl <- fread("subsetKmerLocs.csv", key = c("chr", "start", "end"))

#in this instance, progenitor origin (So- Saccharum officinarum; Ss- Saccharum spontaneum)
#are detected from each kmer
kl[,hap := ifelse(grepl("So", kmerID), "So", "Ss")]


#function to calculate RLEs
add_rle <- function(x, which = "n"){
  if (which == "n") {
    rep(rle(x)$lengths, rle(x)$lengths)
  }else{
    rep(1:length(rle(x)$lengths), rle(x)$lengths)
  }
}
#in order for a block to be assigned to one progenitor or another,
#the RLE needs to detect 100-consecutive calls of the same progenitor with a break length of 10 calls 

kl[,rl := add_rle(hap), by = "chr"]
for(i in seq(from = 10, to = 100, by = 10)){
  kl <- subset(kl, rl > i)
  kl[,rl := add_rle(hap), by = "chr"]
  print(nrow(kl))
}

hapCoords <- kl[,list(start = min(start), end = max(end)),
                by = c("chr", "hap", "rl")]
setkey(hapCoords, chr, start)

fwrite(hapCoords, file = "progenitorRLEs.min100kmerRun.csv")

