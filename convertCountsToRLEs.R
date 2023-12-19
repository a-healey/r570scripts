rm(list = ls())
library(ggplot2)
library(data.table)

#this script converts minimap2 PAF read counts to haplotype collapsed run-length equivalents
#uses coverage to determine which regions contain zero unique mapped reads, and 1-4X copy regions

#this script operates on the output of combinePAFsAndCount.R

# -- read in the data
d <- fread("combinedPAF_counts.dat")

# -- get the chromosome number by striping text)
d[,chrn := as.numeric(gsub("[^0-9.-]", "", chr))]

# -- assign sequences to primary or alternate based on their names
#in R570, alternate sequences either have 'alt' in their names, or
#their scaffold number is greater than 144
#this is pretty specific to R570
d[,hap := ifelse(grepl("Chr", chr) & !grepl("alt", chr), "primary",
                 ifelse(grepl("alt", chr) | chrn > 144, "alt", "primary"))]

# -- calculate median depth and coverage relative to median
d[,`:=`(depth = cov, cov = cov/median(cov))]

# -- drop any high coverage regions where coverage is 5x times greater than the median depth
ds <- subset(d, cov >= 0 & cov < 5)

# -- get list of chromosomes for plotting
chrs <- unique(ds$chr[!grepl("scaff", ds$chr)])

# -- plot coverage peaks and classify into coverage bins
#best way to this is to plot the coverage, apply the thresh ablines to the figure, then refine..
thresh <- c(-Inf, .25, 1.4, 2.3, 3.5, Inf)
ds[,class := as.numeric(as.character(cut(cov, thresh, c(.5, 1:4))))]

ds[,lcov := ifelse(cov < 4.5, cov, 4.5)]
hist(ds$cov, breaks = 1000)
abline(v = thresh, col = "red", lty = 2)

with(subset(ds, cov > 0 & cov < 6), hist(cov, breaks = 1000))
abline(v = thresh, col = "red", lty = 2)

# -- order values by position
setkey(ds, hap, chr, start, end)

# -- calculate run lengths of strings of consecutive identical classes
ds[,runLen := GENESPACE::add_rle(class, which = "n"), by = "chr"]

# -- iteratively drop runs to a min run of 10
dsc <- subset(ds, runLen > 1)
for(i in 1:10){
  dsc <- subset(dsc, runLen > i)
  dsc[,runLen := GENESPACE::add_rle(class, which = "n"), by = "chr"]
}

# -- name the runs
dsc[,run := GENESPACE::add_rle(class, which = "id"), by = "chr"]

# -- get run coordinates
dscc <- dsc[,list(start = min(start),
                  end = max(end)), by = c("chr", "run", "class", "hap")]
dscc[,size := end - start]
setkey(dscc, hap, class)

# -- get total sequence in each category
dscc[,list(nbpInRun = sum((end - start) * class)), by = "hap"]
setkey(dscc, chr, start, end)
dscc[,jump := c(NA, start[-1]-end[-.N]), by = "chr"]
dscc[,list(nbpInGap = sum(jump * class, na.rm = T)), by = "hap"]


with(subset(ds, hap == "primary"), table(class) * 0:4)

# -- output coordinates of the coverage regions
fwrite(
  dscc,
  file = "ccsCoverage_10kbStep_manualClass4_blkCoords_with0_10windowRuns.txt.gz",
  sep = "\t")


