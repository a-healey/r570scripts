rm(list = ls())
library(ggplot2)
library(data.table)
library(GENESPACE)

#the purpose of this script is to read through a dataframe and convert values to run-length equivalents (RLEs)
#RLEs are designed to find the start and stop regions of features in a genome,
#while ignoring outlier data points

#for example, if you wanted to determine a regions of the genome where read coverage is greater 
#than 30:
'''
df=
position,coverage
1,56
2,45
3,55
4,46
6,23
12,34
19,35
22,19
24,7
28,15

the RLE script with a run lenght of 3 will tell you that coverage > 30 from 1-19
position 6 with a coverage of 3 will be ignored because there are not 3 consecutive positions 
where coverage is below 30.


'''
#

# -- read in the data
#in this example, we are examining and obed file from the NucFreq analysis pipeline
d <- fread("nucFreqOutput.obed")

#name the columns in the dataframe
colnames(d) <- c("Chrom","start","end","primaryCov","AltCov")

# -- order values by position
setkey(d, Chrom, start, end)

#add true/false column the run-length equivalent operates on..
d[,tf := primaryCov > 30]
d[,rle := add_rle(tf, which = "n"), by = "Chrom"]

#next, determine the RLE on the true/false column
RLElength = 10
Out <- subset(d, rle >=RLElength)

#next, add blockIds to the dataframe to determine start/stop of each RLE
Out[,blkid := add_rle(tf, which = "id"), by = "Chrom"]


#calculate start/stop for each RLE block and sort
dscc <- Out[,list(start = min(start),
                  end = max(end)), by = c("Chrom", "blkid","tf")]
dscc[,size := end - start]
setkey(dscc, Chrom, blkid,tf)

#write the output
fwrite(
  dscc,
  file = "nucFreqOutput.obed.RLE10.dat",
  sep = "\t")


###########################