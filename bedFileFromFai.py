
from sys import argv
#==============================================================
def real_main():
    #USEAGE= use samtools in faidx reference then use fai file to generate windows
    #mersize = size of kmer you want from kmer generator
    # stepsize = sliding window size for next kmer
    #if you want non-overlapping windows, set stepsize = mersize+1
    
    script, faiFile, merSize, stepSize = argv
    
    
    merSize = int(merSize)
    stepSize = int(stepSize)
    mainDict = {}
    for line in open(faiFile):
        s = line.split(None)
        chrom = s[0]
        length = int(s[1])
        try:
            mainDict[chrom].append(length)
        except KeyError:
            mainDict[chrom] = [length]
            #####
        #####
    #####
    
    for k, v in mainDict.items():
        counter = 0
        for x in range (0, v[0]):
            start = counter
            end = counter + merSize
            if end <= v[0]:
                print("%s\t%d\t%d" % (k, start, end))
                counter = start + stepSize
            #####
        #####
    #####
  

#==============================================================
if ( __name__ == '__main__' ):
    real_main()