__author__="Adam Healey, ahealey@hudsonalpha.org"
__date__="Created: 12/14/23"
import gzip
from sys import argv
from Bio import SeqIO
#==============================================================
def real_main():
    script, inputFasta, merSize, stepSize = argv
    merSize = int(merSize)
    stepSize = int(stepSize)
    countDict = {}
    nonUniqSet =set()

    
    for record in SeqIO.parse(inputFasta, "fasta"):
        startCounter = 0
        endCounter = startCounter + merSize
        #get the length of the sequence to calculate stop point
        seqLen = len(record.seq)
        #create list of all possible kmers from the scaffold
        result = [record.seq[i:i+merSize] for i in range(0,len(record.seq),stepSize)]
        #counter = 1
        for kmer in result:
            if len(kmer) >= merSize:
                positionKey = "|".join((record.id,str(startCounter),str(endCounter)))
                startCounter += stepSize
                endCounter = startCounter + merSize
                #check first to see if that sequence has already been viewed before
                if kmer in nonUniqSet: continue
                #next, try and delete the key from the dictionary
                try:
                    del countDict[kmer]
                    nonUniqSet.add(kmer)
                except KeyError:
                    countDict[kmer] = [positionKey]
                    #####
                #####
            #####
        #####
    #####
	#####

    for k,v in countDict.items():
        chrom,start, end = v[0].split("|")
        print(">%s:%s-%s\n%s" % (chrom,start,end,k))
        
#==============================================================
if ( __name__ == '__main__' ):
    real_main()

