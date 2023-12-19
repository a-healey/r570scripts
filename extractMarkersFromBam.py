
from string import maketrans

from os.path import split

from sys import stdin, stdout

from sys import argv

import subprocess
from collections import Counter
#==============================================================
def real_main():
    script, bedFile, bam_File, merSize = argv
    
    #Setup
    merLength = int(merSize)
    bedDict = {}
    
    #Read in the positions from the bedFile
    for line in open(bedFile):
        s 		= line.split(None)
        chrom 	= s[0]
        start 	= int(s[1])
        end 	= int(s[2])
        key = '_'.join((chrom, str(start)))
        try:
            bedDict[key].append((chrom, start, end))
        except KeyError:
            bedDict[key] = [(chrom, start, end)]
            #####
        #####
	#####
    #Extract positions bed dictionary
    for k, v in bedDict.items():
        kmerDict  	= {}
        countDict 	= {}
        key         = k
        chrom 		= v[0][0]
        windowStart = v[0][1]
        windowEnd 	= v[0][2]

        cmd = 'samtools view %s -q 50 "%s:%d-%d" | cut -f4,10' % (bam_File, chrom, windowStart, windowEnd)
        #print cmd
        p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        for line in p.stdout:
            x               = line.split(None)
            readStart     	= int(x[0])
            read        	= x[1]
            newRead   = read[(windowStart - readStart):]
            if len(newRead) >= merLength:
                mer = newRead[:(merLength)]
                try:
                    countDict[mer] += 1
                except KeyError:
                    countDict[mer] = 1

        p.poll()
        
        #only consider regions with a read depth greater than 4
        #can be adjusted as needed
        newDict = dict((k,v) for k,v in countDict.items() if v > 4)
        
        #use a decorated sort to find the most frequent kmer at that position
        DSU = [ (v,k) for k, v in newDict.items()]
        DSU.sort(reverse=True)
        depth = float(sum(newDict.values()))
        
        #if there is only 1 kmer present and it is counted more than once
        if ((len(DSU) == 1) and (depth > 1)):
            Hap1 = DSU[0][1]
            print("%s_%d\t%s\tNONE\t%f\t0\t100\t100\t0\t%f" % (chrom, windowStart, Hap1, depth, depth))
        
        # if there are multiple kmers present at a locus
        if len(DSU) > 1:
            #assign kmers to het1 (most frequent) and het2 (second most frequent)
            #could be modified to capture more kmers if needed
            #eg. het3 = DSU[2][1]
            het1 = DSU[0][1]
            het2 = DSU[1][1]
            #get the frequency counts for each kmer
            het1_count  = float(DSU[0][0])
            het2_count  = float(DSU[1][0])
            #calculate how many kmers total are present at this locus, and what proportion is represented by het1 and 2
            summed      = float(het1_count + het2_count)
            rel_abund   = (summed / depth) * 100
            fract1      = float(het1_count / depth) * 100
            fract2      = float(het2_count/ depth) * 100
            het_fract1  = float(het1_count / summed) * 100
            het_fract2  = float(het2_count / summed) * 100
            print("%s_%d\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f" % (chrom, windowStart, het1, het2, het1_count, het2_count, rel_abund, het_fract1, het_fract2, depth))

#==============================================================
if ( __name__ == '__main__' ):
    real_main()