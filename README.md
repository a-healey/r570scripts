# r570scripts
This repository contains scripts to build and explore the R570 hybrid sugarcane genome assembly. These scripts are not a pipeline to re-create the assembly, but instead are presented here to complement the genome assembly methods presented in the supplementary material of the associated manuscript. 
**Note** A link to the article and the supplementary information will be provided here upon publication. 

When constructing chromosomes for R570, we found it most useful to use a standard set of scripts and use generalized plots to visualize patterns and interpret.  Therefore, to make this Github repository useful for the widest audience, we have included both:
  
  1) Useful scripts for general analysis, and
  2) Annotated analysis scripts used for R570.

Also, we have included a conda environment yml file (```R570_analysisEnv.yml```) to enable others to run python scripts, as well as other useful softwares frequently used during chromosome construction and analysis (minimap2, mummer, samtools, bedtools, biopython, etc..)

##############################

script: ```extractMarkersFromBam.py```

Explanation:

This is a useful code that generates genetic markers directly from an Illumina reads aligned to a reference genome.  It is especially useful for making markers from diverse genotypes, as the markers are extracted from the reads themselves, not the reference sequence.

The script operates on a bed file, so you can customize marker extraction to regions of interest in the genome.  Also, if you would like all possible markers from the alignment, then you can generate an exhaustive bed file from a samtools .fai index file using bedFileFromFai.py

```

bedFileFromFai.py usage:
  #use samtools faidx to index the assembly
    samtools faidx reference.fasta
  #generate bed file
    python bedFileFromFai.py reference.fasta.fai markerSize stepSize > bedFile.bed
    #markerSize: length of the marker would you like to generate.  (Based on experience, I would not go beyond 80 bp for Illumina reads (otherwise observing the markers in 150 bp reads becomes rare)
    #stepSize: length of basepair walk before generating a new marker.

  #example: for non-overlapping 80mers:
    python bedFileFromFai.py reference.fasta.fai 80 81 > bedFile.bed
  #for 60mers with a 20bp overlap:
    python bedFileFromFai.py reference.fasta.fai 60 20 > bedFile.bed


```

With the newly formed bed file (or a different one with specific regions of interest), use extractMarkersFromBam.py to extract genetic markers from a bam alignment.
```
useage:
  python makeMarkersFromBam.py bedFile bamFile kmerSize > markerFile.dat
  #note, kmerSize should be the same as the length or the marker that was previously specified
example:
  python makeMarkersFromBAM.py regionsOfInterest.bed output.bam 80 > markerFile.dat
```
##############################

Script: ```uniqueMarkersFromFasta.py```

Explanation:

During the chromosome assembly process, we frequently encountered situations where contigs aligned to multiple places equally well in the genome.  Some instances, this was due to sugarcane genome biology, while others are caused by simply aligning repetitive sequences to a complex polyploid genome.  

To solve this problem, we developed a script:```uniqueMarkersFromFasta.py``` that parses a fasta file and generates only unique sequences for alignment, with a fasta header that provides the locus information for where the sequence was generated from.  We found this extremely useful during the assembly process, because 
  1) repetitive markers are excluded, focusing only on unique regions of the genome that should align well
  2) Once repetitive regions were removed, alignments completed more quickly.

The script outputs a fasta file, so it can be aligned to another fasta file using the appropriate aligner (BWA; minimap2, pblat, blast, etc..) of your choice. Additionally, if you are interested in exact matches only to another genome, you can use: ```matchSeqs.R``` to quickly output exact matches between your marker file and a genome assembly.

There are effectively 2 options for running ```uniqueMarkersFromFasta.py```.  Because you can choose both your marker size and basepair walk length, we recommend either:
  1) Selecting a short marker (example: 80bp) and a walk length of 1 bp, or
  2) Selecting larger size markers (2000bp) with a larger walk length (5000 bp)

Short markers with a walk length of 1 ensures that every sequence in the genome is considered and every marker is unique.  Use this option for fine-grain analyses; however, this takes longer to run and has larger memory requirement.  

Larger markers are useful for course-grain analysis when you need faster analysis.  Although you don't consider all possible regions of the genome when generating markers, given sufficient length, all markers are likely unique.

```
Example:
  python uniqueMarkersFromFasta.py reference.fasta 80 1 > fineGrainMarkers.fasta
    #80bp marker size
    #1 bp walk length

Example:
  python uniqueMarkersFromFasta.py reference.fasta 1000 500 > courseGrainMarkers.fasta
    #1000bp marker size
    #500 bp walk length (each marker will have a 500bp overlap with each other)
```

#######################

script: ```parseMummerToPlot.py```

Explanation:

As mentioned in the Detailed Assembly Document, chromosome construction required the manual inspection of thousands of plots.  The simple reason for this was it was easier to generate a simple plot and quickly look at it, rather than try to bioinformatically code for all complexities that arose from sugarcane genome biology.  Most of the plots in the  Detailed Assembly Methods are simply _ad hoc_ seaborn scatterplots, where the color hue is derived from marker information. However, mummer alignments and plots were particularly useful during the assembly process as they: 
1) Provided standardized outputs for easy analysis and interpretation
2) Have good online documentation and tutorials
3) Run quickly with parameters that can be tuned to maximize efficiency.
  
To visualize mummer alignments between two sequences, use ```parseMummerToPlot.py``` and a .coord file from mummer.  While the example below uses a pairwise alignment between sequences, you can also align two multi-sequence fastas, and then specify only sequence pairs (target and query) of interest.

```
#Example mummer alignment commands:
  nucmer -l 100 --maxmatch -b 400 -p alignment seq1.fasta seq2.fasta

  dnadiff -d alignment.delta -p alignment

#get coords for visualization
  show-coords -c alignment.1delta > alignment.1coords

#Generate alignment plots:
  python parseMummerToPlot.py alignment.1coords seq1 seq2 minAlignLen
  #alignment.1coords - coords file of only 1-1 mappings between sequences.  You can also use a regular .coords file, just ensure you use -c with show-coords to generate.
  #seq1- name of sequence 1 to plot
  #seq2- name of sequence 2 to plot
  #minAlignLen = minimum alignment length to consider.  Adjust to cut down on noise.

#Example (used to make Supplemental Figure 4):
  python parseMummerToPlot.py pairwiseAlignment.1coords Chr6E Chr6E_alt 10000
```
#######################

Script: ```matchSeqs.R```

Explanation:

Pretty straightforward use case for this script.  This R code uses string matches to quickly match a query fasta sequences to a reference fasta, outputting exact match locations in the reference fasta. Please read the R script comments for exact usage.

######################

Script: ```RLEruns.R```

Explanation:

This is an example of a very versatile R script that was frequently used for R570 analyses.  ```RLEruns.R``` is the simplest version of the code, but please see analysis scripts: ```progenitorBlockRLEs.R``` and ```convertCountsToRLEs.R``` for more complex use cases.  Its purpose is to convert individual data points into run-length equivalents (or RLEs) within a genome to simplify outputs and ignore noisy regions.  

For example, to calculate progenitor blocks across chromosomes, individual exact kmer matches from progenitor specific repeats (derived from ```matchSeqs.R``` output) were converted into simple start and stop blocks using RLE script: ```progenitorBlockRLEs.R```.  Please read the R file for exact usage. 


###################

Script: ```PID_calc.R```

Explanation:

This script generates alignment scores between two protein sequences.  This script was frequently used to compare: variation within and among progenitor sequences within R570.

```
Input: alignmentFile.dat:
	Format (tab separated):
	OG	fastaHeader1	proteinSeq1 fastaHeader2 	proteinSeq2
  #OG- orthogroup number
	
example alignmentFile.dat
	10005	protein1	MKPFKQEFTFDERLQESAAM		protein2	MKPFKQEFTFDERLQESAMM

Usage (command line):
  Rscript PID_calc.R alignmentFile.dat

Output:
	10005	protein1	protein2	19	1	20	95.0%
	  #19- number of matches
	  #1- number of mismatches
	  #20- alignment length
	  #95.0%- pairwise identity 

```

####################

Script: ```genespaceCommands.R```

Explanation:

GENESPACE and its output are incredibly useful for analyzing complex polyploid genes.  The script: ```genespaceCommands.R``` houses the run commands used to compare annotations between R570 and other reference genomes, but please read the GENESPACE GitHub for useful information and tutorials for how to run the software. 

####################

script: ```GENESPACE_orthogroupParsing.R```

Explanation:

This script is used to parse GENESPACE output files and assess gene presence/absence, copy number and allelic variation within syntenic orthogroups.  The output can be highly customizable, but see comments within the script itself for more explanation. 


####################

Script: ```progenitorBlockRLEs.R```

Explanation: 

Script is a more complex use case that combines ```matchSeqs.R``` and ```RLEruns.R``` to generate blocks of progenitor derived sequences (_S. officinarum_ or _S.spontaneum_) across the R570 genome.  

In brief, exact sequence matches of kmers from progenitor specific repeats are searched, and their positions are outputed to a CSV file.  The CSV file is then used to calculate run-length equivalents (RLEs) to find the start and stop coordinates for each progenitor block, requiring 100 consecutive calls from one progenitor or the other. 

###################

Scripts: ```combinePAFsAndCount.R``` and ```convertCountsToRLEs.R```

Explanation: 

These scripts are used in conjunction with each other to determine regions of haplotype collapse across the R570 genome.  

In brief, once PacBio HiFi reads are aligned back to the genome assembly (using minimap2; see manuscript methods), ```combinePAFsAndCount.R``` gathers all PAF alignment files within a directory, combines them in memory then calculates the mean read depth per non-overlapping 10kb window.

```convertCountsToRLEs.R``` parses the output from combinePAFsAndCount.R to detect regions of haplotype collapse (example 2-4X depth coverage) across the genome.  The 10kb mean depth windows are parsed and converted into their run-length equivalents (RLEs) to detect start and stop coordinates for each region of haplotype collapse.

###################

Jupyter Notebook: ```r570_orthogroupProgenitorAnalysis_forSupp.ipynb```

Explanation:

This jupyter notebook contains the commands used to compare peptide pairwise identity of S.officinarum and S.spontaneum peptide sequences within synthetic orthogroups.  It parses the combined output of ```PID_calc.R``` for each progenitor and performs a non-parametric Mann-Whitney U-Test to compare within orthogroup variation.



