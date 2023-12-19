#R code to assess gene presence/absence, copy number and allele variation from GENESPACE analyses
#These are the codes used for analysis on the R570 genome, but the outputs can be highly customized


library(data.table)

setwd("~/Desktop/")
d <- fread("gffWithOgs.txt.gz")
#subset to include only orthogroups that contain R570 genes
r570 <- subset(d, genome == "r570")

#add
d[,hasSorghum := sum(genome == "Sbicolor" & isArrayRep) == 1, by = "og"]
d[,nR570pos := sum(genome == "r570" & isArrayRep), by = "og"]

#subset to only
ds <- subset(d, hasSorghum & nR570pos >= 1 & genome =="r570")

library(Biostrings)
#read in the R570 primary peptide file that includes flag to denote where the protein is 
#derived from, based on progenitor blocks
peps <- readAAStringSet("r570_primaryPeptides.progenitorAssignments.pipeSep.fasta")
tmp <- tstrsplit(names(peps), "|", fixed = T)
pepn <- data.table(id = tmp[[1]], subg = tmp[[2]], seq = as.character(peps))

pepd <- merge(ds, pepn, by = "id")

#columns (please note that the term 'allele' below denotes proteins that are not
#identical to one another)
#example: 10 R570 genes in an orthogroup, but 2 are identical:
#total number of genes = 10
#total number of alleles = 9 (1 identical copy)

#total genes in a syntenic orthogroup (totGenes)
#number of S.spontaneum derived genes in the syntenic orthogroup (nssGenes)
#number of S.officinarum derived genes in the syntenic orthogroup (nsoGenes)
#total number of R570 alleles in the syntenic orthogroup (totAlleles)
#total number of S.spontaneum alleles in the syntenic orthogroup (ssAlleles)
#total number of S.officinarum alleles in the syntenic orthogroup (soAlleles)
#Names of S.spontaneum alleles in the syntenic orthogroup (ssGenes)
#Names of S.officinarum alleles in the syntenic orthogroup (soGenes)
#Names of alleles in the syntenic orthogroup not assigned to a progenitor (unkGenes)
out <- pepd[,list(totGenes = .N, 
                  nssGenes = sum(subg == "Ss"),
                  nsoGenes = sum(subg == "So"),
                  totAlleles = uniqueN(seq), 
                  ssAlleles = uniqueN(seq[subg == "Ss"]),
                  soAlleles = uniqueN(seq[subg == "So"]),
                  ssGenes = list(id[subg == "Ss"]),
                  soGenes = list(id[subg == "So"]),
                  unkGenes = list(id[subg %in% c("UN", "NA")])),
            by = "og"]

fwrite(out, file = "orthogroupAlleles_sorghumAnchored.allAllowed.dat", sep = '\t')

#orthogroups that contain at least one r570 gene but not more than 8 and does 
#not need to be anchored to a sorghum gene
dt <- subset(d, nR570pos >= 1 & nR570pos <= 8 & genome =="r570")
pept <- merge(dt, pepn, by = "id")
outt <- pept[,list(totGenes = .N, 
                   nssGenes = sum(subg == "Ss"),
                   nsoGenes = sum(subg == "So"),
                   totAlleles = uniqueN(seq), 
                   ssAlleles = uniqueN(seq[subg == "Ss"]),
                   soAlleles = uniqueN(seq[subg == "So"]),
                   ssGenes = list(id[subg == "Ss"]),
                   soGenes = list(id[subg == "So"]),
                   unkGenes = list(id[subg %in% c("UN", "NA")])),
             by = "og"]
fwrite(outt, file = "orthogroupAlleles_noSorghumAnchor_1alleleAllowed.dat", sep = '\t')


#code to extract syntenic orthogroup info 
#useful if you you only want to perform multiple sequence alignment or pairwise identity
#within GENESPACE syntenic orthogroups
#
#extract regions with 6-7 homeologs in R570 for comparative genomics.  
#read in pangenome output from GENESPACE 
d <- fread("Sbicolor_pangenomeDB.txt.gz")
#extract syntenic orthogroups that are anchored to Sorghum bicolor and contain 6-7 R570 orthologs
d[,`:=`(hasSb = any(genome == "Sbicolor"), nr570 = uniqueN(ofID[genome == "r570" & isDirectSyn & isArrayRep]), nSb = uniqueN(ofID[genome == "Sbicolor" & isDirectSyn & isArrayRep])), by = c("pgChr", "pgOrd", "pgID")]
ds <- subset(d, nr570 %in% c(6,7) & nSb == 1 & hasSb & genome %in% c("Sbicolor", "r570") & isArrayRep & isDirectSyn)
ds <- subset(ds, !duplicated(paste(pgChr, pgOrd, pgID, ofID)))
dsout <- ds[,c("pgID", "id")]
fwrite(dsout, file = "pangenomeHomeologsForAlignment.dat")



