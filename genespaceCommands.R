#R570-GENESPACE CODE
#please the the GENESPACE REPO for detailed tutorials
#this script just records what what commands were used to generate the 
#syntenic orthogroups.
devtools::install_github("jtlovell/GENESPACE@v0.9.4", upgrade = F)
library(GENESPACE)
library(data.table)

gpar <- init_genespace(
  # -- unique genome IDs. these can't start with a number of special character
  genomeIDs = c( "Sbicolor","r570", "Spont","STP","Sviridis"), 
  
  # -- subdirectory in the rawGenomeDir matching genomeIDs
  speciesIDs = c( "Sbicolor","r570", "Spont","STP","Sviridis"),
  
  # -- subdirectory in rawGenomeDir/speciesIDs containing the genome to use. 
  # -- for example, the human genome is here:
  # -- ~/Desktop/mammalGenespace/rawRepo/Homo_sapiens/GRCh38.p13/annotation
  
  versionIDs = c("v3.1", "v2.1", "v1.1", "v1.1", "v2.1"),
  
  # -- the genome assemblies are haploid, even though the species are diploid
  ploidy = c( 1, 7,4,1,1),
  
  # -- using a default orthofinder run as the primary database. 
  # -- alternatively could use the "fast" method, which is a bit less accurate
  orthofinderMethod = "default",
  
  # -- where should the run be stored? make sure you have read/write privs. 
  wd = "path/to/working/directory",
  
  # -- how many parallel processes should be run?
  nCores = 8, 
  
  # -- what are the unique strings that specify gff and peptide fasta files?
  gffString = "gff", 
  pepString = "protein",
  
  # -- how do you call orthofinder, diamond and MCScanX?
  path2orthofinder = "orthofinder",
  path2diamond = "diamond",
  path2mcscanx = "/pathTo/MCScanX",
  
  # -- where are the raw genomes stored? See above. 
  rawGenomeDir = "path/to/rawRepo")

#parse and format annotations

parse_phytozome(gsParam = gpar,  genomeIDs = "r570", overwrite = T)
parse_phytozome(gsParam = gpar,  genomeIDs = "Sbicolor", overwrite = T)
parse_annotations(gsParam = gpar, genomeIDs = "Spont", gffIdColumn = "ID", gffEntryType = "gene", headerSep = " ", headerEntryIndex = 1, overwrite = T)
parse_annotations(gsParam = gpar, genomeIDs = "STP", gffIdColumn = "ID", gffEntryType = "gene", headerSep = " ", headerEntryIndex = 1, overwrite = T)

####run genespace
gpar <- set_syntenyParams(
  gsParam = gpar)

gpar <- run_orthofinder(
  gsParam = gpar, 
  overwrite = F)

blks <- synteny(gsParam = gpar, overwrite = F)

pangenome(gpar, 
          genomeIDs = c( "Sbicolor","Spont","r570","STP"), 
          refGenome = "Sbicolor")

#convert the pangenome database output to wide format
pgout<-fread("Sbicolor_pangenomeDB.txt.gz")

wh <- which(pgout$isNSOrtho)
pgout$id[wh] <- paste0(pgout$id[wh], "*")

# -- flag array reps
wh <- which(!pgout$isArrayRep & !pgout$isNSOrtho)
pgout$id[wh] <- paste0(pgout$id[wh], "+")

# -- reshape to wide format
pgID <- pgChr <- pgOrd <- genome <- NULL
pgw <- dcast(
  pgout,
  pgID + pgChr + pgOrd ~ genome,
  value.var = "id",
  fun.aggregate = function(x) list(x))

# -- order by pg position
setorder(pgw, pgChr, pgOrd, na.last = T)
head(pgw)

fwrite(pgw,file="Sbicolor_pangenomeDB_wide.txt.gz", quote = F, sep='\t')

#make riparian plot (Figure 2C from manuscript)
pdf("fullRiparian_Sb_ref.pdf", height = 10, width = 24)
plot_riparianHits(
  gsParam = gpar,
  refGenome = "Sbicolor",
  genomeIDs = c( "Sbicolor","Spont","r570","STP"), 
  useOrder = T,
  labelChrBiggerThan = 100,
  minGenes2plot = 200,
  chrLabCex = 1.2,chrRectBuffer = 2,
  genomeLabCex = 1.0,
  annotatePlot = T,
  braidAlpha = 0.5,
  gapProp = 0.00005)

dev.off()