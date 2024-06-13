# script for 16S dada 2 pipeline
#following DADA2 tutorial at https://benjjneb.github.io/dada2/tutorial.html

#clear workspace
rm(list = ls())


#install code for dada 2

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly=TRUE))
  BiocManager::install("dada2")


#load packages and functions
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")


# set path to fastq files
path <- "sequencing_results/16S"  
list.files(path)

#split into forward and reverse
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE))

AACMGGATTAGATACCCKG

CMGGATTAGATACCCKGG

# read in primer sequences to check for primers
FWD <- "CMGGATTAGATACCCKGG"  ## 799F
REV <- "AGGGTTGCGCTCGTTG"  ## 1115R

# function to find complements of DNA strings
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

# find all complements of primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# filter out files with ambiguous bases before checking for primers
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)





# check how many times the primers appear 
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))




# cut out remaining primers ####
# there are a few primers. lets remove these using Cutadapt - which for windows requires a C++ compiler

cutadapt <- "C:/Users/obiew/AppData/Local/Programs/Python/Python312/Scripts/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R


# and now trim primers
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags,
                             "-m", 20, # drop reads shorter than 20 to remove 0 length reads
                             "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#check that cutting worked
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
    primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
    fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))



# now that primers are cut, we pair forward and reverse reads
cutFs <- sort(list.files(path.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2.fastq.gz", full.names = TRUE))


# Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


# inspect read quality
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#plotQualityProfile(cutFs, aggregate = TRUE)
#plotQualityProfile(cutRs, aggregate = TRUE)

# cut looking at qualityscore, forward everything is above 30
# slight drop to ~45 at 220

#for reverse, its ~40 until around 150

ShortRead::readFastq(cutFs[1])
ShortRead::readFastq(cutRs[1])


# filter and trim ####

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))


# if 220 and 105, this makes 325bp
# 799 to 1115 = 316bp minus primers of 34bp
# so overlap of ~43


#jordan trimmed at 175 and 145
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
              truncLen=c(175,145),
              maxN=0,
              maxEE=c(2,2),
              truncQ=2,
              rm.phix=TRUE,
              compress=TRUE,
              multithread=FALSE) # On Windows set multithread=FALSE
head(out)



# learn error rate ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#then plot error rate
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)



# sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")

dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")


# inspect that it went okay
dadaFs[[1]]
dadaRs[[1]]


#Merge paird reads ####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])




# construct ASV table ####

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

table(nchar(getSequences(seqtab)))




# remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# percent of merged sequence reads that are chimeras
sum(seqtab.nochim)/sum(seqtab)




# track reads through pipeline ####

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#save track
write.csv(track, file="sequencing_results/16S/track_through_pipe.csv")

#assign taxonomy
# using SILVA 138.1 @ https://zenodo.org/records/4587955


taxa <- assignTaxonomy(seqtab.nochim, "sequencing_results/16s/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#then add species based on exact taxonomic matching
taxa <- addSpecies(taxa, "sequencing_results/16s/tax/silva_species_assignment_v138.1.fa.gz")


# look at taxanomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)







#then for now we'll save them
saveRDS(taxa, "input/16s/taxa_16s.rds")

saveRDS(seqtab.nochim, "input/16s/seqtab_nochim_16s.rds")



# then polish and write out fasta file, count table, taxonomy table
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste("ASV", i, sep="_")
}



# fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "input/16s/ASVs_16s.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- asv_headers
write.csv(asv_tab, "input/16s/ASVs_counts_16s.csv")


#taxa table
asv_taxa<-taxa
row.names(asv_taxa) <- asv_headers
write.csv(asv_taxa, "input/16s/asv_taxonomy_16s.csv")

