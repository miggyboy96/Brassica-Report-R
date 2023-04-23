#### SECTION 1 ####

rm(list=ls())
setwd("M:/w2k/BigDataBio/Brassica Report")  
load("brassicaData.rda")
load("candidateAnalysis.rda")
library(seqinr)

load("sequenceAnalysis.rda")

#### SECTION 2 ####

## 2.1 Importing and Manipulating sequence data ##
raw_seqs <- read.delim("http://www-users.york.ac.uk/~ah1309/BigData/data/genes.txt")
seqs <- raw_seqs[,-1]
names(seqs) <- raw_seqs[,1]

## 2.2 Subset to only include candidates ##
seqs <- seqs[which(names(seqs)%in%candidates)]

## 2.3 Write to .txt and .fasta ##
write.table(seqs, "seqs.txt")
write.fasta(sapply(seqs, s2c), names=names(seqs), file.out="seqs.fasta")

seqs <- as.data.frame(read.table("seqs.txt"))

#### Sequences of related genes ####

raw_At_seqs <- read.delim("At_seqs.txt",sep="\t", header=F)
At_seqs <- raw_At_seqs[,-1]
rownames(At_seqs) <- raw_At_seqs[,1]
colnames(At_seqs) <- c("Description", "Sequence")

#### Hit Table Alignment Results ####

results.blast <- read.csv("Alignment-HitTable.csv", header=F)
colnames(results.blast) <- c("Query (Candidate)",
                             "Subject Accession",
                             "% Identity",
                             "Alignment length",
                             "Mismathes",
                             "Gap Opens",
                             "Query start",
                             "Query end",
                             "Subject start",
                             "Subject end",
                             "E-value",
                             "Bit score")
accession <- unique(results.blast$`Subject Accession`)

save(seqs, results.blast, accession, file="sequenceAnalysis.rda")
