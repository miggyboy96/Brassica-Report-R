rm(list=ls())
setwd("M:/w2k/BigDataBio/Brassica Report")

#####
#FULL OSR_101 RPKM DATASET#

#1. importing

raw_RPKM <- read.table(url("https://learn-eu-central-1-prod-fleet01-xythos.content.blackboardcdn.com/62b95ee8698e9/2228803?X-Blackboard-S3-Bucket=learn-eu-central-1-prod-fleet01-xythos&X-Blackboard-Expiration=1678147200000&X-Blackboard-Signature=WhsfmXHmhIhdVaWxiOhr8A5CTqXrPyKFb6U3qsRIYVk%3D&X-Blackboard-Client-Id=301607&X-Blackboard-S3-Region=eu-central-1&response-cache-control=private%2C%20max-age%3D21600&response-content-disposition=inline%3B%20filename%2A%3DUTF-8%27%27OSR101_RPKM%25284%2529.txt&response-content-type=text%2Fplain&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEGMaDGV1LWNlbnRyYWwtMSJIMEYCIQD54u6yH037vgyMgNQCDX115XwgpmH%2BZ3pC3c0%2F6j%2F2pQIhANGDhUxKuixbTSTOFmkoNxoxvklEsJFvOV1cM5AHq7bYKr0FCBwQAhoMNjM1NTY3OTI0MTgzIgyOM9aGBr4kE1utAPcqmgWVym5KRTHeT78GAsG7r8PzQfT8bJw5z%2BIGEnSXFYms9M7R3Pw7CMY6Jqe7r21Uz87jjKaa7vSWEwo3UNO6uoCwwGx1y59AqGjEsWWkW%2BGqsOR1Wp4ldqezyatcTwEVwBIqVplU8Ltnk5UTxC1quLvN%2BtcUIuSuBZ8KKdPpO%2FoKAnxXcr3On657InHjQN4nrs%2FM8i79jjtZvpa95jGpeqhiVb8AK5s%2FWGW6eLyt1mTF6oL0eSBzEjmfCC4EvmqacbkIT3wYzkERB9CjsJX5Z%2BhW9AkxU7Xc1dpDXmeD0eYq%2B%2F20uoYlLCaYs8iYSU%2BIn%2BApCqLxLdHY%2FfLIeQzSJf%2BD7W17ivL3C4eSE%2FjWKzzOyNKF0qplHQFJd%2FRiSJUuf2e09g5lSTQYwWluqBuNp1Oas2TURb5Zh5PEQeMms3rTfwfXl68VzfAO3jRjXyGE9Cf3W0vblmngObk8tYP2qctvHFLxT%2FlBvpcWb3ulIJhQksW4r4JpQ1DiomyPHkkmDzjZkG87gpbRXF%2F6znc47EIGl6J2kA%2F7W1rqwLbRXN479QNJp9hDhKinyT7YKfETr8B6RoD4Q1Pbc77ZmPhvGvUv5C9YOxeS0geEBDejfqfrCFZTkCEfC%2FKvvNCDCdC3WBOTJ0KWr0qcj%2FONXAQ%2FDQMOcPwz7kjrp6vncrM0BrhIp2u4LrLaGs2gZyQhvoyROAJ6OXAQ01mbw7t8%2F1WNHqNnS8EHR6QBknGozK0Vzz%2FrPViKtUov7f8V8mn3dGsEPk4Js1TDHlfmhz1vdSNyQFYgxasfD0mUbQYi7KqY69WWYBFiJ5ug6IHioWArliAKgjX%2F9v2aCsPplcb1wNxHTMd2N3HZfOSI69ty0y7a8zgTVG%2B%2BdsTCjRusAJAw0O%2BYoAY6sAGxPyzhNQMjFDQb58SXNwG9X7ro6a%2F3blhqa9EiMFxap9mYbqdIV%2FkonePjhsrRgbeTi%2BVd4%2Ffo1yprjFH8DkOI4B6vPDxcr%2FILUdoThCo99vyljr1RRbtc5PdEE23BpGpfBNST%2BQAxsOOX1wyYPDAH9Y4sayVI49CwMIT3k%2FsSJcqZNeZhv4FyLRuf0K2vzjE4PQKAJv%2FbKX1%2FXj2nrJ5yk0Qz%2FAYUXNZIxQLVbmtK4w%3D%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230306T180000Z&X-Amz-SignedHeaders=host&X-Amz-Expires=21600&X-Amz-Credential=ASIAZH6WM4PLU7NMY76F%2F20230306%2Feu-central-1%2Fs3%2Faws4_request&X-Amz-Signature=47b0756ac30a4d65552b0ac9c5c7fb7ea3d4f7391c8c0143af6a32ce35ee8f40"))
write.table(raw_RPKM, "raw_RPKM.txt", quote = F, sep = "\t")
raw_RPKM <- read.table("raw_RPKM.txt", dec=".")

#2. manipulating

RPKM <- raw_RPKM[-1,-1] #remove headers and rownames
RPKM <- as.data.frame(lapply(RPKM, as.numeric)) #conversion to numeric df

#setting row and column names and storing as variables
rownames(RPKM) <- raw_RPKM[-1,1]
colnames(RPKM) <- raw_RPKM[1,-1] 
geneNames <- rownames(RPKM)
plantLines <- colnames(RPKM)

#Saving variables

save(RPKM, geneNames, plantLines, file = "brassicaData.rda")

#dim(RPKM) # 8-15 genes and 101 plant lines in data set

#####
#LOCI DATA#

#1. importing
raw_loci <- read.table("osr_dir.txt")

#2. manipulating
loci <- raw_loci[-1,-1] #removing headers and rownames
loci <- as.data.frame(lapply(loci, as.numeric))

#setting row and column names
rownames(loci) <- raw_loci[-1,1]
colnames(loci) <- raw_loci[1,-1]

#merging loci data with RPKM data
merge.loci <- merge(RPKM, loci, by='row.names')
rownames(merge.loci) <- merge.loci[,1]
merge.loci <- merge.loci[,-1]

save(merge.loci, file = "brassicaData.rda")

#####
#TRAIT DATA#

#1. importing


#2. manipulating
trait <- raw_trait[-1,]
trait <- as.data.frame(as.numeric(raw_trait[-1,-1]))

#setting names
rownames(trait) <- raw_trait[-1,1]
colnames(trait) <- raw_trait[1,-1]


#Merging trait with RPKM
merge.trait <- merge(t(RPKM), trait, by="row.names")
rownames(merge.trait) <- merge.trait[,1]
merge.trait <- merge.trait[,-1]

#dim(merge.trait) #53 plant lines #8016 genes

#saving variables to data

save(merge.trait, file = "brassicaData.rda")

#####
#GENE SEQUENCE DATA#

#1. Importing
raw_seq <- read.table(url("https://www-users.york.ac.uk/~ah1309/BigData/data/genes.txt"))

#2. Manipulating
seq <- as.data.frame(raw_seq[-1,-1])

#Setting names
rownames(seq) <- raw_seq[-1,1]
colnames(seq) <- raw_seq[1,-1]

#merging seq with RPKM
merge.seq <- merge(RPKM, seq, by='row.names')
rownames(merge.seq) <- merge.seq[,1]
merge.seq <- merge.seq[,-1]

#saving to variables data

save(merge.seq, file = "brassicaData.rda")

#### PLOTS #####

dev.new(width = 4.13, height = 5.83, unit = "in")

plots <- recordPlot
pdf(file="GLSplots.pdf")
par(mfrow=c(1,2))
barplot(merge.trait$Trait,
        ylab="GLS %",
        xlab="Plant line")
hist(merge.trait$Trait, breaks = "FD",
     xlab="GLS %")
dev.off
