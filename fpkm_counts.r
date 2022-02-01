#!/usr/bin/env Rscript
###calculate the gene length of A.lyrata
#Alyr.gene <- read.table("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/reference/Arabidopsis_lyrata.v.1.0.37.chr.gff3", sep="\t")
Alyr.gene <- read.gff("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/reference/Arabidopsis_lyrata.v.1.0.37.chr.gff3", GFF3=TRUE)
ID <- as.character(Alyr.gene$attributes)
a <- strsplit(ID, "=|,|;")
b <- rep(NA,length(a))
for(i in 1:length(a)) b[i] <- a[[i]][2]
#c <- sub("-Protein", "", b)
#c <- sub("gene", "", b)
c <- strsplit(b, ":")
#c <- c[[1]][2]
for(i in 1:length(c)) ID[i] <- c[[i]][2]
Alyr.gene <- cbind(Alyr.gene, ID) #add ID for each row
Alyr.exon <- Alyr.gene[Alyr.gene$type=="exon",]
transcript <- unique(as.vector(Alyr.exon$ID))
#transcript <- Alyr.gene[Alyr.gene$V3=="gene",10]
transcript.length <- rep(NA, length(transcript))
for(i in 1:length(transcript))
{
  unit <- Alyr.exon[Alyr.exon$ID == transcript[i],]
  transcript.length[i] <- sum(abs(unit$end - unit$start)+1)
  if(i == trunc(i/100)*100) cat(i, "\n")
}
#length <- transcript.length
#genes <- substr(transcript, 1, nchar(transcript)-2)
#genes <- paste0("gene:", genes)
genes <- paste0("gene:", transcript)
Aly.gene.average.length = round(unlist(tapply(transcript.length, genes, mean)))

##
###read the table of read count
setwd("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/stats/")

csvprint <- function(x, y, row.names=FALSE, col.names=TRUE){
  write.table(x, file=y, append=FALSE, eol='\n', sep="\t", na = "NA", dec='.', row.names=FALSE, col.names=TRUE)
}

filelist <- list.files(pattern="counts.txt")
filelist
FPKMsum <- NULL
for (u in filelist){
  #the table must have the same names of gene length object, you can use match function
  read.table <- read.csv(u, sep="\t", header=F)
  read.number = read.table[,2]
  FPKM <- rep(NA, length(read.number))
  for (i in 1:length(read.number)){FPKM[i] <- read.number[i]/Aly.gene.average.length[read.table[i,1]]/sum(read.number)*1000*1000000}
  FPKMsum <- rbind(FPKMsum, cbind(sum(FPKM, na.rm=T), u))
  read.table <- cbind(read.table, FPKM)
  colnames(read.table) <- c("Gene", "Read Count", "FPKM")
  name = paste(u, "Count.csv", sep = "_")
  csvprint(read.table, name)
}

FPKMsum1 <- as.data.frame(FPKMsum)
colnames(FPKMsum1) <- c("FPKM", "file")
FPKMsum1$FPKM <- as.numeric(as.character(FPKMsum1$FPKM))
pdf("FPKMperSample.pdf")
barplot(FPKMsum1$FPKM)
dev.off()
hist(FPKMsum1$FPKM)
summary(FPKMsum1$FPKM)
sd(FPKMsum1$FPKM)