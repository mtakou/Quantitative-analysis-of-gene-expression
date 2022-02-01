###load all related functions and environment. Prepatre the datasetes.
setwd("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/GWAS/")
load("../dominance/counts.RData") #where my phenotypes are
source("/home/mtakou/Documents/GWAS-master/emma.r")
source("/home/mtakou/Documents/GWAS-master/gwas.r")
source("/home/mtakou/Documents/GWAS-master/plots_gwas.r")

library(ggplot2)
library(gridExtra)
library(reshape2)

####prepare GWAS dataset####
snpMatrix <- read.table("Alyrata_GWAS.012", header=F) #columns are the SNPs, rows are the IDs
snpMatrix <- snpMatrix[,2:127586]
indvOrder <- read.table("Alyrata_GWAS.012.indv", header=F)
indvOrder <- as.vector(as.character(indvOrder$V1))
snpOrder <- read.table("Alyrata_GWAS.012.pos", header=F)
snpOrder <- paste(snpOrder$V1, snpOrder$V2, sep="- ")
rownames(snpMatrix) <- as.numeric(indvOrder) 
colnames(snpMatrix) <- snpOrder
snpMatrix <- as.matrix(snpMatrix)
snpMatrix.1 <- t(snpMatrix)

#create the kinship matrix. Function from the emma.r file which was sourced above
kinshipMatrix <- emma.kinship(snpMatrix.1, method="additive", use="all")
rownames(kinshipMatrix) <-indvOrder
colnames(kinshipMatrix) <-indvOrder

#merge counts with names
samples <- read.table("infoFile.txt", header=T)
samples$PlantID <- as.character(samples$PlantID)
temp <- counts
temp$ID <- rownames(counts)
pheno <- merge(samples, temp, by.x="SampleName", by.y="ID") #the IDs from the plant number are in the same order as in the counts

csvprint <- function(x, nameP, row.names=FALSE, col.names=TRUE){
  write.table(x, file=nameP, append=FALSE, eol='\n', sep="\t", na = "NA", dec='.', row.names=FALSE, col.names=TRUE)
}

##function to run GWAS in a loop####
gwasLoop <- function(pheno, start, end){
  outputTable <- NULL
  #pdf(paste("gwasplots", start, end, "pdf", sep="."))
  c <- 0
  for (i in start:end){
    c <- c + 1
    pheno.1 <- pheno[,c(2,i)] #pick the ID column and one gene
    pheno.1[,2] <- log2(as.numeric(as.character(pheno.1[,2])) + 1) #same transformation as for the other one
    geneName <- colnames(pheno)[i]
    pheno.1$PlantID <- as.numeric(as.character(pheno.1$PlantID))
    rownames(pheno.1) <- (pheno.1$PlantID)
    pheno.1 <- as.matrix(pheno.1)
    print(geneName)
    gwas.gene <- amm_gwas(Y = pheno.1, X = snpMatrix, K = kinshipMatrix, calculate.effect.size = F)
    #plot_gwas(gwas.gene, lower.limit = 0.50, name=geneName)
    temp <- gwas.gene[gwas.gene$Pval <= 0.05,]
    temp <- cbind(geneName, temp)
    outputTable <- rbind(outputTable, temp)
    #if (c == 50){save(outputTable, file=paste("gwashits", i - c + 1, i, "RData", sep=".")); outputTable <- NULL; c <- 0}
  }
  #dev.off()
  save(outputTable, file=paste("gwashitsch8", end, "RData", sep="."))
}

gwasLoop(pheno, 4, 6)

####Get the top hits.######
bonf.pval <- 0.05 / 127585

allgenes <- read.table("../reference/all_genes.txt")
colnames(allgenes) <- c("Chr", "Start", "End", "Length", "Gene")
allgenes.1 <- allgenes[allgenes$Chr == 1,]
allgenes.2 <- allgenes[allgenes$Chr == 2,]
allgenes.3 <- allgenes[allgenes$Chr == 3,]
allgenes.4 <- allgenes[allgenes$Chr == 4,]
allgenes.5 <- allgenes[allgenes$Chr == 5,]
allgenes.6 <- allgenes[allgenes$Chr == 6,]
allgenes.7 <- allgenes[allgenes$Chr == 7,]
allgenes.8 <- allgenes[allgenes$Chr == 8,]

###function based a specific cut off
##file, is the RData file with the saved associations
#pval is the pval cutoff. Smaller or equal to this
#min is the minimum number of snps allowed
#the function keeps only the genes with significant associations
summariseGenes <- function(file, pval){
  load(file)
  output.01 <- outputTable[outputTable$Pval <= pval,] #get only the SNPs with the association bellow 
  genes <- levels(as.factor(as.character(output.01$geneName))) ##get the gene names
  ass.genes <- vector(mode="list", length=length(genes))
  names(ass.genes) <- genes
  ass.snp <- vector(mode="list", length=length(genes))
  names(ass.snp) <- genes
  for (g in genes){ #iterate over the chromosomes of the files
    co <- NULL
    sn <- NULL
    #print(g)
    temp <- NULL
    temp <- output.01[output.01$geneName == g,]
    temp <- temp[,c("Chr", "Pos")]
    for (l in 1:length(temp$Chr)){
      ge <- NULL
      allg <- get(paste("allgenes", temp[l,1], sep=".")) #get the correct chromosome, no need to iterate over everything
      ge <- which(allg$Start <= temp[l,2] & allg$End >= temp[l,2])
      if(length(ge)==0){co <- c(co,NA); sn <- c(sn,paste(temp[l,1], temp[l,2], sep = "- "))}
      else{co <- c(co, as.character(allg[ge,5])); sn <- c(sn, paste(temp[l,1], temp[l,2], sep = "- "))}
       }
    ass.genes[[g]] <- unique(co)
    ass.snp[[g]] <- sn
  }
  flist <- list("genes" = ass.genes, "snps" = ass.snp)
  return(flist)
}

start_time <- Sys.time()
temp <- summariseGenes("gwashits.92.102.RData", 0.01)
end_time <- Sys.time()
end_time - start_time

###iterate over all the files
##only once the output table will be kept
filelist <-  list.files(path="/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/GWAS/", pattern="gwashits*")
filelist <-  list.files(path="/home/mtakou/Dropbox/margarita/Transcriptomes/GWAS/", pattern="gwashits*")
start_time <- Sys.time()
assoc.list <- NULL
assoc.snps <- NULL
for (file in filelist){
  temp <- summariseGenes(file, bonf.pval)
  assoc.list <- c(assoc.list, temp$genes)
  assoc.snps <- c(assoc.snps, temp$snps)
}
end_time <- Sys.time()
end_time - start_time

#get how many associations per gene
assoc.nums <- lengths(assoc.list)
summary(assoc.nums)

#get with how many expressed genes each association is related to
str(assoc.list)
temp <- NULL
for (i in assoc.list){temp <- c(temp, i)}
expr.nums <- as.data.frame(table(temp)) #get how many times the associations are there
summary(expr.nums)

####Combine with other datasets##########
#keep only the ones with at least one association
assoc.nums <- assoc.nums[assoc.nums > 0]

sign.assoc <- data.frame(Gene=names(assoc.nums), Assoc.Nums=assoc.nums, row.names=NULL)
sign.assoc <- sign.assoc[!duplicated(as.character(sign.assoc$Gene)),]

###read in additional files to check information
var.all <- read.table("../dominance/gene.all.variances.csv", header=T)

sign.assoc <- merge(sign.assoc, var.all, by="Gene", all.y=FALSE) 
summary(sign.assoc)
sign.assoc$Gene <- as.character(sign.assoc$Gene)

####Get the values of the variances for the genes with the significant association
#give the assoc.list, the table, and the column number of the information I need to look into. return a table to merge with the rest
getAssocValues <- function(assoc.list, var.all, indx){
  genes <- sign.assoc$Gene
  genes <- as.character(genes)
  sumStat <- NULL
  for (g in genes){
    assocs <- assoc.list[[g]] #get the associations
    if (length(assocs) == 1){tmp <- c(g, var.all[which(var.all$Gene == assocs),indx], NA, NA)}
    else if (length(assocs) > 1){
      tmp1 <- NULL 
      for (i in assocs){tmp1 <- c(tmp1, var.all[which(var.all$Gene == i),indx])}; tmp <- c(g, median(tmp1, na.rm=TRUE), mean(tmp1, na.rm=T), var(tmp1, na.rm=TRUE))}
    else{print("No assoc."); tmp <- c(g, NA, NA, NA)}
    tryCatch({sumStat <- rbind(sumStat, tmp)},warning=function(cond){message(paste("Warning message: ",g))})
  }
  return(sumStat)
}

#Vd
temp <- getAssocValues(assoc.list, var.all, 2)
colnames(temp) <- c("Gene", "med.Vd", "mean.Vd", "var.Vd")
temp <- as.data.frame(temp)
temp$med.Vd <- as.numeric(as.character(temp$med.Vd))
temp$mean.Vd <- as.numeric(as.character(temp$mean.Vd))
temp$var.Vd <- as.numeric(as.character(temp$var.Vd))
temp$Gene <- as.character(temp$Gene)
summary(temp)
sign.assoc <- merge(sign.assoc, temp, by='Gene')

#Va 
temp <- getAssocValues(assoc.list, var.all, 3)
colnames(temp) <- c("Gene", "med.Va", "mean.Va", "var.Va")
temp <- as.data.frame(temp)
temp$med.Va <- as.numeric(as.character(temp$med.Va))
temp$mean.Va <- as.numeric(as.character(temp$mean.Va))
temp$var.Va <- as.numeric(as.character(temp$var.Va))
temp$Gene <- as.character(temp$Gene)
summary(temp)
sign.assoc <- merge(sign.assoc, temp, by='Gene')

#Vg
temp <- getAssocValues(assoc.list, var.all, 4)
colnames(temp) <- c("Gene", "med.Vg", "mean.Vg", "var.Vg")
temp <- as.data.frame(temp)
temp$med.Vg <- as.numeric(as.character(temp$med.Vg))
temp$mean.Vg <- as.numeric(as.character(temp$mean.Vg))
temp$var.Vg <- as.numeric(as.character(temp$var.Vg))
temp$Gene <- as.character(temp$Gene)
summary(temp)
sign.assoc <- merge(sign.assoc, temp, by='Gene')

sign.assoc <- sign.assoc[!duplicated(sign.assoc$Gene),]
str(sign.assoc)
summary(sign.assoc)

control.Vd <- var.all[!(var.all$Gene %in% sign.assoc$Gene),]
control.Vd <- control.Vd[!(control.Vd$Gene %in% assoc.list),]

ks.test(control.Vd$Vd.prop, sign.assoc$med.Vd)
ks.test(control.Vd$Va.prop, sign.assoc$med.Va)


####Top association, SNP frequences#####
##first for the hybrids
freqPlots <- function(gene, assocgene){
  indx <- which(assoc.list[[gene]] == assocgene) #get the position of the association gene in the snps list
  snp <- assoc.snps[[gene]][indx] #get the associating snp 
  phn <- pheno[,c("PlantID",gene)] #get the phenotype
  colnames(phn) <- c("id","counts")
  #mns <- NULL
  if (length(snp) == 1){
    idx <- which(rownames(snpMatrix.1) == snp) #get which row is in the snpmatrix
    alls <- snpMatrix.1[idx,] #get all the individuals values
    tmp <- as.data.frame(alls)
    tmp$id <- rownames(tmp)
    tmp <- merge(tmp, phn, by="id") #merge with reads
    tmp$alls <- as.factor(as.character(tmp$alls))
    tmp$counts <- as.numeric(as.character(tmp$counts))
    print(boxplot(log2(tmp$counts+1) ~ tmp$alls, main=gene))
    #print(1)
    #mns <- rbind(mns, c())
  }
  else if(length(snp) > 1){
    for (s in snp){
      idx <- which(rownames(snpMatrix.1 == s)) #get which row is in the snpmatrix
      alls <- snpMatrix.1[idx,] #get all the individuals values
      tmp <- data.frame(alls)
      tmp$id <- rownames(tmp)
      tmp <- merge(tmp, phn, by="id") #merge with reads 
      tmp$alls <- as.factor(as.character(tmp$alls))
      tmp$counts <- as.numeric(as.character(tmp$counts))
      print(boxplot(log2(tmp$counts+1) ~ tmp$alls, main=gene))
      #print(2)
    }}
    else if (length(snp) == 0){print(3)}
}

freqPlots(expr.list[["gene:fgenesh2_kg.8__2505__AT5G64740.1"]][1], "gene:fgenesh2_kg.8__2505__AT5G64740.1") 
##check the SNP frequencies
filelist <- list.files(path="/home/mtakou/Dropbox/margarita/Transcriptomes/GWAS/", pattern="bial.all_genes*")
freq.sppl <- NULL
for (i in filelist){
  temp <- read.table(i, header = TRUE)
  freq.sppl <- rbind(freq.sppl, temp)
}
temp <- strsplit(as.character(freq.sppl$Freq1), ":")
temp <- do.call(rbind, temp)
freq.sppl <- cbind(freq.sppl, temp)
temp <- strsplit(as.character(freq.sppl$Freq1.1), ":")
temp <- do.call(rbind, temp)
freq.sppl <- cbind(freq.sppl, temp)
colnames(freq.sppl)[9:12] <- c("REF", "Ref.freq", "ALT", "Alt.freq")
freq.sppl$Ref.sppl <- as.numeric(as.character(freq.sppl$Ref.freq))
freq.sppl$Alt.freq <- as.numeric(as.character(freq.sppl$Alt.freq))
freq.sppl$Pos.1 <- paste(freq.sppl$Chr, freq.sppl$Pos, sep="- ")
summary(freq.sppl)

###get only the associating SNPs the frequencies
assoc.snps.hyb.freq <- NULL
for (i in assoc.snps){
  tmp <- freq.sppl[freq.sppl$Pos.1 %in% i,]
  assoc.snps.hyb.freq <- rbind(assoc.snps.hyb.freq, tmp)
}
assoc.snps.hyb.freq <- assoc.snps.hyb.freq[!(duplicated(assoc.snps.hyb.freq)),]
summary(assoc.snps.hyb.freq)
assoc.snps.hyb.freq$Ref.freq <- as.numeric(as.character(assoc.snps.hyb.freq$Ref.freq))
assoc.snps.hyb.freq$Alt.freq <- as.numeric(as.character(assoc.snps.hyb.freq$Alt.freq))

boxplot(log10((assoc.snps.hyb.freq$Ref.freq / assoc.snps.hyb.freq$Alt.freq)) ~ assoc.snps.hyb.freq$Population)
summary(assoc.snps.hyb.freq)

######Check how far away are the associated genes####
vals <- NULL
for (tmp in assoc.list){
  if (length(tmp) == 1){}
  else if (length(tmp > 1)){
    tmp1 <- poss[poss$Gene %in% tmp, ]
    if (length(tmp1$Gene) == 1){} #there are NAs
    else if(length(tmp1$Gene > 1)){
      for (i in 2:length(tmp1$Gene)){
        #print(tmp1)
        if(tmp1[i,2] == tmp1[(i-1),2]){vals <- c(vals, tmp1[i,3]-tmp1[(i-1),4])} #if they are physically on the same chromosome check the distance
        else {} #NA if not in same chromosome  
    }}}}

dist.genes <- abs(vals)
hist(dist.genes)
length(dist.genes)
summary(dist.genes)
length(dist.genes[dist.genes < 2200])

###get the genes with more than 10 associations
head(assoc.nums)
temp <- assoc.nums[assoc.nums > 10]
temp <- names(temp)
tmp <- pheno[,temp]
tmp[,1:8] <- apply(tmp[,1:8], 2, function(x)as.numeric(as.character(x)))
tmp <- cbind(pheno[,1:2], tmp)
tmp <- as.data.frame(tmp)

###
###create a plot of the associating significantly####
poss <- read.csv("../dominance/gene.orde.pos.csv", header=T) #ordered it in Excel
poss <- poss[!duplicated(poss$Gene),]
poss$Gene <- as.character(poss$Gene)

sign.assoc.matrix <- NULL
for (g in 1:length(poss$Gene)){
  tmp <- NULL
  tmp <- assoc.list[[poss[g,1]]]
  tmp1 <- NULL
 tryCatch({if (length(tmp) > 0){for (i in tmp){
   tmp1 <- c(poss[g,1], paste(poss[g,2], poss[g,3], sep="."), paste(poss[which(poss$Gene == i),2], poss[which(poss$Gene == i),3], sep="."), poss[g,2], poss[which(poss$Gene == i),1], "Sign.Assoc")
   sign.assoc.matrix <- rbind(sign.assoc.matrix, tmp1)}}
   else if (length(tmp)==0) {}
   },warning=function(cond){})
}

colnames(sign.assoc.matrix) <- c("Expr.Gene.Name", "Expressed.Gene", "Assoc.Gene", "Chr", "Assoc.Gene.Name", "Feature")
sign.assoc.matrix <- as.data.frame(sign.assoc.matrix)
summary(sign.assoc.matrix)
  
#add information about the start position of all the genes we run a GWAS
temp <- cbind(poss$Gene, paste(poss$Chr, poss$Start, sep="."), 0.5, poss$Chr, NA, "Expressed.Genes")
temp <- as.data.frame(temp)
colnames(temp) <- c("Expr.Gene.Name", "Expressed.Gene", "Assoc.Gene", "Chr", "Assoc.Gene.Name", "Feature")
sign.assoc.matrix <- rbind(sign.assoc.matrix, temp)
summary(sign.assoc.matrix)

###add information about which genes have at least one SNP in the GWAS analysis
temp <- NULL
for (s in snpOrder){
  tmp <- strsplit(s, "- ")
  ge <- which(tmp[[1]][1] == poss$Chr & as.numeric(as.character(tmp[[1]][2])) < poss$End & as.numeric(as.character(tmp[[1]][2])) > poss$Start)
  if (length(ge)==1){
    tmp1 <- c(NA, 0.5, paste(poss[ge,2], poss[ge,3], sep="."), poss[ge,2], poss[ge,1], "Genes.GWAS")
    temp <- rbind(temp, tmp1)
  }
  else{}
}


#get info of how many snps each gene has in the GWAS
gwas.genes <- temp[,5]
gwas.genes <- table(gwas.genes)
head(gwas.genes)

#remove the duplicates
colnames(temp) <- c("Expr.Gene.Name", "Expressed.Gene", "Assoc.Gene", "Chr", "Assoc.Gene.Name", "Feature")
temp <- as.data.frame(temp)
temp <- temp[!duplicated(temp$Assoc.Gene.Name),]
sign.assoc.matrix <- rbind(sign.assoc.matrix, temp)

##
sign.assoc.matrix$Assoc.Gene <- as.numeric(as.character(sign.assoc.matrix$Assoc.Gene))
sign.assoc.matrix$Expressed.Gene <- as.numeric(as.character(sign.assoc.matrix$Expressed.Gene))
sign.assoc.matrix$Chr <- as.factor(as.character(sign.assoc.matrix$Chr))
sign.assoc.matrix$Feature <- as.factor(as.character(sign.assoc.matrix$Feature))

##correct for the chromosome length
sign.assoc.matrix.2 <- sign.assoc.matrix[sign.assoc.matrix$Feature == "Sign.Assoc",]
head(sign.assoc.matrix.2)
cum.len <- c(0, 33132539, 52453403, 76917950, 100246287,121468233, 146581821, 171231018) #cummmulative length from previous chromosome
for (g in 1:length(sign.assoc.matrix.2$Expr.Gene.Name)){
  tmp <- strsplit(as.character(sign.assoc.matrix.2[g,2]), "[.]") #this is the expressed gene
  sign.assoc.matrix.2[g,2] <- as.numeric(tmp[[1]][2]) + cum.len[as.numeric(tmp[[1]][1])] #add the length of the previous chromosome to get the correct position
  tmp <- strsplit(as.character(sign.assoc.matrix.2[g,3]), "[.]") #this is the associated gene
  sign.assoc.matrix.2[g,3] <- as.numeric(tmp[[1]][2]) + cum.len[as.numeric(tmp[[1]][1])] #divide the start position per chromosome length
}

##add information on where the cis are expected to be. If its own gene would be associated with itself
temp <- sign.assoc.matrix.2
temp$Assoc.Gene <- temp$Expressed.Gene
temp <- temp[!duplicated(temp$Expr.Gene.Name),]
temp$Feature <- "Expected.Cis"
summary(temp)
sign.assoc.matrix.2 <- rbind(temp, sign.assoc.matrix.2)

mat.plot2 <- ggplot(sign.assoc.matrix.2, aes(x=Assoc.Gene, y=Expressed.Gene, col=Feature)) + geom_point(size=0.7, alpha=0.5) +
  scale_color_manual(values=c("bisque2", "cadetblue")) +
  geom_hline(yintercept=c(cum.len, sum(chr.len)), col="grey") + geom_vline(xintercept=c(cum.len, sum(chr.len)), col="grey") + 
  ylim(c(0.485,9.002)) + xlim(c(0.485,9.002)) + xlab("Associated Genes") + ylab("Expressed Genes") +
  theme(legend.position="top", text = element_text(size=15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  scale_y_continuous(breaks = cum.len + chr.len/2, labels=c("chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8")) +
  scale_x_continuous(breaks = cum.len + chr.len/2, labels=c("chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8"))

ggplotly(mat.plot2)

##
pdf("sign.associations.genome.wide.pdf")
mat.plot2
dev.off()

htmlwidgets::saveWidget(ggplotly(mat.plot2), "sign.assoc.genome.html")


####Check the CV of the phenotypes and the association and rerun the analysis for those.#####
##check one gnee
medFreqCheck <- function(gene, assocgene){
  indx <- which(assoc.list[[gene]] == assocgene) #get the position of the association gene in the snps list
  snp <- assoc.snps[[gene]][indx] #get the associating snp 
  phn <- pheno[,c("PlantID",gene)] #get the phenotype
  colnames(phn) <- c("id","counts")
  mns <- NULL
  if (length(snp) == 1){
    idx <- which(rownames(snpMatrix.1) == snp) #get which row is in the snpmatrix
    alls <- snpMatrix.1[idx,] #get all the individuals values
    tmp <- as.data.frame(alls)
    tmp$id <- rownames(tmp)
    tmp <- merge(tmp, phn, by="id") #merge with reads
    tmp$alls <- as.factor(as.character(tmp$alls))
    tmp$counts <- as.numeric(as.character(tmp$counts))
    n1 <- table(tmp$alls)
    if (length(n1) == 2){
      t1 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2]))
      mns <- rbind(mns, c(gene, assocgene, t1, n1[1], n1[2], n1[3]))
    }
    else if (length(n1) == 3){
      t1 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2]))
      t2 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3]))
      t3 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2]))
      #print(c(gene, assocgene, mean(c(t1,t2,t3), na.rm=T), n1[1], n1[2], n1[3]))
      mns <- rbind(mns, c(gene, assocgene, mean(c(t1,t2,t3), na.rm=T), n1[1], n1[2], n1[3]))
    }
  }
  else if(length(snp) > 1){
    for (s in snp){
      idx <- which(rownames(snpMatrix.1 == s)) #get which row is in the snpmatrix
      alls <- snpMatrix.1[idx,] #get all the individuals values
      tmp <- data.frame(alls)
      tmp$id <- rownames(tmp)
      tmp <- merge(tmp, phn, by="id") #merge with reads 
      tmp$alls <- as.factor(as.character(tmp$alls))
      tmp$counts <- as.numeric(as.character(tmp$counts))
      n1 <- table(tmp$alls)
      if (length(n1) == 2){
        t1 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2]))
        mns <- rbind(mns, c(gene, assocgene, t1, n1[1], n1[2], n1[3]))
      }
      else if (length(n1) == 3){
        t1 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2]))
        t2 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[1] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3]))
        t3 <-  (tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3] - tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2])/(0.5*(tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[3] + tapply(log2(tmp$counts+1), tmp$alls, mean, na.rm=T)[2]))
        #print(c(gene, assocgene, mean(c(t1,t2,t3), na.rm=T), n1[1], n1[2], n1[3]))
        mns <- rbind(mns, c(gene, assocgene, mean(c(t1,t2,t3), na.rm=T), n1[1], n1[2], n1[3]))
      }
    }}
  return(mns)
}

tmp <- medFreqCheck(expr.list[["gene:fgenesh2_kg.8__448__AT5G43745.1"]][20], "gene:fgenesh2_kg.8__2505__AT5G64740.1") 


cvFreq.Assoc <- NULL
for (u in names(expr.list)){
  if (is.na(u)){print(u)}
  else{for (i in 1:length(expr.list[[as.character(u)]])){
    temp <- medFreqCheck(expr.list[[as.character(u)]][i], as.character(u))
    cvFreq.Assoc <- rbind(cvFreq.Assoc, temp)
  }}}

cvFreq.Assoc <- as.data.frame(cvFreq.Assoc)
colnames(cvFreq.Assoc) <- c("Gene", "Assoc", "CV", "All1", "All2", "All3")
summary(cvFreq.Assoc)
cvFreq.Assoc$All1 <- (as.numeric(as.character(cvFreq.Assoc$All1)))
cvFreq.Assoc$All2 <- (as.numeric(as.character(cvFreq.Assoc$All2)))
cvFreq.Assoc$All3 <- (as.numeric(as.character(cvFreq.Assoc$All3)))
cvFreq.Assoc$CV <- (as.numeric(as.character(cvFreq.Assoc$CV)))
quantile(abs(as.numeric(as.character(cvFreq.Assoc[,3]))), c(0.25,0.40, 0.93))

hist(abs(as.numeric(as.character(cvFreq.Assoc[,3]))))
abline(v=0.02116464, col="red")

temp <- cvFreq.Assoc[cvFreq.Assoc$All1 > 7,]
temp <- temp[temp$All2 > 7,]
temp <- temp[temp$All3 > 7 | is.na(temp$All3),] #need to keep the NAs, just means they do not have a third group
length(temp$All1)

summary(sign.assoc2)
summary(sign.assoc.matrix.2)
sign.assoc2.matrix.2 <- sign.assoc.matrix.2[sign.assoc.matrix.2$Expr.Gene.Name %in% temp[,1],]
temp <- sign.assoc2.matrix.2$Expr.Gene.Name
temp <- !duplicated(temp) #405 genes

mat.plot3 <- ggplot(sign.assoc2.matrix.2, aes(x=Assoc.Gene, y=Expressed.Gene, col=Feature)) + geom_point(size=0.7, alpha=1) +
  scale_color_manual(values=c("bisque2", "cadetblue")) +
  geom_hline(yintercept=c(cum.len, sum(chr.len)), col="grey") + geom_vline(xintercept=c(cum.len, sum(chr.len)), col="grey") + 
  ylim(c(0.485,9.002)) + xlim(c(0.485,9.002)) + xlab("Associated Genes") + ylab("Expressed Genes") +
  theme(legend.position="top", text = element_text(size=10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),
        panel.background = element_blank(), axis.line = element_blank()) +
  scale_y_continuous(breaks = cum.len + chr.len/2, labels=c("chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8")) +
  scale_x_continuous(breaks = cum.len + chr.len/2, labels=c("chr1", "chr2", "chr3", "chr4", "chr5","chr6", "chr7", "chr8"))

pdf("associ.genome.wide.filter.pdf")
mat.plot3
dev.off()

htmlwidgets::saveWidget(ggplotly(mat.plot3), "sign.assoc.genome.html")