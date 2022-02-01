library(topGO)
library("org.At.tair.db")

####RANKED GENE ONTOLOGY ANALYSIS####
load("variancesFilt.RData") #RData object with minimum information the gene names, proportion of  Vd, Va, Vg, Vr, Vm

univ <- read.table("../reference/universeAT.txt")
univ <- univ[,1]

####ranked for all genes
##fam0102 is the dataset
#indx is the column index with the values to test
#name.str --> specify name prefix to save the object
#d is whether I want to run on sorted decreasing values ("d"), increasing values("i"), or both ("b")
callGOranked <- function(fam0102, indx, name.str, d){
  genes2 <- NULL
  values <- NULL
  genesUp <- as.character(fam0102[,1])
  for (i in 1:length(genesUp)){ #split the names and get only the orthologues
    tmp <- strsplit(genesUp[i], "_")
    tmp <- tmp[[1]][length(tmp[[1]])]
    if (length(grep("A", tmp)) == 0){}
    else{genes2 <- c(genes2, tmp); values <- c(values, fam0102[i,indx])}
  }
  genes2 <- substr(genes2, 1, 9)
  print(length(genes2))
  #prepare the universe for the topGO
  univv <-values #get the values I need to sort
  univv <- as.numeric(as.vector(univv))#make them numbers, then vector
  names(univv) <- genes2 #the universe is the expressed genes 
  #run the topGo, for decreasing
  if (d %in% c("d", "b")){  
    print("Decreasing values...")
    univv <- sort(univv, decreasing = T) ##rank them according to the fst value. from biggest ot smallest
    ranked = function(x)return(x)
    tGO <-  new("topGOdata", description = "desc", ontology = "BP", geneSel=ranked, allGenes = univv, nodeSize =10, mapping = "org.At.tair.db", annot = annFUN.org)
    resGO <- runTest(tGO, algorithm = "elim", statistic = "ks",cutoff = 0.05, scoreOrder="decreasing", ranked=TRUE)
    GO.Original <- GenTable(tGO,resultKS=resGO,topNodes=length(genesInTerm(tGO)))
    assign(paste("decrGO", name.str, sep="."), GO.Original, envir=.GlobalEnv)
    assign(paste("decrGO.genes", name.str, sep="."), genesInTerm(tGO), envir=.GlobalEnv)
  }
  if (d %in% c("i", "b")){
    print("Increasing values...")
    univv <-values #get the values I need to sort
    univv <- as.numeric(as.vector(univv))#make them numbers, then vector
    names(univv) <- genes2 #the universe is the expressed genes 
    #run the topGo, for increasing
    univv <- sort(univv, decreasing = F) ##rank them according to the fst value. from smallest to biggest
    ranked = function(x)return(x)
    tGO1 <-  new("topGOdata", description = "desc", ontology = "BP", geneSel=ranked, allGenes = univv, nodeSize =10, mapping = "org.At.tair.db", annot = annFUN.org)
    resGO1 <- runTest(tGO1, algorithm = "elim", statistic = "ks",cutoff = 0.05, scoreOrder="increasing", ranked=TRUE)
    GO.Original1 <- GenTable(tGO1,resultKS=resGO1, topNodes=length(genesInTerm(tGO1)))
    assign(paste("incrsGO", name.str, sep="."), GO.Original1, envir=.GlobalEnv)
    assign(paste("incrGO.genes", name.str, sep="."), genesInTerm(tGO), envir=.GlobalEnv)
  }
}

callGOranked(variancesFilt, 22, "Vd.prop", "d")
callGOranked(variancesFilt, 23, "Va.prop", "d")
callGOranked(variancesFilt, 29, "Vg.prop", "d")
variancesFilt$Vrm <- variancesFilt$Vm.prop + variancesFilt$Vr.prop

variancesFilt$Ratio <- variancesFilt$Vd.prop / variancesFilt$Va.prop
callGOranked(variancesFilt, 27, "Vd.Va.ratio", "d")

####Get the pval for the GOs. Permute the gene identity.
permuteGOranked <- function(fam0102, indx, name.str, d){
  genes2 <- NULL
  values <- NULL
  genesUp <- as.character(fam0102[,1])
  for (i in 1:length(genesUp)){ #split the names and get only the orthologues
    tmp <- strsplit(genesUp[i], "_")
    tmp <- tmp[[1]][length(tmp[[1]])]
    if (length(grep("A", tmp)) == 0){}
    else{genes2 <- c(genes2, tmp); values <- c(values, fam0102[i,indx])}
  }
  genes2 <- substr(genes2, 1, 9)
  print(length(genes2))
  #prepare the universe for the topGO
  univv1 <-values #get the values I need to sort
  univv1 <- as.numeric(as.vector(univv1))#make them numbers, then vector
  alll <- data.frame("GO.ID"=NA, "Term"=NA, "Annotated"=NA, "Significant"= NA, "Expected"=NA, "resultKS"=NA)
  #mix them and then put the names so they are mixed the labels 
  for (indx in 1:1000){
    print(indx)
    univv <- sample(univv1)
    names(univv) <- genes2 #the universe is the expressed genes
    #print(head(univv))
  #run the topGo, for decreasing
    tryCatch({
      if (d %in% c("d", "b")){  
        print("Decreasing values...")
        univv <- sort(univv, decreasing = T) ##rank them according to the fst value. from biggest ot smallest
        ranked = function(x)return(x)
        tGO <-  new("topGOdata", description = "desc", ontology = "BP", geneSel=ranked, allGenes = univv, nodeSize =10, mapping = "org.At.tair.db", annot = annFUN.org)
        resGO <- runTest(tGO, algorithm = "elim", statistic = "ks",cutoff = 0.05, scoreOrder="decreasing", ranked=TRUE)
        GO.Original <- GenTable(tGO,resultKS=resGO,topNodes=length(genesInTerm(tGO)))
        alll <- rbind(alll, GO.Original)
        message("All ok.")
      }
      if (d %in% c("i", "b")){
        print("Increasing values...")
        univv <-values #get the values I need to sort
        univv <- as.numeric(as.vector(univv))#make them numbers, then vector
        names(univv) <- genes2 #the universe is the expressed genes 
        #run the topGo, for increasing
        univv <- sort(univv, decreasing = F) ##rank them according to the fst value. from smallest to biggest
        ranked = function(x)return(x)
        tGO1 <-  new("topGOdata", description = "desc", ontology = "BP", geneSel=ranked, allGenes = univv, nodeSize =10, mapping = "org.At.tair.db", annot = annFUN.org)
        resGO1 <- runTest(tGO1, algorithm = "elim", statistic = "ks",cutoff = 0.05, scoreOrder="increasing", ranked=TRUE)
        GO.Original1 <- GenTable(tGO1,resultKS=resGO1, topNodes=length(genesInTerm(tGO1)))
        alll <- rbind(alll, GO.Original1)
      }},
      error=function(cond) {
     message("Error in Ks.test")
      todel <-  data.frame("GO.ID"=NA, "Term"=NA, "Annotated"=NA, "Significant"= NA, "Expected"=NA, "resultKS"=NA)
      alll <- rbind(alll, todel)
    },
      finally={
      message("ok")
    })
    }  
 
  ##save all the values
  assign(paste("permGO", name.str, sep="."), alll, envir=.GlobalEnv)
  ##get the pvalue directly
  alll <- alll[!is.na(alll$resultKS),] #remove the NAs from the analysis
  pval <- alll[alll$resultKS < 0.05,]
  pval <- length(pval$resultKS) / length(alll$resultKS) #pval for cutoff when the 5% of permutations are lower
  print(pval)
  return(pval)
}

pvalVd.GOs <- permuteGOranked(variancesFilt, 22, "Vd.prop", "d")


####Extract the genes from the signifcant GOs and check for the transcript length
trans.length.GOs <- function(countsVariancesFilt, GO.table, GO.genes, pval){
  #deconstruct the gene names so it is only for the AT genes
  genes2 <- NULL
  values <- NULL
  genesUp <- as.character(countsVariancesFilt[,1])
  for (i in 1:length(genesUp)){ #split the names and get only the orthologues
    tmp <- strsplit(genesUp[i], "_")
    tmp <- tmp[[1]][length(tmp[[1]])]
    if (length(grep("A", tmp)) == 0){}
    else{genes2 <- c(genes2, tmp); values <- c(values, countsVariancesFilt[i,21])}
  }
  genes2 <- substr(genes2, 1, 9)
  dataT <- cbind(genes2, values)
  colnames(dataT) <- c("Genes", "Length")
  #get the significant GOs
  signVals <- GO.table[GO.table[,6] < pval,] 
  print(paste("Significant GOs:", length(signVals[,1]), sep= " "))
  trans.l <- NULL
  for (i in 2:length(signVals[,1])){ #ignore the biological_process
    temp <- NULL
    rt <- GO.genes[[signVals[i,1]]]
    temp <- cbind(dataT[dataT[,1] %in% rt,], signVals[i,1])  #get the lengths from the table above. Add the ID of the GO 
    trans.l <- rbind(trans.l, temp)
  }
  colnames(trans.l) <- c("Genes", "Length", "GO.ID")
  return(trans.l)
}

decr.GO.length.Va <- trans.length.GOs(countsVariancesFilt, decrGO.Va.prop, decrGO.genes.Va.prop, pvalVa.GOs)
decr.GO.length.Va <- as.data.frame(decr.GO.length.Va)
decr.GO.length.Va$Length <- as.numeric(as.character(decr.GO.length.Va$Length))
decr.GO.length.Va$GO.ID <- as.factor(decr.GO.length.Va$GO.ID)
summary(decr.GO.length.Va)

decr.GO.length.Vd <- trans.length.GOs(countsVariancesFilt, decrGO.Vd.prop, decrGO.genes.Vd.prop, pvalVd.GOs)
decr.GO.length.Vd <- as.data.frame(decr.GO.length.Vd)
decr.GO.length.Vd$Length <- as.numeric(as.character(decr.GO.length.Vd$Length))
decr.GO.length.Vd$GO.ID <- as.factor(decr.GO.length.Vd$GO.ID)
summary(decr.GO.length.Vd)

decr.GO.length.Vd$Variance <- "Vd"
decr.GO.length.Va$Variance <- "Va"
decr.GO.length <- rbind(decr.GO.length.Vd, decr.GO.length.Va)

ks.test(decr.GO.length.Va$Length, decr.GO.length.Vd$Length)

###what is the gene and GO overlap between the enriched groups?
drawVenn2D <- function(t1, t2, name.str, cols){
  venn.diagram(x = list(t1, t2),
               category.names = c("Va" , "Vd"),
               filename = name.str,
               output = TRUE ,
               imagetype="png" ,
               height = 1000 ,
               width = 1000 ,
               resolution = 1000,
               compression = "lzw",
               lwd = 2,
               lty = 'blank',
               fill = cols,
               cex = 0.4,
               cat.cex=0.4,
               print.mode="percent",
               sigdigs = 2
  )
}

drawVenn2D(decr.GO.length.Va$Genes, decr.GO.length.Vd$Genes, 'genesInEnrichedVaVd.png', c('black', 'plum4'))
drawVenn2D(decr.GO.length.Va$GO.ID, decr.GO.length.Vd$GO.ID, 'GosEnrichedVaVd.png', c('black', 'plum4'))

tmp <- decr.GO.length.Va[decr.GO.length.Va$GO.ID %in% decr.GO.length.Vd$GO.ID,]
summary(tmp)

###how does the values compare with random non enriched GOs (in neither of the two categories?)
u <- 0
temp <- NULL
tmp <- sample(decrGO.Vd.prop$GO.ID, replace=F)
for (i in 1:length(tmp)){
  if(u < 55){
    if (!(tmp[i] %in% decr.GO.length$GO.ID)){
      temp <- c(temp, tmp[i])
      u <- u + 1
    }
  }
  else{break}
}

length(temp)

###function to extract the transcript lengths from a specific list of GOs
trans.length.GO2s <- function(countsVariancesFilt, GO.genes,  signVals){
  #deconstruct the gene names so it is only for the AT genes
  genes2 <- NULL
  values <- NULL
  genesUp <- as.character(countsVariancesFilt[,1])
  for (i in 1:length(genesUp)){ #split the names and get only the orthologues
    tmp <- strsplit(genesUp[i], "_")
    tmp <- tmp[[1]][length(tmp[[1]])]
    if (length(grep("A", tmp)) == 0){}
    else{genes2 <- c(genes2, tmp); values <- c(values, countsVariancesFilt[i,21])}
  }
  genes2 <- substr(genes2, 1, 9)
  dataT <- cbind(genes2, values)
  colnames(dataT) <- c("Genes", "Length")
  #get the genes in GOs
  trans.l <- NULL
  for (i in signVals){ 
    temp <- NULL
    rt <- GO.genes[[i]]
    temp <- cbind(dataT[dataT[,1] %in% rt,], i)  #get the lengths from the table above. Add the ID of the GO 
    trans.l <- rbind(trans.l, temp)
  }
  colnames(trans.l) <- c("Genes", "Length", "GO.ID")
  return(trans.l)
}

randomGOs.length <- trans.length.GO2s(countsVariancesFilt, decrGO.genes.Vd.prop, temp)
randomGOs.length <- as.data.frame(randomGOs.length)
randomGOs.length$Length <- as.numeric(as.character(randomGOs.length$Length))
randomGOs.length$GO.ID <- as.factor(randomGOs.length$GO.ID)
summary(randomGOs.length)
ggplot(randomGOs.length, aes(x=GO.ID, y=Length, by=GO.ID)) + geom_boxplot()

randomGOs.length$Variance <- "Random"
decr.GO.length2 <- rbind(decr.GO.length, randomGOs.length)
ggplot(decr.GO.length2, aes(x=Variance, y=Length, fill=Variance)) + geom_boxplot() + scale_fill_manual(name=NULL, values=c("mediumseagreen","black", "plum4"))
ggplot(decr.GO.length2, aes(x=Length, fill=Variance)) + geom_density(alpha=0.7) + scale_fill_manual(name=NULL, values=c("mediumseagreen","black", "plum4"))

ks.test(decr.GO.length.Va$Length, randomGOs.length$Length)
ks.test(randomGOs.length$Length, decr.GO.length.Vd$Length)
ks.test(decr.GO.length.Va$Length, decr.GO.length.Vd$Length)

###Run a lmer model
hist(log2(decr.GO.length2$Length))
length.model <- lmer (log2(Length) ~ Variance +(1|GO.ID), REML=F, data=decr.GO.length2)
qqnorm(resid(length.model))
plot(length.model)
summary(length.model)
Anova(length.model)
cld(glht(length.model, mcp(Variance="Tukey")))$mcletters$Letters

drawVenn2D(randomGOs.length$GO.ID, decr.GO.length.Vd$GO.ID, 'GosEnrichedVdRandom.png', c('mediumseagreen', 'plum4'))
drawVenn2D(randomGOs.length$GO.ID, decr.GO.length.Va$GO.ID, 'GosEnrichedVaRandom.png', c('mediumseagreen', 'black'))
