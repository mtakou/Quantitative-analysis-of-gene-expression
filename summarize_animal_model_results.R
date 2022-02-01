#Required libraries
library(reshape2)

####LOAD THE GENE VARIANCES PRODUCED WITH MCMCglmm.R.#####
###Open the files to read them####
##Note the paths are hardcoded##
csvlist <- list.files("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/", pattern="variances")

##open the files and make a list
variancesAll <- NULL
for (f in csvlist){
  temp <- read.csv(f, header=T, sep='\t')
  variancesAll <- rbind(variancesAll, temp)
}

summary(variancesAll)
variancesAll <- variancesAll[!(duplicated(variancesAll$Gene)),]

###
csvlist <- list.files("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/", pattern="flag")
flags <- NULL
for (f in csvlist){
  temp <- read.csv(f, header=T, sep='\t')
  flags <- rbind(flags, temp)
}

summary(flags)
colnames(flags) <- c("Gene", "Flag")

csvlist <- list.files("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/", pattern="gelmanValues")
gelman <- NULL
for (f in csvlist){
  temp <- read.csv(f, header=T, sep='\t')
  gelman <- rbind(gelman, temp)
}

summary(gelman)

variancesAll <- merge(variancesAll, gelman, by="Gene")

######COMPARE THE RESULTS OF MODELS, WHICH HAVE ONE FACTOR REMOVED#######
#######Compare the models, significance of Vd######
##not saved under the normal object
setwd("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/significance/")

filelist <- list.files(path="/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/significance/", pattern="modelComparison*")
model.comp <- NULL
for (i in filelist){
  temp <- read.table(i, header=TRUE)
  model.comp <- rbind(model.comp, temp)
}

model.comp$Pref.Model <- model.comp$DIC.full < model.comp$DIC.red #if it is the prefered model including Dominance then the value should be true 

var.prop <- read.table("../gene.all.variances.csv", header =T)
model.comp <- merge(var.prop, model.comp, by="Gene")
length(model.comp$Gene)
summary(model.comp)

red.model <- model.comp[model.comp$Pref.Model == FALSE,]
full.model <- model.comp[model.comp$Pref.Model == TRUE,]

plot(density(full.model$Vr.prop), col="blue")
lines(density(red.model$Vr.prop),col="red")
lines(density(full.model$Vd.prop), col="blue", lty=3)
lines(density(red.model$Vd.prop),col="red", lty=3)

#pval of fixed effect
hist(full.model$pval.Pop.full)
hist(full.model$pval.Inter.full)
summary(full.model)
sign.mother <- full.model[full.model$pval.Pop.full < (0.05/20106),]
length(sign.mother$Gene)
summary(sign.mother)
hist(sign.mother$pval.Pop.full)

###
#####Compare the models, significance of Vm#####
setwd("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/maternalComp/")

filelist <- list.files(path="/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/maternalComp/", pattern="modelComparison*")
model.comp <- NULL
for (i in filelist){
  temp <- read.table(i, header=TRUE)
  model.comp <- rbind(model.comp, temp)
}

model.comp$Pref.Model <- model.comp$DIC.full < model.comp$DIC.red #if it is the prefered model including Dominance then the value should be true 

var.prop <- read.table("../gene.all.variances.csv", header =T)
model.comp <- merge(var.prop, model.comp, by="Gene")
length(model.comp$Gene)
summary(model.comp)

red.model <- model.comp[model.comp$Pref.Model == FALSE,]
full.model <- model.comp[model.comp$Pref.Model == TRUE,]

plot(density(full.model$Vr.prop), col="blue")
lines(density(red.model$Vr.prop),col="red")
lines(density(full.model$Vd.prop), col="blue", lty=3)
lines(density(red.model$Vd.prop),col="red", lty=3)

#pval of fixed effect
hist(full.model$pval.Pop.full)
hist(full.model$pval.Inter.full)
summary(full.model)
sign.mother <- full.model[full.model$pval.Pop.full < (0.05/20106),]
length(sign.mother$Gene)
summary(sign.mother)
hist(sign.mother$pval.Pop.full)
save.image("modelComparison.RData")

####
##########Check the permutations if the distributions are different#########
setwd("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/permutations/")

csvlist.p <- list.files("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/permutations/", pattern="variances")

##open the files and make a list
variancesAll.p <- NULL
for (f in csvlist.p){
  temp <- read.csv(f, header=T, sep='\t')
  variancesAll.p <- rbind(variancesAll.p, temp)
}

summary(variancesAll.p)
variancesAll.p <- variancesAll.p[!(duplicated(variancesAll.p$Gene)),]

###
csvlist.p <- list.files("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/permutations/", pattern="flag")
flags.p <- NULL
for (f in csvlist.p){
  temp <- read.csv(f, header=T, sep='\t')
  flags.p <- rbind(flags.p, temp)
}

summary(flags.p) 
colnames(flags.p) <- c("Gene", "Flag")

##
csvlist <- list.files("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/permutations/", pattern="gelmanValues")
gelman.p <- NULL
for (f in csvlist){
  temp <- read.csv(f, header=T, sep='\t')
  gelman.p <- rbind(gelman.p, temp)
}

summary(gelman.p)
variancesAll.p <- merge(variancesAll.p, gelman.p, by="Gene")

##Filter with establish criteria
variancesPerm <- variancesAll.p[variancesAll.p$Ef.A > 499 & variancesAll.p$Ef.D > 499 & variancesAll.p$Ef.M > 499 & variancesAll.p$Ef.U > 499 & variancesAll.p$Ef.R > 499 & variancesAll.p$Inter < 1.1 & variancesAll.p$PopSp < 1.1, ]
variancesPerm <- variancesPerm[variancesPerm$Vd < 2 & variancesPerm$Va < 2 & variancesPerm$Vm < 2 & variancesPerm$Vr < 2, ]
length(variancesPerm$Gene)
summary(variancesPerm)

variancesPerm$Vd.prop <- variancesPerm$Vd / (variancesPerm$Vd + variancesPerm$Va + variancesPerm$Vr + variancesPerm$Vm)
variancesPerm$Va.prop <- variancesPerm$Va / (variancesPerm$Vd + variancesPerm$Va + variancesPerm$Vr + variancesPerm$Vm)
variancesPerm$Vr.prop <- variancesPerm$Vr / (variancesPerm$Vd + variancesPerm$Va + variancesPerm$Vr + variancesPerm$Vm)
variancesPerm$Vm.prop <- variancesPerm$Vm / (variancesPerm$Vd + variancesPerm$Va + variancesPerm$Vr + variancesPerm$Vm)
variancesPerm$Vg.prop <- (variancesPerm$Va + variancesPerm$Vd) / (variancesPerm$Vd + variancesPerm$Va + variancesPerm$Vr + variancesPerm$Vm)

###pvalue
sum(variancesPerm$Vd.prop >= mean(variancesFilt$Vd.prop)) / length(variancesPerm$Gene)
sum(variancesPerm$Va.prop >= mean(variancesFilt$Va.prop)) / length(variancesPerm$Gene)
sum(variancesPerm$Vg.prop >= mean(variancesFilt$Vg.prop)) / length(variancesPerm$Gene)
sum(variancesPerm$Vr.prop >= mean(variancesFilt$Vd.prop)) / length(variancesPerm$Gene)
sum(variancesPerm$Vm.prop >= mean(variancesFilt$Vd.prop)) / length(variancesPerm$Gene)

###FILTERING#####
variancesFilt <- variancesAll[variancesAll$Ef.A > 499 & variancesAll$Ef.D > 499 & variancesAll$Ef.M > 499 & variancesAll$Ef.U > 499 & variancesAll$Ef.R > 499 & variancesAll$Inter < 1.1 & variancesAll$PopSp < 1.1, ]
variancesFilt <- variancesFilt[variancesFilt$Vd < 2 & variancesFilt$Va < 2 & variancesFilt$Vm < 2 & variancesFilt$Vr < 2, ]
length(variancesFilt$Gene)
summary(variancesFilt)

variancesFilt$Vd.prop <- variancesFilt$Vd / (variancesFilt$Vd + variancesFilt$Va + variancesFilt$Vr + variancesFilt$Vm)
variancesFilt$Va.prop <- variancesFilt$Va / (variancesFilt$Vd + variancesFilt$Va + variancesFilt$Vr + variancesFilt$Vm)
variancesFilt$Vg.prop <- (variancesFilt$Va  + variancesFilt$Vd) / (variancesFilt$Vd + variancesFilt$Va + variancesFilt$Vr + variancesFilt$Vm)
variancesFilt$Vr.prop <- variancesFilt$Vr / (variancesFilt$Vd + variancesFilt$Va + variancesFilt$Vr + variancesFilt$Vm)
variancesFilt$Vm.prop <- variancesFilt$Vm / (variancesFilt$Vd + variancesFilt$Va + variancesFilt$Vr + variancesFilt$Vm)

####Include the RowSums filter#####
countsVariancesFilt <- load("/media/mtakou/a65c0d3c-457e-4768-916e-2773e5d7cd65/mtakou/Transcriptomes_Alyrata/dominance/countsStats.csv", header=T, sep=",")
tempFilt <- countsVariancesFilt[countsVariancesFilt$RowSums > 130,]
variancesFilt <- variancesFilt[variancesFilt$Gene %in% tempFilt$Gene,]


####CLUSTERING BASED ON GENE EXPRESSION####
###KMEANS clustering of expression
clusterMeans <- countsFinal[,2:133]
for (i in 1:132){
  clusterMeans[,i] <- log2(clusterMeans[,i]+1)
}

#look for a bend or elbow
wss <- (nrow(clusterMeans)-1) * sum(apply(clusterMeans,2,var))
for (i in 2:100) wss[i] <- sum(kmeans(clusterMeans,
                                     centers=i)$withinss)
plot(1:100, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

#get the clusters
cluster.hier <- cutree(hr.rho, 200)
clusterGenes.hier <- countsFinal
clusterGenes.hier$Cluster <- cluster.hier
clusterGenes.hier <- clusterGenes.hier[,c("Cluster", "Gene")]
clusterGenes.hier <- merge(clusterGenes.hier,countsVariancesFilt , by="Gene")
clusterGenes.hier$Cluster <- as.factor(clusterGenes.hier$Cluster)
levels(clusterGenes.hier$Cluster)

clusterGenes.hier$Vg.prop <- clusterGenes.hier$Vd.prop+clusterGenes.hier$Va.prop

par(mfrow=c(2,2))
plot(1:300, wss.h, type="b",xlab="Number of Clusters",
     ylab="Within groups sum of squares")
hist(table(clusterGenes.hier$Cluster), breaks=100, xlab="Genes per Cluster", ylab="count", main="")

#get number of genes per cluster.
table(clusterGenes.hier$Cluster)
tapply(clusterGenes.hier$Vd.prop,  clusterGenes.hier$Cluster,	 median , na.rm=T)
summ.cluster <- cbind(1:200,table(clusterGenes.hier$Cluster), tapply(clusterGenes.hier$Vd.prop,  clusterGenes.hier$Cluster,	 median , na.rm=T),tapply(clusterGenes.hier$Va.prop,  clusterGenes.hier$Cluster,	 median , na.rm=T), tapply(clusterGenes.hier$Vd.prop,  clusterGenes.hier$Cluster,	 sd , na.rm=T))
summ.cluster <- as.data.frame(summ.cluster)
colnames(summ.cluster) <- c("Cluster","Genes", "Vd", "Va","Sd")