###Script to permute and then partition the phenotypic variance to its components. Can be run iteratively for a specific number of genes on the terminal/cluster.
#example of command line:
# Rscript cluster_MCMCvPerm.R 1 600
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
start_time <- Sys.time()

####Arguements are the array number, and number of genes to be run in each file
library(MCMCglmm)
library(pedigreemm)
library(nadiv)
library(coda)

load("counts.RData")
#load("genes100.RData")

#give the table with the counts, the first gene index, last gene index
#return a table with the VA, Vd, Vr, h2, diagnostics
variancesCalc <- function(tre, st, end){
  out1 <- NULL
  flags <- NULL
  gelman <- NULL
  prior4 <- list(R=list(V=1, nu=0.02), G = list(G1 = list(V=1, nu=0.02), G2 = list(V=1, nu=0.02), G3 = list(V=1, nu=0.02)))
  #counts <- tre[st:end,] #pick a subset to iterate over
  counts <- tre[,st:end] #already transposed dataset
  pdf(paste(st, end, "gelman", "pdf", sep = "."))
  par(mfrow=c(2,2), mar=c(2, 1, 1, 1))
  for (iter in 1:length(colnames(counts))){
    #print(counts[,iter])
    temp <- sample(counts[,iter]) #for each gene permute the sample identity
    temp <- as.data.frame(temp)
    #print(temp)
    gene <- colnames(counts)[iter]
    print(gene)
    colnames(temp) <- "counts"
    #print(rownames(temp))
    print(names(counts[,iter]))
    temp$animal <- names(counts[,iter])
    temp$dom <- names(counts[,iter])
    #print(temp$animal)
    temp <- merge(temp,  toMergeFam, by="animal") 
    #print(summary(temp))
    temp$counts <- as.numeric(as.character(temp$counts))
    print(summary(temp))
    print(hist(log2(temp$counts)))
    f <- "o"
    tryCatch({
      gene.model <- MCMCglmm(log2(counts+1) ~ pop, random= ~ animal + dom + Dam, ginverse = list(dom=Dinv), start=list(QUASI=FALSE), singular.ok = FALSE, family="gaussian", prior=prior4,  pedigree=fam_matrix2, data=temp, nitt = 2200000, burnin = 200000, thin=2000, verbose = FALSE)
      print(summary(gene.model))
      EfR <- effectiveSize(gene.model$Sol)[1] #effective size of fixed effects. More than 1000, close or above 10000 is good
      EfA <- effectiveSize(gene.model$VCV)[1] #effective size of random effects
      EfD <- effectiveSize(gene.model$VCV)[2] 
      EfM <- effectiveSize(gene.model$VCV)[3]
      EfU <- effectiveSize(gene.model$VCV)[4]
      Va <- median(gene.model$VCV[,"animal"])
      Vd <- median(gene.model$VCV[,"dom"])
      Vr <- median(gene.model$VCV[,"units"])
      Vm <- median(gene.model$VCV[,"Dam"])
      h2 <- Va / (Va + Vm + Vr + Vd)
      a.cor <- autocorr.diag(gene.model$VCV)[2]
      d.cor <- autocorr.diag(gene.model$VCV)[7]
      m.cor <- autocorr.diag(gene.model$VCV)[12]
      u.cor <- autocorr.diag(gene.model$VCV)[17]
      sdA <- sd(gene.model$VCV[,"animal"])
      sdD <- sd(gene.model$VCV[,"dom"])
      sdU <- sd(gene.model$VCV[,"units"])
      sdM <- sd(gene.model$VCV[,"Dam"])
      outp <- c(as.character(gene), Va, Vd, Vm, Vr, h2, sdA, sdD, sdM, sdU, a.cor, d.cor,m.cor, u.cor,EfR, EfA, EfD, EfM, EfU)
      out1 <- rbind(out1, outp)
      f <- "e"
    }, error = function(cond){print("Error in first model.")})
    ###run the other chain
    tryCatch({
      gene.model1 <- MCMCglmm(log2(counts+1) ~ pop, random= ~ animal + dom + Dam, ginverse = list(dom=Dinv), start=list(QUASI=FALSE), singular.ok = TRUE, family="gaussian", prior=prior4,  pedigree=fam_matrix2, data=temp, nitt = 2200000, burnin = 200000, thin=2000, verbose = FALSE)
      chains <- mcmc.list(gene.model$Sol, gene.model1$Sol) #collect the chain
      gelman <- rbind(gelman, c(as.character(gene), gelman.diag(chains)$psrf[1], gelman.diag(chains)$psrf[2]))#the actual diagnostic. Has to be close to 1
      print(plot(chains, ask=F, auto.layout=F)) #plot them on top of each other
      f <- "y"
    }, error = function(cond){print("Error in second model.")})
    if (f == "y"){flags <- rbind(flags, c(as.character(gene), "pass"))} #both chains run properly
    else if (f == "e"){flags <- rbind(flags, c(as.character(gene), "fail2nd"))} #only the first passed
    else if (f == "o"){flags <- rbind(flags, c(as.character(gene), "fail"))} #both chains failed
  }
  colnames(out1) <- c("Gene", "Va", "Vd","Vm", "Vr", "h2", "Va.sd", "Vd.sd","Vm.sd", "Vr.sd", "Va.cor", "Vd.cor","Vm.cor", "Vr.cor", "Ef.R","Ef.A", "Ef.D","Ef.M", "Ef.U")
  colnames(gelman) <- c("Gene", "Inter", "PopSp")
  dev.off()
  assign("gelmanValues", gelman, envir =.GlobalEnv)
  assign("flags", flags, envir =.GlobalEnv)
  return(out1)
}

csvprint <- function(x, nameP, row.names=FALSE, col.names=TRUE){
  write.table(x, file=nameP, append=FALSE, eol='\n', sep="\t", na = "NA", dec='.', row.names=FALSE, col.names=TRUE)
}

#pick the index based on the array iteration on the cluster
#args[1] is the array being run. 
#args[2] is the number of genes I want to run in each array
end = as.numeric(as.character(args[1])) * as.numeric(as.character(args[2]))
start = end - (as.numeric(as.character(args[2])) - 1)
print(start)
print(end)
variances1.2k <- variancesCalc(counts, start, end)
csvprint(variances1.2k, paste("variances", start, end, "Perm.csv", sep="."))
csvprint(gelmanValues, paste("gelmanValues", start, end, "Perm.csv", sep="."))
csvprint(flags, paste("flags", start, end, "Perm.csv", sep="."))
end_time <- Sys.time()
end_time - start_time