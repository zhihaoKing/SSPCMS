rm(list = ls())
setwd("~/Nutstore Files/coad/old/new_data")
library(Biobase) # if input is ExpressionSet
library(CMScaller)
library(dplyr)
library(irlba)
library(Matrix)
library(parallel)
library(preprocessCore)
library(tidyr)
library(threejs)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(IOBR)
library(data.table)
library(ggplot2)
library(cowplot)
library(table1)
library(ggVennDiagram)
library(estimate)
library(ComplexHeatmap)
library(circlize)

estimate <- function(dat){
  require(estimate)
  require(data.table)
  input.f='tmp_estimate_input.txt'
  output.f='tmp_estimate_gene.gct'
  output.ds='tmp_estimate_score.gct'
  write.table(dat,file = input.f,sep = '\t',quote = F)
  require(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f,id="EntrezID")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")  
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}

###new_estimata(36 genes)#######
estimateScore1 <- function (input.ds, output.ds, platform = c("affymetrix", "agilent", 
                                                              "illumina")) 
{
  stopifnot(is.character(input.ds) && length(input.ds) == 1 && 
              nzchar(input.ds))
  stopifnot(is.character(output.ds) && length(output.ds) == 
              1 && nzchar(output.ds))
  platform <- match.arg(platform)
  ds <- read.delim(input.ds, header = TRUE, sep = "\t", skip = 2, 
                   row.names = 1, blank.lines.skip = TRUE, as.is = TRUE, 
                   na.strings = "")
  descs <- ds[, 1]
  ds <- ds[-1]
  row.names <- row.names(ds)
  names <- names(ds)
  dataset <- list(ds = ds, row.names = row.names, descs = descs, 
                  names = names)
  m <- data.matrix(dataset$ds)
  gene.names <- dataset$row.names
  sample.names <- dataset$names
  Ns <- length(m[1, ])
  Ng <- length(m[, 1])
  temp <- strsplit(input.ds, split = "/")
  s <- length(temp[[1]])
  input.file.name <- temp[[1]][s]
  temp <- strsplit(input.file.name, split = ".gct")
  input.file.prefix <- temp[[1]][1]
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  idx <- c(1,11,23,25,31,41,52,59,62,64,67,69,70,76,80,86,91,93,100,102,104,
           106,107,110,112,119,122,123,125,130,132,133,135,137,139,141,142)
  SI_geneset1 <- SI_geneset[,idx]
  
  gs <- as.matrix(SI_geneset1[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(SI_geneset1)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", 
                length(gene.overlap)))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    }
    else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG/Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector/sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        }
        else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES = ES, arg.ES = arg.ES, 
                                RES = RES, indicator = TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  estimate.score <- apply(score.data, 2, sum)
  if (platform != "affymetrix") {
    score.data <- rbind(score.data, estimate.score)
    rownames(score.data) <- c("StromalScore", "ImmuneScore", 
                              "ESTIMATEScore")
  }
  else {
    convert_row_estimate_score_to_tumor_purity <- function(x) {
      stopifnot(is.numeric(x))
      cos(0.6049872018 + 0.0001467884 * x)
    }
    est.new <- NULL
    for (i in 1:length(estimate.score)) {
      est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
      est.new <- rbind(est.new, est_i)
      if (est_i >= 0) {
        next
      }
      else {
        message(paste(sample.names[i], ": out of bounds", 
                      sep = ""))
      }
    }
    colnames(est.new) <- c("TumorPurity")
    estimate.t1 <- cbind(estimate.score, est.new)
    x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 
      0
    estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
    score.data <- rbind(score.data, t(estimate.t1))
    rownames(score.data) <- c("StromalScore", "ImmuneScore", 
                              "ESTIMATEScore", "TumorPurity")
  }
  outputGCT(score.data, output.ds)
}


estimate1 <- function(dat){
  require(estimate)
  require(data.table)
  input.f='tmp_estimate_input.txt'
  output.f='tmp_estimate_gene.gct'
  output.ds='tmp_estimate_score.gct'
  write.table(dat,file = input.f,sep = '\t',quote = F)
  require(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f,id="EntrezID")
  estimateScore1(input.ds = output.f,
                 output.ds=output.ds,
                 platform="illumina")  
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}


process <- function(x,isCounts=FALSE)
{
  
  x <- as.matrix(x)
  if(isCounts)
    x <- limma::normalizeQuantiles(log2(x + 0.25))
  keep <- rowMeans(x>=2)>=0.5
  x <- x[keep,]
  x
}

f <- c("TCGA(COAD)","TCGA(READ)","GSE14333","GSE143985",
       "GSE161158","GSE17536","GSE38832","GSE87211",
       "GSE17537","GSE39582","GSE171680")##3 RNA-seq datasets, 8 microarray datasets
title <- as.list(c("TCGA-COAD","TCGA-READ","GSE14333","GSE143985","GSE161158","GSE17536",
                   "GSE38832","GSE87211","GSE17537","GSE39582","GSE171680"))
                   
                   
f2 <- paste0("exp_",f,".csv")
names(f2) <- f
ex <- lapply(f2,fread)

for(i in seq(ex))
{
  ex[[i]][,V1:=as.character(V1)]
  ex[[i]][,V1:=strsplit(V1,"///",fixed=TRUE)]
  ex[[i]] <- ex[[i]] %>% unnest(V1) %>% as.data.table(.)
}

for(i in c(3,6,7,9:11))
{
  ex[[i]][,V1:=mapIds(org.Hs.eg.db,V1,'ENTREZID', 'SYMBOL')]
}
#for(i in c(3))
#{
  #ex[[i]][,V1:=mapIds(org.Hs.eg.db,V1,'ENTREZID', 'ENSEMBL')]
#}

for(i in seq(ex))
  ex[[i]] <- na.omit(ex[[i]])

for(i in seq(ex))
{
  ex[[i]] <- ex[[i]][,lapply(.SD,sum),by=V1,.SDcols=names(ex[[i]])[!names(ex[[i]])=="V1"]]
  nms <- ex[[i]][,V1]
  ex[[i]] <- as.matrix(ex[[i]][,-"V1"])
  rownames(ex[[i]]) <- nms
}


for(i in c(1,2))
  ex[[i]] <- process(ex[[i]],isCounts=TRUE)

for(i in c(3:11))
  ex[[i]] <- process(ex[[i]],isCounts=FALSE)


####ESTIMATE results#####
#es <- lapply(ex,function(x) estimate(x)[,2])  

####CMScaller results#####
res <- lapply(ex,CMScaller,RNAseq=FALSE,FDR=0.05)
res <- lapply(res,function(x) data.table(V1=rownames(x),cms=x[,"prediction"]))

############CIBERSORT#######################

for (i in seq(res)) {
  res[[i]]$cms <- addNA(res[[i]]$cms)
}

for (i in seq(res)) {
  levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}

#for(i in seq(res))
  #res[[i]]$immune <- es[[i]] 

f0 <- as.list(f)
file0 <- list()

for(i in c(1:11)){
  file0[[i]] <- paste0("exp_",f0[[i]],".txt")
}

for(i in c(1:11)){
  ex[[i]] <- as.data.frame(ex[[i]])
  #file0 <- paste0("exp_",f0[[i]],".txt")
  write.table(ex[[i]],file = file0[[i]],sep = '\t',quote = F)
}

#####CIBERSORT result######
source("CIBERSORT.R")
ciber <- list()
for (i in c(1:11)) {
  ciber[[i]] <- CIBERSORT("LX22_EntrezID.txt", file0[[i]] , perm = 1000, QN = T)
}
names(ciber) <- f
#save(ciber,ex,res,es,file="main_result.RData")


load("~/Nutstore Files/coad/old/new_data/CIBERSORT_result.RData")


#####CIBERSORT boxplot######
for (i in c(1:11)) {
  ciber[[i]] <- as.data.frame(ciber[[i]])
}

cres <- list()
for (i in c(1:11)) {
  cres[[i]] <- cbind(res[[i]],ciber[[i]])
  cres[[i]] <- cres[[i]][,c(1,2,3:24)]
  cres[[i]] <- as.data.frame(cres[[i]])
  rownames(cres[[i]]) <- cres[[i]]$V1
  cres[[i]] <- cres[[i]][,-1]
}
names(cres) <- f

ccres <- list()
for (i in c(1:11)) {
  ccres[[i]] <- cres[[i]] %>%
    pivot_longer(cols =2:23, names_to = "celltype", values_to = "Fraction")
}

#for (i in c(1:11)) {
  #ccres[[i]]$cms <- addNA(ccres[[i]]$cms)
#}
names(ccres) <- f

p_ci <- list()###CIBERSORT plot###
for (i in c(1:11)) {
  p_ci[[i]] <- ggplot(data=ccres[[i]],aes(x=celltype,y=Fraction,fill=cms))+
               geom_boxplot()+
               theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
               labs(title = title[[i]])+
               scale_fill_manual(values=c("red","yellow","green","blue","purple")) 
}
names(p_ci) <- f


pdf("Fig_all_ciber.pdf",w=30,h=25)
plot_grid(p_ci[[1]],p_ci[[2]],p_ci[[3]],p_ci[[4]],p_ci[[5]],p_ci[[6]],p_ci[[7]],p_ci[[8]],
          p_ci[[9]],p_ci[[10]],p_ci[[11]],ncol = 3)
dev.off()





#pdf("Fig2.pdf",w=18,h=12)
#fig_ci <- plot_grid(p_ci[[1]],p_ci[[3]],p_ci[[5]],p_ci[[10]],
                    #labels=c("(a)","(b)","(c)","(d)"),ncol = 2)
#dev.off()

#cccres <- rbindlist(ccres)
#p_ci_all <- ggplot(data=cccres,aes(x=celltype,y=Fraction,fill=cms))+
            #geom_boxplot()+
            #theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
            #labs(title = "Bulk dataset (all)")+
            #scale_fill_manual(values=c("red","yellow","green","blue")) 

#pdf("Fig_all_ciber.pdf",w=30,h=25)
#plot_grid(p_ci_all,p_ci[[1]],p_ci[[2]],p_ci[[3]],p_ci[[4]],p_ci[[5]],p_ci[[6]],p_ci[[7]],p_ci[[8]],
         #p_ci[[9]],p_ci[[10]],p_ci[[11]],ncol = 3)
#dev.off()
