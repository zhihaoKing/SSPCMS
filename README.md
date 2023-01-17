colorectal cancer CMS single sample predictor
==========
## Set working directory
```r
rm(list=ls())
setwd(dir_data)
```
## Load R package
```r
library(org.Hs.eg.db)
library(CMScaller)
library(estimate)
library(circlize)
library(data.table)
library(dplyr)
library(Biobase) # if input is ExpressionSet
library(irlba)
library(Matrix)
library(parallel)
library(preprocessCore)
library(tidyr)
library(threejs)
library(clusterProfiler)
library(DOSE)
library(IOBR)
library(data.table)
library(ggplot2)
library(cowplot)
library(table1)
library(ggVennDiagram)
library(estimate)
library(ComplexHeatmap)
library(CMSclassifier)
library(corrplot)
library(gridExtra)
library(enrichplot)
library(ggupset)
library(RISmed)
library(WGCNA)
library(caret)
```
## Data process
```r
subHeatmap1 <- function (emat, res, templates, keepN = TRUE, 
    classCol = getOption("subClassCol"), featureCol = getOption("subClassCol"), heatCol = NULL, ...) 
{
    if (!all(colnames(emat) == rownames(res))) 
        stop("check if emat/res match")
    emat <- emat[, keepN]
    res <- res[keepN]
    class <- res
    K <- nlevels(class)
    K1 <- nlevels(templates$class)
    #templates <- templates[!duplicated(templates$probe), ]
    id_templates <- templates$probe
    N <- order(res)
    P <- intersect(rownames(emat), id_templates)
    templates <- templates[probe%in%P, ]
    templates <- templates[order(templates$class), ]
    emat <- ematAdjust(emat[templates[,probe], N])
    if (is.null(heatCol)) 
        heatCol <- CMScaller:::subData[["hmCol"]]
    pMax = 3
    breaksLeft = seq(-pMax, -0.25, length = length(heatCol)/2)
    breaks <- c(breaksLeft, 0, abs(rev(breaksLeft)))
    emat[emat > pMax] <- pMax
    emat[emat < -pMax] <- -pMax
    xx <- seq(0, 1, length.out = ncol(emat) + 1)
    yy <- seq(0, 1, length.out = nrow(emat) + 1)
    graphics::image(x = xx, y = yy, z = t(emat), yaxt = "n", 
        xaxt = "n", useRaster = TRUE, col = heatCol, breaks = breaks, 
        xlab = "CMS class predictions", ylab = "", 
        ...)
    xb <- cumsum(c(0, (sapply(split(res[N], res[N]), 
        length)/length(N))))
    xl <- xb[-(K + 1)]
    xr <- xb[-1]
    yy <- CMScaller:::line2user(line = 1:0, side = 1)
    yy <- seq(yy[1], yy[2], length = 3)
    graphics::rect(xleft = xl, xright = xr, ybottom = yy[1], 
        ytop = yy[2], col = classCol, xpd = TRUE, border = FALSE, 
        lwd = 0.75)
     segments(xl,0,xl,1)
    xxl <- xl + (xr - xl)/2
    graphics::text(xxl, yy[1], pos = 1, levels(class), adj = 0, 
        xpd = TRUE, cex = 0.75)
    xx <- CMScaller:::line2user(line = 1:0, side = 2)
    xx <- seq(xx[1], xx[2], length = 3)
    yy <- c(0, cumsum(sapply(split(templates$probe, templates$class), 
        length)/length(templates$probe)))
    
    yb <- yy[-(K1 + 1)]
    yt <- yy[-1]
    graphics::rect(xleft = xx[1], xright = xx[2], ybottom = yb, 
        ytop = yt, col = featureCol, xpd = NA, lwd = 0.75)
    segments(0,yb,1,yb)
    yyy <- yb + (yt - yb)/2
    nms <- paste0(levels(templates$class)," genes")
    graphics::text(xx[1], yyy, pos = 2,nms, adj = c(0,0), 
        xpd = TRUE, cex = 0.75,srt=70)
}
```
## Function estScore
```r
estScore <- function(x,ref) 
{
  m <- x
  gene.names <- rownames(m)
  sample.names <- colnames(m)
  Ns <- dim(m)[2]
  Ng <- dim(m)[1]
  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  N.gs <- length(ref)
  gs.names <- names(ref)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
        gene.set <- ref[[gs.i]]
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
  score.data
}
```
##  Function foldchange
```r
foldchange <- function(Y,x,s)
{
   x <- plyr::revalue(x,c("CMS1"="1","CMS2"="0","CMS3"="0","CMS4"="1"))
   Y <- scale(Y)
   id <- !x=="Unknown"
   x <- x[id]
   x <- droplevels(x)
   x <- relevel(x,"0")
   Y <- Y[,id]
   fc <- list()
   for(i in seq(s))
   {
      si <- s[i]
      if(si %in% rownames(Y))
      {
         y <- Y[si,]
         m <- lm(y~x)
         fc[[i]] <- summary(m)$coef[2,1]
      }
      else fc[[i]] <- 0
    }
    fc <- unlist(fc)
}
```
##  Function process 
```r
process <- function(x,isCounts=FALSE)
{
  
  x <- as.matrix(x)
  if(isCounts)
    x <- limma::normalizeQuantiles(log2(x + 0.25))
  keep <- rowMeans(x>=2)>=0.5
  x <- x[keep,]
  x
}
```
## Data process
```r
###data processing###
f0 <- c("TCGA(COAD)","TCGA(READ)","GSE14333","GSE143985",
       "GSE161158","GSE17536","GSE38832","GSE87211",
       "GSE17537","GSE39582","GSE171680")##3 RNA-seq datasets, 8 microarray datasets
title <- as.list(c("TCGA-COAD","TCGA-READ","GSE14333","GSE143985","GSE161158","GSE17536",
                   "GSE38832","GSE87211","GSE17537","GSE39582","GSE171680"))
f1 <- paste0("meta_",f0,"_dfs.csv")
names(f1) <- f0
meta <- lapply(f1,fread)

#meta <- rbindlist(meta,use.names=TRUE)
f2 <- paste0("exp_",f0,".csv")
names(f2) <- f0
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
```
## CMScaller
```r
  ####CMScaller#####
res <- lapply(ex,CMScaller,RNAseq=FALSE,FDR=0.05)
res <- lapply(res,function(x) data.table(V1=rownames(x),cms=x[,"prediction"]))

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}
```
## Heatmap
```r
f <- fread("LX22_EntrezID.txt")

s1 <- as.character(f[[1]])
est <- SI_geneset[2,-1] %>% as.character
s2 <- mapIds(org.Hs.eg.db,est,'ENTREZID','SYMBOL')
ref <- list(cibersort=s1,est_immune=s2)
es <- lapply(ex,estScore,ref)

for(i in seq(res))
   res[[i]] <- cbind(res[[i]],t(es[[i]]))


x1 <- lapply(res,function(x) split(x$est_immune,x$cms))
x2 <- lapply(res,function(x) split(x$cibersort,x$cms))


########################################
########################################

cms <- templates.CMS %>% data.table

est <- SI_geneset[2,-1] %>% as.character
s <- mapIds(org.Hs.eg.db,est,'ENTREZID','SYMBOL')
est <- data.table(probe=s,class="EST immune",symbol=est)

f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
g <- mapIds(org.Hs.eg.db,s,'SYMBOL','ENTREZID')
ciber <- data.table(probe=s,class="CIBERSORT",symbol=g)
template <- rbind(cms,est,ciber)

#pdf("heatmap_coad.pdf",w=10,h=8)
#layout(mat=matrix(c(1, 1, 2, 3),nrow=2,ncol=2))
subHeatmap1(ex[[1]],res[[1]]$cms,template)
boxplot(x1[[1]],main="ESTIMATE immunescore (TCGA-COAD)")
boxplot(x2[[1]],main="modified ESTIMATE immunescore (TCGA-COAD)")
#dev.off()

#setwd(dir_result)
pdf("Figure1.pdf",w=10,h=8)
   layout(mat=matrix(c(1, 1, 2, 3),nrow=2,ncol=2))
   par(oma=c(0,0,3,0))
   subHeatmap1(ex[[1]],res[[1]]$cms,template,main="Heatmap",cex.main=1.2)
   boxplot(x1[[1]],main="ESTIMATE immune scores",cex.main=1)
   boxplot(x2[[1]],main="Modified ESTIMATE immune scores \n of CIBERSORT signiture genes",cex.main=1)
   title("TCGA (COAD)", outer = TRUE,cex.main=1.5)
dev.off()


pdf("Supp_figure1.pdf",w=10,h=8)
   for(i in 2:length(ex))
   {
   layout(mat=matrix(c(1, 1, 2, 3),nrow=2,ncol=2))
   par(oma=c(0,0,3,0))
   nmi <- names(ex)[i]
   subHeatmap1(ex[[i]],res[[i]]$cms,template,main="Heatmap",cex.main=1.2)
   boxplot(x1[[i]],main="ESTIMATE immune scores",cex.main=1)
   boxplot(x2[[i]],main="Modified ESTIMATE immune scores \n of CIBERSORT signiture genes",cex.main=1)
   title(nmi, outer = TRUE,cex.main=1.5)
   }
dev.off()
```
Figure1/Supp_figure1 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure1.pdf)
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Supp_figure1.pdf)
## Fold change calculating
```r
###estimate fold change calculating####
s <-SI_geneset[2,-1] %>% unlist
s <- mapIds(org.Hs.eg.db,s,'ENTREZID', 'SYMBOL')

y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(-y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r1 <- r[o[1:20],]
s1 <- rownames(r1)


f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])

y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(-y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r2 <- r[o[1:30],]
s2 <- rownames(r2)

r1 <- as.data.frame(r1)
r1$entrizID <- rownames(r1)
r1$symbol <- mapIds(org.Hs.eg.db,s1,'SYMBOL','ENTREZID')
l1 <- r1[,c(12,13)]
R1 <- cbind(l1,r1)[,-c(14,15)]

r2 <- as.data.frame(r2)
r2$entrizID <- rownames(r2)
r2$symbol <- mapIds(org.Hs.eg.db,s2,'SYMBOL','ENTREZID')
l2 <- r2[,c(12,13)]
R2 <- cbind(l2,r2)[,-c(14,15)]

setwd(dir_result)
write.csv(R1,file="est_top_genes.csv",row.names = FALSE)
write.csv(R2,file="ciber1_top_genes.csv",row.names = FALSE)

setwd(dir_data)
f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
#rr <- r[,1]
o <- order(rr)
r3 <- r[o[1:50],]
s3 <- rownames(r3)

r3 <- as.data.frame(r3)
r3$entrizID <- rownames(r3)
r3$symbol <- mapIds(org.Hs.eg.db,s3,'SYMBOL','ENTREZID')
l3 <- r3[,c(12,13)]
R3 <- cbind(l3,r3)[,-c(14,15)]

setwd(dir_result)
write.csv(R3,file="ciber2_top_genes.csv",row.names = FALSE)


cms <- templates.CMS %>% data.table
est <- mapIds(org.Hs.eg.db,s1,'SYMBOL','ENTREZID')
est <- data.table(probe=s1,class="EST immune",symbol=est)
ciber1 <- mapIds(org.Hs.eg.db,s2,'SYMBOL','ENTREZID')
ciber1 <- data.table(probe=s2,class="CIBERSORT1",symbol=ciber1)
ciber2 <- mapIds(org.Hs.eg.db,s3,'SYMBOL','ENTREZID')
ciber2 <- data.table(probe=s3,class="CIBERSORT2",symbol=ciber2)
template <- rbind(cms,est,ciber1,ciber2)
template <- unique(template,by="probe")
template <- template[-c(390:481),]
write.csv(template,file="template.csv",row.names = FALSE)


pdf("Figure4.pdf",w=10,h=8)
subHeatmap1(ex[[1]],res[[1]]$cms,template)
dev.off()


setwd(dir_data)
f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA


f <- fread("LX22_EntrezID.txt")
g <- f[[1]]
id <- g%in%rownames(ex[[1]])
f <- f[,-1]
f[!id,] <- NA
r=cor(f,y,use="pairwise.complete.obs")
rr <- apply(r,2,rank)
setwd(dir_result)
library(corrplot)

pdf("Figure5.pdf",w=10,h=20)
corrplot(round(r,2),is.corr=F,addCoef.col = 'black',tl.srt = 45,tl.col='black')
dev.off()
```
Figure4/Figure5 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure4.pdf)
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure5.pdf)
## Circlize
```r
#### circlize plot
templates <- CMScaller::templates.CMS

cms <- templates.CMS %>% data.table

est <- SI_geneset[2,-1] %>% as.character
s <- mapIds(org.Hs.eg.db,est,'ENTREZID','SYMBOL')
est <- data.table(probe=s,class="ESTIMATE immune",symbol=est)

f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
g <- mapIds(org.Hs.eg.db,s,'SYMBOL','ENTREZID')
ciber <- data.table(probe=s,class="CIBERSORT",symbol=g)
template <- rbind(cms,est,ciber)

l <- split(template,by="class")
l <- lapply(l,function(x) x$symbol)

mat=crossprod(table(stack(l)))
names( attr(mat, "dimnames")) <- NULL

setwd(dir_result)
pdf("Figure2.pdf",w=5,h=5)
chordDiagram(mat, grid.col = 1:6, symmetric=TRUE,keep.diagonal=TRUE)
dev.off()
```
Figure2 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure2.pdf)
## Benchmark 
```r
res <- rbindlist(res,idcol="ID")   




library(WGCNA)
library(caret)
template <- fread("template.csv")
id <- template[,probe]
x <- lapply(ex,function(x) {
   idx <- match(id,rownames(x))
   x[idx,]})
x <- do.call("cbind",x) 
nms <- names(f1)

Y <- list()
for(i in seq(nms))
{
   nm1 <- nms[-i]
   nm2 <- nms[i]
   x1 <- x[,res$ID%in%nm1]
   x2 <- x[,res$ID%in%nm2]
   res1 <- res[ID%in%nm1,]
   idk <- !res1$cms %in% "Unknown"
   res1 <- res1[idk,]
   x1 <- x1[,idk]
   res2 <- res[ID%in%nm2,]
   #y <- bicor(x2,x1,use="pairwise.complete.obs")
   y <- cor(x2,x1,use="pairwise.complete.obs")
   r <- apply(-y,1,rank)  %>% t
   r <- apply(r,1,function(x) which(x<=10)) %>% t
   r <- apply(r,1,function(x) res1$cms[x]) %>% t
   label <- res2$cms
   pred <- apply(r,1,function(x) c(mean(x=="CMS1"),mean(x=="CMS2"),mean(x=="CMS3"),mean(x=="CMS4"),mean(x=="Unknown"))) %>% t
   pred <- data.frame(pred)
   #pred <- apply(pred,1,function(x) ifelse(max(x)>0.3,which.max(x),0))
   pred <- apply(pred,1,which.max)
   pred[pred==5] <- "Unknown"
   pred[pred==1] <- "CMS1"
   pred[pred==2] <- "CMS2"
   pred[pred==3] <- "CMS3"
   pred[pred==4] <- "CMS4"
   pred <- factor(pred,levels=levels(label))
   res2[,pred:=pred]
   tmp <- confusionMatrix(res2$pred,res2$cms)
   Y[[i]] <- c(Accuracy=tmp$overall[["Accuracy"]],colMeans(tmp$byClass,na.rm=TRUE)[c(1,2,5,6,7)])
}
names(Y) <- nms
y1 <- do.call("rbind",Y)



template <- templates.CMS %>% data.table
id <- template[,probe]
x <- lapply(ex,function(x) {
   idx <- match(id,rownames(x))
   x[idx,]})
x <- do.call("cbind",x) 
nms <- names(f1)
#nms <- nms[-c(4,8)]
Y <- list()
for(i in seq(nms))
{
   nm1 <- nms[-i]
   nm2 <- nms[i]
   x1 <- x[,res$ID%in%nm1]
   x2 <- x[,res$ID%in%nm2]
   res1 <- res[ID%in%nm1,]
   idk <- !res1$cms %in% "Unknown"
   res1 <- res1[idk,]
   x1 <- x1[,idk]
   res2 <- res[ID%in%nm2,]
   #y <- bicor(x2,x1,use="pairwise.complete.obs")
   y <- cor(x2,x1,use="pairwise.complete.obs")
   r <- apply(-y,1,rank)  %>% t
   r <- apply(r,1,function(x) which(x<=10)) %>% t
   r <- apply(r,1,function(x) res1$cms[x]) %>% t
   label <- res2$cms
   pred <- apply(r,1,function(x) c(mean(x=="CMS1"),mean(x=="CMS2"),mean(x=="CMS3"),mean(x=="CMS4"),mean(x=="Unknown"))) %>% t
   pred <- data.frame(pred)
   #pred <- apply(pred,1,function(x) ifelse(max(x)>0.3,which.max(x),0))
   pred <- apply(pred,1,which.max)
   pred[pred==5] <- "Unknown"
   pred[pred==1] <- "CMS1"
   pred[pred==2] <- "CMS2"
   pred[pred==3] <- "CMS3"
   pred[pred==4] <- "CMS4"
   pred <- factor(pred,levels=levels(label))
   res2[,pred:=pred]
   tmp <- confusionMatrix(res2$pred,res2$cms)
   Y[[i]] <- c(Accuracy=tmp$overall[["Accuracy"]],colMeans(tmp$byClass,na.rm=TRUE)[c(1,2,5,6,7)])
}
names(Y) <- nms
y2 <- do.call("rbind",Y)

y1 <- melt(y1)
y1$id <- "Modified CMS model"
y2 <- melt(y2)
y2$id <- "original CMS model"
y <- rbind(y1,y2)
names(y) <- c("Dateset","Metric","Value","Model")

setwd(dir_result)
write.csv(y,file="Table4.csv",row.names=FALSE)

p <- ggplot(y, aes(x=Metric, y=Value,fill=Model)) + 
  geom_boxplot() +
  ylab("")
pdf("Figure6.pdf",w=8,h=6)
   p
dev.off()   
```
Figure6 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure6.pdf)
## CIBERSORT
```r
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
```
