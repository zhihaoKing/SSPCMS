Colorectal Cancer Consensus molecular subtype(CMS) single sample predictor(SSP)
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
process <- function(x,isCounts=FALSE)
{
  
  x <- as.matrix(x)
  if(isCounts)
    x <- limma::normalizeQuantiles(log2(x + 0.25))
  keep <- rowMeans(x>=2)>=0.5
  x <- x[keep,]
  x
}

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
## CMS predictor compare (TCGA_COAD)
```r
res1 <- CMScaller(emat=ex[[1]],	RNAseq=FALSE,	FDR=0.05)
res1 <- res1[,c(1,7)]
colnames(res1) <- c("CMS_prediction","FDR")

Rfcms1 <- CMSclassifier::classifyCMS(ex[[1]],method="RF")[[3]]
Rfcms1 <- Rfcms1[,c(5,6)]

SScms1 <- CMSclassifier::classifyCMS(ex[[1]],method="SSP")[[3]]
SScms1 <- SScms1[,c(13,14)]

cms_compare1 <- cbind(res1,Rfcms1,SScms1)
```
##  CMScaller
```r
res <- lapply(ex,CMScaller,RNAseq=FALSE,FDR=0.05)
res <- lapply(res,function(x) data.table(V1=rownames(x),cms=x[,"prediction"]))

save(ex,meta,res,title,f0,f1,f2,i,nms,file="res_ex.RData")
```
##  Heatmap 
```r
rm(list=ls())
setwd(dir_data)
load("~/Nutstore Files/coad/data/res_ex.Rdata")

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}

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
    nms <- levels(templates$class)
    graphics::text(xx[1], yyy, pos = 2,nms, adj = c(0,0), 
        xpd = TRUE, cex = 0.75,srt=70)
}


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

cms <- templates.CMS %>% data.table

est <- SI_geneset[2,-1] %>% as.character
s <- mapIds(org.Hs.eg.db,est,'ENTREZID','SYMBOL')
est <- data.table(probe=s,class="I141",symbol=est)

f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
g <- mapIds(org.Hs.eg.db,s,'SYMBOL','ENTREZID')
ciber <- data.table(probe=s,class="LM22",symbol=g)
template <- rbind(cms,est,ciber)

setwd(dir_result)
pdf("Figure 3.pdf",w=10,h=8)
   layout(mat=matrix(c(1, 1, 2, 3),nrow=2,ncol=2))
   par(oma=c(0,0,3,0))
   subHeatmap1(ex[[1]],res[[1]]$cms,template,main="Heatmap",cex.main=1.2)
   mtext("(a)", side=3, line=1, cex=1.5, padj=-0.45,adj=0.02)
   boxplot(x1[[1]],main="OICI scores",cex.main=1)
   mtext("(b)", side=3, line=1, cex=1.5, padj=-0.45,adj = 0.02)
   boxplot(x2[[1]],main="ESTIMATE score (LM22)",cex.main=1)
   mtext("(c)", side=3, line=1, cex=1.5, padj=-0.45,adj = 0.02)
   title("TCGA (COAD)", outer = TRUE,cex.main=1.5)
dev.off()

pdf("Supp_figure 1-10.pdf",w=10,h=8)
   for(i in 2:length(ex))
   {
   layout(mat=matrix(c(1, 1, 2, 3),nrow=2,ncol=2))
   par(oma=c(0,0,3,0))
   nmi <- names(ex)[i]
   ii <- i-1
   nnmi <- paste0("Supp_figure ",ii,"   ",nmi)
   subHeatmap1(ex[[i]],res[[i]]$cms,template,main="Heatmap",cex.main=1.2)
   mtext("(a)", side=3, line=1, cex=1.5, padj=-0.45,adj=0.02)
   boxplot(x1[[1]],main="OICI scores",cex.main=1)
   mtext("(b)", side=3, line=1, cex=1.5, padj=-0.45,adj = 0.02)
   boxplot(x2[[1]],main="ESTIMATE score (LM22)",cex.main=1)
   mtext("(c)", side=3, line=1, cex=1.5, padj=-0.45,adj = 0.02)
   title(nnmi, outer = TRUE,cex.main=1.5)
   }
dev.off()
```
## Heatmap
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
## Circlize
```r
templates <- CMScaller::templates.CMS

cms <- templates.CMS %>% data.table

est <- SI_geneset[2,-1] %>% as.character
s <- mapIds(org.Hs.eg.db,est,'ENTREZID','SYMBOL')
est <- data.table(probe=s,class="I141",symbol=est)

f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
g <- mapIds(org.Hs.eg.db,s,'SYMBOL','ENTREZID')
ciber <- data.table(probe=s,class="LM22",symbol=g)
template <- rbind(cms,est,ciber)

l <- split(template,by="class")
l <- lapply(l,function(x) x$symbol)

mat=crossprod(table(stack(l)))
names( attr(mat, "dimnames")) <- NULL

chordDiagram(mat, grid.col = 1:6, symmetric=TRUE,keep.diagonal=TRUE)
```
## K-M plot
```r
rm(list=ls())
setwd(dir_data)
load("~/Nutstore Files/coad/data/res_ex.Rdata")

rm <- list()
for (i in c(1:9)) {
  rm[[i]] <- cbind(res[[i]],meta[[i]])
}
names(rm) <- f

for(i in c(1,2)){
  rm[[i]]$dfs_days <- ifelse(rm[[i]]$dfs_days=='NA'|rm[[i]]$dfs_days==0,NA,rm[[i]]$dfs_days)
  rm[[i]] <- na.omit(rm[[i]])
}

rm[[3]] <- na.omit(rm[[3]])

rm[[4]]$vital_status <- ifelse(rm[[4]]$vital_status=='NA',NA,rm[[4]]$vital_status)
rm[[4]]$dfs_days <- ifelse(rm[[4]]$dfs_days=='0',NA,rm[[4]]$dfs_days)
rm[[4]] <- na.omit(rm[[4]])

rm[[6]] <- rm[[6]]
rm[[6]]$dfs_days <- ifelse(rm[[6]]$dfs_days=='NA',NA,rm[[6]]$dfs_days)
rm[[6]] <- na.omit(rm[[6]])

for(i in c(5,7:9)){
  rm[[i]]$dfs_days <- ifelse(rm[[i]]$dfs_days=='0',NA,rm[[i]]$dfs_days)
  rm[[i]] <- na.omit(rm[[i]])
}

rmm <- rbindlist(rm,use.names=TRUE)

rmm$group <- ifelse(rmm$cms=='CMS1'|rmm$cms=='CMS4','CMS1&4','CMS2&3') 
rmm$event <- ifelse(rmm$vital_status=='Alive',0,1)
rmm$DFS <- as.numeric(rmm$dfs_days)

rmm <- rmm[-1105,] #remove sample whose dfs is 6030 day
rmm$CMS <- rmm$cms
metainfo <- rmm[,c(3:9,12)]
colnames(metainfo) <- c("dataset","sample_id","vital_status","age","dfs_days",
                        "metastasis","group","CMS")
write.csv(metainfo,"Supp_table1.csv")

sfit <- survfit(Surv(DFS, event)~group, data=rmm)
print(sfit)

ggsurv <- ggsurvplot(sfit, data = rmm,pval = TRUE,conf.int = TRUE,
           palette = c("#E7B800", "#2E9FDF"),xlab = "DFS time in days")

customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

ggsurv$plot <- ggsurv$plot + labs(
  title    = "Kaplan-Meier plotter",
)

ggsurv <- customize_labels(
  ggsurv,
  font.title    = c(16, "bold", "darkblue"),
  font.x        = c(14, "bold.italic", "red"),
  font.y        = c(14, "bold.italic", "darkred"),
  font.xtickslab = c(12, "plain", "darkgreen"),
  font.ytickslab = c(12, "plain", "darkgreen")
)

pdf("Figure 4.pdf",w=5,h=5)
ggsurv
dev.off()
```
Figure1/Supp_figure1 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure1.pdf)
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Supp_figure1.pdf)
## Customized framework
```r
rm(list=ls())
setwd(dir_data)
load("~/Nutstore Files/coad/data/res_ex.Rdata")

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}

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

#estimate fold change calculating
s <-SI_geneset[2,-1] %>% unlist
s <- mapIds(org.Hs.eg.db,s,'ENTREZID', 'SYMBOL')

y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(-y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r1 <- r[o[1:20],]
s1 <- rownames(r1)##estimate top20

#cibersort fold change calculating
f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])

y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(-y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r2 <- r[o[1:30],]
s2 <- rownames(r2)##cibersort top30

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
est <- data.table(probe=s1,class="I141",symbol=est)
ciber1 <- mapIds(org.Hs.eg.db,s2,'SYMBOL','ENTREZID')
ciber1 <- data.table(probe=s2,class="LM22(1)",symbol=ciber1)
ciber2 <- mapIds(org.Hs.eg.db,s3,'SYMBOL','ENTREZID')
ciber2 <- data.table(probe=s3,class="LM(2)",symbol=ciber2)
template <- rbind(cms,est,ciber1,ciber2)
template <- unique(template,by="probe")
template <- template[-c(390:481),]
write.csv(template,file="template.csv",row.names = FALSE)

pdf("Supp_Figure 12.pdf",w=10,h=8)
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

pdf("Figure 5.pdf",w=10,h=20)
corrplot(round(r,2),is.corr=F,addCoef.col = 'black',tl.srt = 45,tl.col='black')
dev.off()
```
Figure4/Figure5 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure4.pdf)
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure5.pdf)
## Functional enrichment analysis
```r
rm(list=ls())
setwd(dir_data)
load("~/Nutstore Files/coad/data/res_ex.Rdata")

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}

#estimate fold change calculating
s <-SI_geneset[2,-1] %>% unlist
s <- mapIds(org.Hs.eg.db,s,'ENTREZID', 'SYMBOL')

y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(-y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r1 <- r[o[1:20],]
s1 <- rownames(r1)##estimate top20

#cibersort fold change calculating
f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])

y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(-y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r2 <- r[o[1:30],]
s2 <- rownames(r2)##cibersort top30

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

setwd(dir_data)
f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
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
template <- template[(391:482),]
m1 <- template$probe

template <- rbind(est,ciber1,ciber2)
template <- unique(template,by="probe")
m1 <- template$probe

go <- enrichGO(gene          = m1,
               OrgDb         = org.Hs.eg.db,
               ont           = "ALL",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.1,
               readable      = TRUE)
go <- setReadable(go, 'org.Hs.eg.db', 'ENTREZID')
dot_go <- dotplot(go,showCategory=10,label_format = 30) + ggtitle("GO analysis")
heat_go <- heatplot(go,showCategory=10,label_format = 30) + ggtitle("GO analysis")
plot_grid(dot_bp, heat_bp)
dot_go
heat_go

kegg <- enrichKEGG(m1,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  qvalueCutoff = 0.1,
  use_internal_data = FALSE)
kegg <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
dot_kegg <- dotplot(kegg,showCategory=10,label_format = 30) + ggtitle("KEGG analysis")
heat_kegg <- heatplot(kegg,showCategory=10,label_format = 30) + ggtitle("KEGG analysis")

setwd(dir_result)
pdf("figure 6.pdf",w=14,h=12)
plot_grid(dot_go,heat_go,dot_kegg,heat_kegg,labels = c('(a)', '(b)','(c)','(d)'))
dev.off()

GO <- go@result
write.csv(GO,"Supp_Table 3.csv")

KEGG <- kegg@result
write.csv(KEGG,"Supp_Table 4.csv")
```
Figure2 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure2.pdf)
## 100 novel biomarkers literature mining 
```r
rm(list=ls())
setwd(dir_data)
load("~/Nutstore Files/coad/data/res_ex.Rdata")

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}

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
s1 <- rownames(r1)##estimate top20

#cibersort fold change calculating####
f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])

y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(-y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r2 <- r[o[1:30],]
s2 <- rownames(r2)##cibersort top30

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

setwd(dir_data)
f <- fread("LX22_EntrezID.txt")
s <- as.character(f[[1]])
y <- mapply(function(u,v) foldchange(u,v$cms,s=s),u=ex,v=res)
y[rowMeans(y)==0,] <- NA
r <- apply(y,2,rank)
rownames(r) <- s
rr <- rowMeans(r)
o <- order(rr)
r3 <- r[o[1:50],]
s3 <- rownames(r3)

r3 <- as.data.frame(r3)
r3$entrizID <- rownames(r3)
r3$symbol <- mapIds(org.Hs.eg.db,s3,'SYMBOL','ENTREZID')
l3 <- r3[,c(12,13)]
R3 <- cbind(l3,r3)[,-c(14,15)]

cms <- templates.CMS %>% data.table
est <- mapIds(org.Hs.eg.db,s1,'SYMBOL','ENTREZID')
est <- data.table(probe=s1,class="EST immune",symbol=est)
ciber1 <- mapIds(org.Hs.eg.db,s2,'SYMBOL','ENTREZID')
ciber1 <- data.table(probe=s2,class="CIBERSORT1",symbol=ciber1)
ciber2 <- mapIds(org.Hs.eg.db,s3,'SYMBOL','ENTREZID')
ciber2 <- data.table(probe=s3,class="CIBERSORT2",symbol=ciber2)
template <- rbind(est,ciber1,ciber2)
m1 <- template$symbol
m1 <- as.list(m1)

search_topic <- list()
search_query <- list()
counts <- list()

for(i in c(1:100)){
  search_topic[[i]] <- paste0("'(",m1[[i]],") AND (colorectal cancer)'")
  search_query[[i]] <- EUtilsSummary(search_topic[[i]],db="pubmed", retmax=10000,
                              datetype='pdat', mindate=2012, maxdate=2022)
  counts[[i]] <- search_query[[i]]@count
  counts[[i]] <- as.data.frame(counts[[i]])
}
names(counts) <- m1

countss <- rbindlist(counts)
countss$symbol <- m1
colnames(countss) <- c("counts","symbol")
countss[, rank1 := frank(counts)]                                       # Create ranking variable for V1
countss <- countss[order(counts,decreasing = TRUE), ]

countss$symbol <- as.character(countss$symbol)
write.csv(countss,"table 3.csv")
```
Figure6 is shown in results
![Image text](https://github.com/zhihaoKing/SSPCMS/blob/main/result/Figure6.pdf)
## 10 key driver genes literature mining
```r
x <- c("(CCL20) AND","(CCL4) AND","(CCL18) AND","(CCL8) AND","(CXCL10) AND",
       "(CXCL11) AND","(CXCL5) AND","(CXCL9) AND","(DAPK2) AND","(S100A8) AND")
y <- c("(B cells naive)","(B cells memory)","(Plasma cells)","(T cells CD8)","(T cells CD4 naive)",
       "(T cells CD4 memory resting)","(T cells CD4 memory activated)",
       "(T cells follicular helper)","(T cells regulatory(Tregs))","(T cells gamma delta)",
       "(NK cells resting)","(NK cells activated)","(Monocytes)","(Macrophages M0)",
       "(Macrophages M1)","(Macrophages M2)","(Dendritic cells resting)","(Dendritic cells activated)",
       "(Mast cells resting)","(Mast cells activated)","(Eosinophils)","(Neutrophils)")

z <- sapply(x,paste,y)

z <- as.list(z)
search_query <- list()
counts <- list()
search_topic <- list()
for(i in c(1:220))
  search_topic[[i]] <- paste0("'",z[[i]],"'")

for(i in c(1:220))
  search_query[[i]] <- EUtilsSummary(search_topic[[i]],db="pubmed", retmax=10000,
                              datetype='pdat', mindate=2012, maxdate=2022)
                              
for(i in c(1:220))
  counts[[i]] <- search_query[[i]]@count
  
for(i in c(1:220))
  counts[[i]] <- as.data.frame(counts[[i]])  
  
names(counts) <- z

countss <- rbindlist(counts)

countss$symbol <- m1
countss$MeSH <- as.character(countss$MeSH)
colnames(countss) <- c("counts","MeSH")

write.csv(countss,"Supp_Table 5.csv")

setwd(dir_data)
lite <- t(fread("literature_mining.csv"))
colnames(lite) <- lite[1,]
lite <- lite[-1,]
lite <- as.matrix(t(lite))
lit <- matrix(as.numeric(lite), ncol = 10)
rownames(lit) <- rownames(lite)
colnames(lit) <- colnames(lite)

setwd("~/Nutstore Files/coad/result")
pdf("Figure 7.pdf",w=10,h=20)
corrplot(round(lit,2),is.corr=F,addCoef.col = 'black',tl.srt = 45,tl.col='black')
dev.off()
```
## Benchmark
```r
rm(list=ls())
setwd(dir_data)
load("~/Nutstore Files/coad/data/res_ex.Rdata")

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}

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
   y <- cor(x2,x1,use="pairwise.complete.obs")
   r <- apply(-y,1,rank)  %>% t
   r <- apply(r,1,function(x) which(x<=10)) %>% t
   r <- apply(r,1,function(x) res1$cms[x]) %>% t
   label <- res2$cms
   pred <- apply(r,1,function(x) c(mean(x=="CMS1"),mean(x=="CMS2"),mean(x=="CMS3"),mean(x=="CMS4"),mean(x=="Unknown"))) %>% t
   pred <- data.frame(pred)
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
   y <- cor(x2,x1,use="pairwise.complete.obs")
   r <- apply(-y,1,rank)  %>% t
   r <- apply(r,1,function(x) which(x<=10)) %>% t
   r <- apply(r,1,function(x) res1$cms[x]) %>% t
   label <- res2$cms
   pred <- apply(r,1,function(x) c(mean(x=="CMS1"),mean(x=="CMS2"),mean(x=="CMS3"),mean(x=="CMS4"),mean(x=="Unknown"))) %>% t
   pred <- data.frame(pred)
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
y1$id <- "Modified SSP model"
y2 <- melt(y2)
y2$id <- "original SSP model"
y <- rbind(y1,y2)
names(y) <- c("Dateset","Metric","Value","Model")

setwd(dir_result)
write.csv(y,file="Supp_Table 7.csv",row.names=FALSE)

p <- ggplot(y, aes(x=Metric, y=Value,fill=Model)) + 
  geom_boxplot() +
  ylab("")

pdf("Figure 8.pdf",w=8,h=6)
   p
dev.off()   
```
## Relative expression of 10 genes 
```r
rm(list=ls())
setwd(dir_data)
load("~/Nutstore Files/coad/data/res_ex.Rdata")

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
}

ex[[4]] <- process(ex[[4]],isCounts=TRUE)

for(i in c(1:11))
{
  rownames(ex[[i]]) <- mapIds(org.Hs.eg.db,rownames(ex[[i]]),'SYMBOL', 'ENTREZID')
}

res1 <- as.data.frame(t(res[[1]]))
colnames(res1) <- res1[1,]
res1 <- res1[-1,]
ex1 <- as.data.frame(ex[[1]])
exp <- rbind(res1,ex1)
g1 <- as.data.frame(t(exp[c("cms","CCL18","CCL20","CCL4","CCL8","CXCL10","CXCL11","CXCL5","CXCL9","DAPK2","S100A8"),]))
g1 <- g1 %>%
  pivot_longer(cols =2:11, names_to = "Gene", values_to = "Expression")
g1$group <- ifelse(g1$cms=='CMS1'|g1$cms=='CMS4','CMS1&4',g1$cms)
g1$group <- ifelse(g1$cms=='CMS2'|g1$cms=='CMS3','CMS2&3',g1$group)
g1$group <- as.factor(g1$group)
g1$Expression <- as.numeric(g1$Expression)
ggplot(data=g1,aes(x=Gene,y=Expression,fill=group))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = title[[1]])+
  scale_fill_manual(values=c("red","blue","purple")) 

g1 <- list()
g_box <- list()
exp <- list()
for(i in c(1:11))
  res[[i]] <- as.data.frame(t(res[[i]]))

for(i in c(1:11)) 
  colnames(res[[i]]) <- res[[i]][1,]

for(i in c(1:11))  
  res[[i]] <- res[[i]][-1,]

for(i in c(1:11))  
  ex[[i]] <- as.data.frame(ex[[i]])

for(i in c(1:11))  
  exp[[i]] <- rbind(res[[i]],ex[[i]])

for(i in c(1:11))
g1[[i]] <- as.data.frame(t(exp[[i]][c("cms","CCL18","CCL20","CCL4","CCL8","CXCL10",
                                        "CXCL11","CXCL5","CXCL9","DAPK2","S100A8"),]))
for(i in c(1:11)) 
  g1[[i]] <- g1[[i]] %>%
  pivot_longer(cols =2:11, names_to = "Gene", values_to = "Expression")

for(i in c(1:11)) 
  g1[[i]]$cms <- as.factor(g1[[i]]$cms)

for(i in c(1:11))  
  g1[[i]]$group <- ifelse(g1[[i]]$cms=='CMS1'|g1[[i]]$cms=='CMS4','CMS1&4',g1[[i]]$cms)

for(i in c(1:11))  
  g1[[i]]$group <- ifelse(g1[[i]]$cms=='CMS2'|g1[[i]]$cms=='CMS3','CMS2&3',g1[[i]]$group)

for(i in c(1:11))  
  g1[[i]]$group <- as.factor(g1[[i]]$group)

for(i in c(1:11))  
  g1[[i]]$Expression <- as.numeric(g1[[i]]$Expression)

for(i in c(1:11))
 g1[[i]] <- na.omit(g1[[i]])

for(i in c(1:11))  
  g_box[[i]] <- ggplot(data=g1[[i]],aes(x=Gene,y=Expression,fill=group))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(title = title[[i]])+
  scale_fill_manual(values=c("red","blue"))
 
setwd(dir_result)
pdf("Supp_figure 13.pdf",w=6,h=30) 
plot_grid(g_box[[1]],g_box[[2]],g_box[[3]],g_box[[4]],g_box[[5]],g_box[[6]],g_box[[7]],
          g_box[[8]],g_box[[9]],g_box[[10]],g_box[[11]],ncol = 1)
dev.off()
```
