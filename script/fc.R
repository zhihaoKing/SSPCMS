rm(list = ls())
dir_data <- "~/Nutstore\ Files/.symlinks/Nutstore/coad/data"
dir_result <- "~/Nutstore\ Files/.symlinks/Nutstore/coad/result"
setwd(dir_data)
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


process <- function(x,isCounts=FALSE)
{
  
  x <- as.matrix(x)
  if(isCounts)
    x <- limma::normalizeQuantiles(log2(x + 0.25))
  keep <- rowMeans(x>=2)>=0.5
  x <- x[keep,]
  x
}


###data processing###
f <- c("TCGA(COAD)","TCGA(READ)","GSE14333","GSE143985",
       "GSE161158","GSE17536","GSE38832","GSE87211",
       "GSE17537","GSE39582","GSE171680")##3 RNA-seq datasets, 8 microarray datasets
title <- as.list(c("TCGA-COAD","TCGA-READ","GSE14333","GSE143985","GSE161158","GSE17536",
                   "GSE38832","GSE87211","GSE17537","GSE39582","GSE171680"))
f1 <- paste0("meta_",f,"_dfs.csv")
names(f1) <- f
meta <- lapply(f1,fread)


#meta <- rbindlist(meta,use.names=TRUE)
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
  
  ####CMScaller results#####
res <- lapply(ex,CMScaller,RNAseq=FALSE,FDR=0.05)
res <- lapply(res,function(x) data.table(V1=rownames(x),cms=x[,"prediction"]))

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
