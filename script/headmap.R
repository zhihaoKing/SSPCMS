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




#f <- fread("CIBERSORT3.csv")
#g <- na.omit(unlist(f[,16]))
#g <- g[!g==""]
#names(g) <- NULL


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
  
  ####ESTIMATE results#####
#es <- lapply(ex,function(x) estimate(x)[,2])

  ####CMScaller#####
res <- lapply(ex,CMScaller,RNAseq=FALSE,FDR=0.05)
res <- lapply(res,function(x) data.table(V1=rownames(x),cms=x[,"prediction"]))

for(i in seq(res))
{
   res[[i]]$cms <- addNA(res[[i]]$cms)
   levels(res[[i]]$cms)[is.na(levels(res[[i]]$cms))] <- "Unknown"
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

setwd(dir_result)
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

