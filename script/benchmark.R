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