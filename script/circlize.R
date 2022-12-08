rm(list=ls())
dir_data <- "~/Nutstore\ Files/.symlinks/Nutstore/coad/data"
dir_result <- "~/Nutstore\ Files/.symlinks/Nutstore/coad/result"
setwd(dir_data)
library(org.Hs.eg.db)
library(CMScaller)
library(estimate)
library(circlize)
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


	