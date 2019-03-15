library(edgeR)
library(DESeq2)
library(gtools)

run_edgeR_qlf <- function(gcount, ncount, filtcpm=10, filtpercent=0.2, perm=F) {
  # perm true of false, if true will perform a permuted version.
  coldata <- data.frame(row.names = c(colnames(gcount),colnames(ncount)),
                        condition=c(rep('G',dim(gcount)[2]),rep('N',dim(ncount)[2])))
  countall <- cbind(gcount,ncount)
  if (perm==T){
    coldata$condition <- permute(coldata$condition)
  }
  y <- DGEList(counts= countall,group=coldata$condition)

  keep <- rowSums(cpm(y)>filtcpm) >= dim(countall)[2] * filtpercent
  y <- y[keep, keep.lib.sizes=FALSE]
  print(dim(y))

  y <- calcNormFactors(y)
  group= y$samples[,"group"]
  design <- model.matrix(~group)
  y <- estimateDisp(y,design)

  fitqlf <- glmQLFit(y,design)
  qlf <- glmQLFTest(fitqlf,coef=2)
  out <- topTags(qlf, n=Inf, adjust.method = "BH")
  return(out)
}

run_ttest <- function(gcount, ncount, filtcpm=10, filtpercent=0.2, perm=F) {
  coldata <- data.frame(row.names = c(colnames(gcount),colnames(ncount)),
                        condition=c(rep('G',dim(gcount)[2]),rep('N',dim(ncount)[2])))
  countall <- cbind(gcount,ncount)
  if (perm==T){
    coldata$condition <- permute(coldata$condition)
  }
  y <- DGEList(counts= countall,group=coldata$condition)

  keep <- rowSums(cpm(y)>filtcpm) >= dim(countall)[2] * filtpercent
  y <- y[keep, keep.lib.sizes=FALSE]
  countfl <- y$counts
  pv <- rep(1,dim(countfl)[1])
  for (i in 1:dim(countfl)[1]){
    a <- t.test(countfl[i,1:dim(gcount)[2]], y=countfl[i,1:dim(ncount)[2]])
    pv[i] <- a$p.value
  }
  pv
}

run_deseq2 <- function(gcount, ncount, filtcpm=10, filtpercent=0.2, perm=F) {
  coldata <- data.frame(row.names = c(colnames(gcount),colnames(ncount)),
                        condition=c(rep('G',dim(gcount)[2]),rep('N',dim(ncount)[2])))
  countall <- cbind(gcount,ncount)
  if (perm==T){
    coldata$condition <- permute(coldata$condition)
  }
  y <- DGEList(counts= countall,group=coldata$condition)

  keep <- rowSums(cpm(y)>filtcpm) >= dim(countall)[2] * filtpercent
  y <- y[keep, keep.lib.sizes=FALSE]
  countfl <- y$counts
  dds = DESeqDataSetFromMatrix(countData = countfl, colData = coldata,design = ~condition)
  dds = estimateSizeFactors(dds)
  ddsWARD = DESeq(dds)
  resWARD = results(ddsWARD)
  return(resWARD)
}
