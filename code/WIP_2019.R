source("code/DE_functions.R")
library(foreach)
library(Matrix)

ncluster <- 5
load("data/DE_input.Rd")
Nperm <- 30
filtcpm=0
filtpercent=0.03
gRNA <- "BAG5_3_gene"

genedm <- dm[1:(dim(dm)[1]-76), ]
gRNAdm0 <- dm[(dim(dm)[1]-75):dim(dm)[1],]
gRNAlocus <- strsplit(gRNA, split = "_")[[1]][1]
gcount <- genedm[,gRNAdm0[gRNA,] >0 & colSums(gRNAdm0[rownames(gRNAdm0) != gRNA,]) ==0]
if (dim(gcount)[2] !=0 & gRNAlocus != "neg") {
    ncount1 <- genedm[, colnames(dm1dfagg)[dm1dfagg["neg",] >0 & nlocus==1]]
    res.edger <- run_edgeR_qlf(gcount,ncount1, filtcpm, filtpercent, perm=F)
    #res.ttest <- run_ttest(gcount,ncount1, filtcpm, filtpercent, perm=F)
    #res.deseq2 <- run_deseq2(gcount,ncount1, filtcpm, filtpercent, perm=F)
  }

summ_pvalues(res.edger$table$PValue)
#summ_pvalues(res.ttest)
#summ_pvalues(res.deseq2$pvalue)

gRNA <- "DPYD_1_gene"
load(paste0("data/gRNA_edgeR-QLF/", gRNA, "_edgeR-qlf_Neg1.Rd"))
summ_pvalues(res1$table$PValue)
summ_pvalues(permres1[[2]]$table$PValue)
summ_pvalues(permres1[[5]]$table$PValue)
summ_pvalues(permres1[[8]]$table$PValue)


gRNAfile <- "BCL11B_2_gene_edgeR-qlf_Neg1.Rd"
gRNA <- strsplit(gRNAfile, "_edgeR")[[1]][1]
gRNAlocus <- strsplit(gRNA, split = "_")[[1]][1]
load(paste0("data/gRNA_edgeR-QLF/", gRNAfile))
outpvalues <- c(res1$table[1:10,"PValue"], pbins)
outpbins <- matrix(NA, nrow=length(pbins), ncol=dim(res1$table)[2])
rownames(outpbins) <- pbins
colnames(outpbins) <- colnames(res1$table)
outm <- rbind(res1$table[1:10,], outpbins)
fdrmetric <- empiricalFDR(outpvalues, res1$table$PValue, unlist(lapply(permres1, function(x) x$table$PValue)))
empiricalPall <- empiricalPvalue(res1$table$PValue, unlist(lapply(permres1, function(x) x$table$PValue)))
empiricalPfdrall <- p.adjust(empiricalPall, method="BH")
empiricalP <- c(empiricalPall[1:10], rep(NA,length(pbins)))
empiricalP.FDR <- c(empiricalPfdrall[1:10], rep(NA,length(pbins)))
outm <- cbind(outm,empiricalP,empiricalP.FDR, fdrmetric)
outm$logFC <- -outm$logFC
print(signif(outm,2))
plotgn <- rownames(outm[1:10,][outm[1:10, "EmpiricalFDR"] <0.2,])

gn <- "PHLDA1"  # CACYBP
gcount <- dm[gn, gRNAdm0[gRNA,] >0 & colSums(gRNAdm0[rownames(gRNAdm0) != gRNA,]) == 0]
ncount1 <- dm[gn, colnames(dm1dfagg)[dm1dfagg["neg",] >0 & nlocus==1]]
a <- rbind(cbind(gcount,1),cbind(ncount1,2))
colnames(a) <- c("count","categ")
par(mfrow=c(1,3))
# boxplot(count ~categ, data = a, lwd = 2, main=paste(gRNA, gn, sep=":"),  outcol="white")
stripchart(count ~categ, data=a, vertical = TRUE,  method = "jitter", add = F, pch = 21, cex=2, col = 'blue')
qqplot( -log10(unlist(lapply(permres1, function(x) x$table$PValue))),-log10(res1$table$PValue), xlab="permuted (-log10 pvalue)", ylab="observed (-log10 pvalue)")
abline(a=0, b=1, col="red")
grid()
qqplot(-log10(runif(length(empiricalPall))),-log10(empiricalPall))
abline(a=0, b=1, col="red")
grid()


