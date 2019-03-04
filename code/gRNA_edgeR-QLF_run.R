source("code/DE_functions.R")
library(foreach)
library(Matrix)

ncluster <- 5
load("data/DE_input.Rd")
Nperm <- 30
filtcpm=50
filtpercent=0.2

genedm <- dm[1:(dim(dm)[1]-76), ]
gRNAdm0 <- dm[(dim(dm)[1]-75):dim(dm)[1],]

# cellpercent <-  apply(genedm0,1,function(x) length(x[x>1])/length(x))
# genedm <- genedm0[cellpercent> filt,]

cl <- parallel::makeCluster(ncluster, outfile="")
doParallel::registerDoParallel(cl)
foreach (gRNAindex=1:dim(gRNAdm0)[1], .packages = c("Matrix", "edgeR","gtools")) %dopar% {
  gRNA <- rownames(gRNAdm0)[gRNAindex]
  print(gRNA)
  gRNAlocus <- strsplit(gRNA, split = "_")[[1]][1]
  gcount <- genedm[,gRNAdm0[gRNA,] >0 & colSums(gRNAdm0[rownames(gRNAdm0) != gRNA,]) ==0]
  if (dim(gcount)[2] !=0) {
    # ncount1 <- genedm[, colnames(dm1dfagg)[dm1dfagg["neg",] >0 & nlocus==1]]
    # res1 <- run_edgeR_qlf(gcount,ncount1, filtcpm, filtpercent, perm=F)
    # permres1 <- list()
    # for (i in 1:Nperm) {
    #   permres1[[i]] <- run_edgeR_qlf(gcount,ncount1,filtcpm, filtpercent, perm=T)
    # }
    # save(res1, permres1, file=paste0("data/gRNA_edgeR-QLF/", gRNA, "_edgeR-qlf_Neg1.Rd"))

    ncount2 <- genedm[, colnames(dm1dfagg)[dm1dfagg[gRNAlocus,] == 0]] # much larger neg control, takes very long time
    res2 <- run_edgeR_qlf(gcount,ncount2, filtcpm, filtpercent, perm=F)
    permres2 <- list()
    for (i in 1:Nperm) {
      permres2[[i]] <- run_edgeR_qlf(gcount,ncount2,filtcpm, filtpercent, perm=T)
    }
    save(res2, permres2, file=paste0("data/gRNA_edgeR-QLF/", gRNA, "_edgeR-qlf_Neg2.Rd"))
  }
}

parallel::stopCluster(cl)



