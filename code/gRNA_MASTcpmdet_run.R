source("code/DE_functions.R")
library(foreach)
library(Matrix)

ncluster <- 3 # for neg2 need 12G/task
load("data/DE_input.Rd")
Nperm <- 0 # for MAST seems no need to permute
filtcpm= -1 # no filtering for MAST
filtpercent=0.2

genedm <- dm[1:(dim(dm)[1]-76), ]
gRNAdm0 <- dm[(dim(dm)[1]-75):dim(dm)[1],]

cl <- parallel::makeCluster(ncluster, outfile="")
doParallel::registerDoParallel(cl)
foreach (gRNAindex=1:dim(gRNAdm0)[1], .packages = c("MAST", "Matrix", "edgeR","gtools")) %dopar% {
#for (gRNAindex in 1:dim(gRNAdm0)[1]){
  gRNA <- rownames(gRNAdm0)[gRNAindex]
  print(gRNA)
  gRNAlocus <- strsplit(gRNA, split = "_")[[1]][1]
  gcount <- genedm[,gRNAdm0[gRNA,] >0 & colSums(gRNAdm0[rownames(gRNAdm0) != gRNA,]) ==0]
  if (dim(gcount)[2] !=0 & gRNAlocus != "neg") {
    # ncount1 <- genedm[, colnames(dm1dfagg)[dm1dfagg["neg",] >0 & nlocus==1]]
    # y <- DE_process(gcount, ncount1, filtcpm=filtcpm, filtpercent=filtpercent, perm=F)
    # res1 <- run_MASTcpmDetRate(y)
    # if (Nperm >0) {
    #   permres1 <- list()
    #   for (i in 1:Nperm) {
    #     y <- DE_process(gcount, ncount1, filtcpm=filtcpm, filtpercent=filtpercent, perm=T)
    #     permres1[[i]] <- run_MASTcpmDetRate(y)
    #   }
    #   save(res1, permres1, file=paste0("data/gRNA_MASTcpmdet/", gRNA, "_MASTcpmdet_Neg1.Rd"))
    # } else {
    #   save(res1, file=paste0("data/gRNA_MASTcpmdet/", gRNA, "_MASTcpmdet_Neg1.Rd"))
    # }

    ncount2 <- genedm[, colnames(dm1dfagg)[dm1dfagg[gRNAlocus,] == 0]] # much larger neg control, takes very long time
    y <- DE_process(gcount, ncount2, filtcpm=filtcpm, filtpercent=filtpercent, perm=F)
    res2 <- run_MASTcpmDetRate(y)
    if (Nperm >0) {
      permres2 <- list()
      for (i in 1:Nperm) {
        y <- DE_process(gcount, ncount2, filtcpm=filtcpm, filtpercent=filtpercent, perm=T)
        permres2[[i]] <- run_MASTcpmDetRate(y)
      }
      save(res2, permres2, file=paste0("data/gRNA_MASTcpmdet/", gRNA, "_MASTcpmdet_Neg2.Rd"))
    } else {
      save(res2, file=paste0("data/gRNA_MASTcpmdet/", gRNA, "_MASTcpmdet_Neg2.Rd"))
    }
  }
}
parallel::stopCluster(cl)



