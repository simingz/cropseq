empiricalFDR <- function(pvalues, realp, nullpall) {
  # will caculate FDR based on permutations
  .empiricalFDR <- function(p, realp, nullpall){
      realc <- length(realp[realp <= p])
      nullc <- length(nullpall[nullpall <= p]) /length(nullpall) * length(realp)
      fdr <- nullc/realc
      return(c(realc, nullc, fdr))
  }
  fdrmetric <- t(sapply(pvalues, .empiricalFDR, realp=realp, nullpall=nullpall))
  colnames(fdrmetric) <- c("#in_real", "#in_perm","EmpiricalFDR")
  return(fdrmetric)
}

empiricalPvalue <- function(pvalues, nullpall) {
  # will caculate FDR based on permutations
  .empiricalPvalue <- function(p, nullpall){
    length(nullpall[nullpall <= p]) /length(nullpall)
  }
  p.emp <- sapply(pvalues, .empiricalPvalue, nullpall=nullpall)
  return(p.emp)
}
