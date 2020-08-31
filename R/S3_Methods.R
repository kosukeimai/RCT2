print.CADE <- function(x) UseMethod("print.CADE")
print.CADE <- function(x){
  if(length(x) == 24){
    variance <- c(x$var.CADE1, x$var.CADE0, x$var.CASE1, x$var.CASE0,x$var.DEY1, x$var.DEY0, x$var.DED1, x$var.DED0, x$var.SEY1, x$var.SEY0, x$var.SED1, x$var.SED0)
    estimate <- c(x$CADE1, x$CADE0, x$CASE1, x$CASE0, x$DEY1, x$DEY0, x$DED1, x$DED0, x$SEY1, x$SEY0, x$SED1, x$SED0)
    name <- c("CADE1", "CADE0", "CASE1", "CASE0", "DEY1", "DEY0", "DED1", "DED0", "SEY1", "SEY0", "DED1", "DED0")
    out <- as.data.frame(cbind(name, estimate, variance))
  }else if(length(x) == 14){
    CIs <- c(x$CADE1.CI, x$CADE0.CI)
    CADEs <- c(x$CADE1, x$CADE2)
    name <- c("CADE1", "CADE0")
    variances <- c("Cluster Robust Variance", "HC2 Variance", "Cluster Robust HC2 Variance", "Individual Variance", "Proposed Variance")
    cluster_robust_variance <- c(x$var1.clu, x$var0.clu)
    HC2_variance <- c(x$var1.hc2, x$var0.hc2)
    cluster_robust_HC2_variance <- c(x$var1.clu.hc2, x$var0.clu.hc2)
    individual_variance <- c(x$var1.ind, x$var0.ind)
    proposed_variance <- c(x$var1.reg, x$var0.reg)
    variance_values <- as.data.frame(rbind(cluster_robust_variance, HC2_variance, cluster_robust_HC2_variance, individual_variance, proposed_variance))
    names(variance_values) <- c("var(CADE1)", "var(CADE0)")
    out <- list(as.data.frame(cbind(name, CADEs, CIs)), variance_values)
  }else if(length(x) == 12){
    names <- c("ITT Direct Effects", "IV Direct Effects", "ITT Spillover Effects", "IV Spillover Effects")
    effects <- rbind(x$ITT.DE, x$IV.DE, x$ITT.SE, x$IV.SE)
    CIs <- rbind(x$ITT.DE.CI, x$ITT.SE.CI, x$IV.DE.CI, x$IV.SE.CI)
    init <- as.data.frame(cbind(names, effects, CIs))
    colnames(init) <- c("Method", "Treatment", "Control", "Treatment CI", "Control CI")
    tstat_pvals <- as.data.frame(cbind(x$ITT.pvals, x$ITT.tstat, x$IV.pvals, x$IV.tstat))
    colnames(tstat_pvals) <- c("ITT pvalues", "ITT tstat", "IV pvalues", "IV tstat")
    out <- list(init, tstat_pvals)
    
  }
  return(out)
}
print.CADE(rand)
print.CADE(reg)
print.CADE(paramreg)
