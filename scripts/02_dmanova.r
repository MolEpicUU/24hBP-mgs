calculateK <- function (G, H, df2) {
  n <- nrow(G)
  dG <- diag(G)
  
  mu1 <- sum(dG) / (n - df2)
  mu2 <- (sum(G^2) - sum(dG^2)) / ((n - df2)^2 + sum(H^4) - 2 * sum(diag(H^2)))
  
  k <- mu1^2 / mu2
  list(K = k, mu1 = mu1, mu2 = mu2)
}

dmanova <- function (formula, data = NULL, positify = FALSE,
                     contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
                     returnG = FALSE) {
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  
  # The first group includes the intercept, nterms count the intercept
  u.grps <- unique(grps)
  nterms <- length(u.grps)
  
  Z <- rhs[, grps %in% u.grps[1:(nterms-1)], drop=F]
  XZ <- rhs[, grps %in% u.grps[1:(nterms)], drop=F]
  
  n <- nrow(XZ)
  if (inherits(lhs, 'dist')) {
    D <- as.matrix(lhs)
  } else {
    D <- lhs
  }
  
  D <- -D^2 / 2
  
  G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol=1) %*% colMeans(D)
  
  if (positify) {
    G <- as.matrix(nearPD(G)$mat)
  }
  
  XZi <- solve(t(XZ) %*% XZ)
  Zi <- solve(t(Z) %*% (Z))
  HZ <- Z %*% Zi %*% t(Z)
  HXZ <- XZ %*% XZi %*% t(XZ)
  HX <- HXZ - HZ
  HIXZ <- diag(n) - HXZ
  HIX <- diag(n) - HX
  
  
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)
  
  MSS <- sum(G * HX)
  RSS <- sum(G * HIXZ)
  TSS <- sum(diag(G))
  
  
  f.stat <- (MSS / df1) / (RSS / df2)
  
  GXZ <- G %*% XZ
  XZXZi <- XZ %*% XZi
  GXZtXZXZi <- GXZ %*% t(XZXZi)
  
  G.tilde <- G + XZXZi %*% (t(GXZ) %*% XZ) %*% t(XZXZi) - GXZtXZXZi - t(GXZtXZXZi)
  
  
  obj <- calculateK(G.tilde, HIXZ, n - df2)
  K <- obj$K
  
  
  p.value <- pchisq(f.stat * K * df1, df = K * df1, lower.tail = F)
  
  
  
  SumsOfSqs <- c(MSS, RSS, TSS)
  
  tab <- data.frame(Df = c(df1, df2, n - 1), SumsOfSqs = SumsOfSqs, 
                    MeanSqs = c(MSS/df1, RSS/df2, NA), F.Model = c(f.stat, 
                                                                   NA, NA), R2 = c(MSS/TSS, NA, NA), 
                    "Pr(>F)" = c(p.value, NA, NA))
  
  rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- c("Pr(>F)")
  class(tab) <- c("anova", class(tab))
  attr(tab, "heading") <- c("F stat and P value of the last term is adjusted by preceding terms!\n")
  
  if (returnG == TRUE) {
    out <- list(aov.tab = tab,  df = K, G = G,  call = match.call())
  } else {
    out <- list(aov.tab = tab, df = K, call = match.call())
  }
  
  out
}

dmanova.overall <- function (formula, data = NULL, positify = FALSE,
                             contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
                             returnG = FALSE) {
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  
  # The first group includes the intercept, nterms count the intercept
  u.grps <- unique(grps)
  nterms <- length(u.grps)
  
  ## this is changed
  Z <- rhs[, 1, drop=F] # was: Z <- rhs[, grps %in% u.grps[1:(nterms-1)], drop=F]
  XZ <- rhs[, , drop=F] # was: XZ <- rhs[, grps %in% u.grps[1:(nterms)], drop=F]
  
  n <- nrow(XZ)
  if (inherits(lhs, 'dist')) {
    D <- as.matrix(lhs)
  } else {
    D <- lhs
  }
  
  D <- -D^2 / 2
  
  G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol=1) %*% colMeans(D)
  
  if (positify) {
    G <- as.matrix(nearPD(G)$mat)
  }
  
  XZi <- solve(t(XZ) %*% XZ)
  Zi <- solve(t(Z) %*% (Z))
  HZ <- Z %*% Zi %*% t(Z)
  HXZ <- XZ %*% XZi %*% t(XZ)
  HX <- HXZ - HZ
  HIXZ <- diag(n) - HXZ
  HIX <- diag(n) - HX
  
  
  df1 <- ncol(XZ) - ncol(Z)
  df2 <- n - ncol(XZ)
  
  MSS <- sum(G * HX)
  RSS <- sum(G * HIXZ)
  TSS <- sum(diag(G))
  
  
  f.stat <- (MSS / df1) / (RSS / df2)
  
  GXZ <- G %*% XZ
  XZXZi <- XZ %*% XZi
  GXZtXZXZi <- GXZ %*% t(XZXZi)
  
  G.tilde <- G + XZXZi %*% (t(GXZ) %*% XZ) %*% t(XZXZi) - GXZtXZXZi - t(GXZtXZXZi)
  
  
  obj <- calculateK(G.tilde, HIXZ, n - df2)
  K <- obj$K
  
  
  p.value <- pchisq(f.stat * K * df1, df = K * df1, lower.tail = F)
  
  
  
  SumsOfSqs <- c(MSS, RSS, TSS)
  
  tab <- data.frame(Df = c(df1, df2, n - 1), SumsOfSqs = SumsOfSqs, 
                    MeanSqs = c(MSS/df1, RSS/df2, NA), F.Model = c(f.stat, 
                                                                   NA, NA), R2 = c(MSS/TSS, NA, NA), 
                    "Pr(>F)" = c(p.value, NA, NA))
  
  rownames(tab) <- c("Model(Adjusted)", "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- c("Pr(>F)")
  class(tab) <- c("anova", class(tab))
  attr(tab, "heading") <- c("F stat and P value of the last term is adjusted by preceding terms!\n")
  
  if (returnG == TRUE) {
    out <- list(aov.tab = tab,  df = K, G = G,  call = match.call())
  } else {
    out <- list(aov.tab = tab, df = K, call = match.call())
  }
  
  out
}
