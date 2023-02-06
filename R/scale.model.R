scale.model = function (Q, constr = NULL, eps = sqrt(.Machine$double.eps)) 
{
  marg.var <- rep(0, nrow(Q))
  Q <- inla.as.sparse(Q)
  g <- inla.read.graph(Q)
  for (k in seq_len(g$cc$n)) {
    i <- g$cc$nodes[[k]]
    n <- length(i)
    QQ <- Q[i, i, drop = FALSE]
    if (n == 1) {
      QQ[1, 1] <- 1
      marg.var[i] <- 1
    }
    else {
      cconstr <- constr
      if (!is.null(constr)) {
        cconstr$A <- constr$A[, i, drop = FALSE]
        eeps <- eps
      }
      else {
        eeps <- 0
      }
      res <- inla.qinv(QQ + Diagonal(n) *  
                         eeps, constr = cconstr)
      fac <- exp(mean(log(diag(res))))
      QQ <- fac * QQ
      marg.var[i] <- diag(res)/fac
    }
    Q[i, i] <- QQ
  }
  return(list(Q = Q, var = marg.var))
}
