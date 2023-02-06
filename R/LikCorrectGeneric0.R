#' Calculated the correction factor needed for the generic0 model within INLA
#' 
#' @param 'Q' The precision matrix 
#' @param 'A' The constraint matrix
#' @param 'eps' The noise terms added to the diagonal
#' 
#' @return The correction factor (for the LOG-marginal likelohood)
#'
#' @export
LikCorrectGeneric0 = function(Q,A,eps)
{
  n = nrow(Q)
  C2 = forceSymmetric(as(Q+diag(n)*eps,"sparseMatrix"))
  QA = solve(C2,t(A))
  AQA = forceSymmetric(solve(t(QA)%*%t(A)))
  Sig.cond = forceSymmetric(solve(C2)-QA%*%AQA%*%t(QA))
  val = eigen(Sig.cond)$val
  np = diff(dim(A))
  logdet = sum(log(val[1:np])) 
  -0.5*logdet
}