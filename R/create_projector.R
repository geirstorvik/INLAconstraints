
#'Create the projector matrix fast
#'This is even fast for grid sizes of 40 000!
#'
#'@export
create_projector_matrix_interactionIV<-function(ns,nt){
  TempPart=as(Matrix::Diagonal(x=0,ns*nt),"TsparseMatrix")
  TempPart=as(TempPart,"TsparseMatrix")
  for(i in seq(1,ns*nt,by=ns)){
    TempPart[cbind(seq(i,ns*nt),seq(1,ns*nt-i+1))]=1/nt
    TempPart[cbind(seq(1,ns*nt-i+1),seq(i,ns*nt))]=1/nt

  }
  return(Matrix::Diagonal(x=1,ns*nt)-TempPart)
}
