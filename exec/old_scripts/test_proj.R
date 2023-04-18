library(SpaceTimePaper)
data(denguedata)

nt=7
ns = 9
n = ns*nt
Q_ICAR=matrix(0,nrow=9,ncol=9)
Q_ICAR[1,c(2,4)] = -1
Q_ICAR[2,c(1,3,5)] = -1
Q_ICAR[3,c(2,6)] = -1
Q_ICAR[4,c(1,5,7)] = -1
Q_ICAR[5,c(2,4,6,8)] = -1
Q_ICAR[6,c(3,5,9)] = -1
Q_ICAR[7,c(4,8)] = -1
Q_ICAR[8,c(5,7,9)] = -1
Q_ICAR[9,c(6,8)] = -1
for(i in 1:ns)
  Q_ICAR[i,i] = -sum(Q_ICAR[i,-i])
ns = ncol(Q_ICAR)
Q_RW2=GMRF_RW(n=nt,order=2)
Q_st=kronecker(Q_RW2,Q_ICAR)
Q.e = eigen(Q_st)
eT = rep(1,nt)
eT = matrix(eT/sqrt(sum(eT^2)),ncol=1)
eS = rep(1,ns)
eS = matrix(eS/sqrt(sum(eS^2)),ncol=1)
k1 = nt
k2 = ns-1
A1 = (Q.e$vectors[,(n-k1+1):(n)])
A2 = (Q.e$vectors[,(n-k1-k2+1):(n-k1)])
k3 = n-ncol(A1)-ncol(A2)
A3 = (Q.e$vectors[,1:k3])
Lambda3 = Q.e$values[1:k3]

eps = 1e-04
Q.eps = Q_st+diag(eps,n)
Q.eps.inv = solve(Q.eps)
CC = cbind(A1,A2)
QC = Q.eps.inv%*%(CC)
Sigma.C = Q.eps.inv-QC%*%solve(t(CC)%*%QC)%*%t(QC)
Sigma.C2 = (A3)%*%diag(1/(Lambda3+eps))%*%t(A3)
show(range(Sigma.C-Sigma.C2))

kappa = 1e05
Z = diag(n)-(A1)%*%t(A1)
Q.J = rbind(cbind(diag(kappa,n),-kappa*Z),
            cbind(-kappa*Z,kappa*Z+Q.eps))
Sigma.star = diag(1/kappa,n)+Z%*%Q.eps.inv%*%t(Z)
Sigma.star = diag(1/kappa,n)+(A2)%*%t(A2)/eps+(A3)%*%diag(1/(Lambda3+eps))%*%t(A3)
Sigma.star2 = solve(Q.J)[1:n,1:n]
show(range(Sigma.star-Sigma.star2))
Sigma.star2 = A3%*%diag(1/(Lambda3+eps))%*%t(A3)+(diag(1,n)-A2%*%t(A2))/kappa
plot(Sigma.C,Sigma.star2);abline(c(0,1))

if(0)
{
A = Q.J[1:n,1:n]
B = Q.J[1:n,-c(1:n)]
C = Q.J[-c(1:n),1:n]
D = Q.J[-c(1:n),-c(1:n)]
A11 = solve(A)+solve(A)%*%B%*%solve(D-C%*%solve(A)%*%B)%*%C%*%solve(A)
}