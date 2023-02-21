#include <TMB.hpp>
#include <math.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace Eigen;

  //Load data--------------
  DATA_VECTOR(y); //The response
  DATA_SPARSE_MATRIX(Q_t);
  DATA_SPARSE_MATRIX(Q_s);
  DATA_SPARSE_MATRIX(Q_st);
  DATA_SPARSE_MATRIX(Z_s);
  DATA_SPARSE_MATRIX(Z_t);
  DATA_SPARSE_MATRIX(Z_st);
  DATA_MATRIX(X);
  DATA_INTEGER(m_s);
  DATA_INTEGER(m_t);
  DATA_INTEGER(m_st);
  DATA_INTEGER(n_data);
  DATA_SCALAR(shape);
  DATA_SCALAR(scale);
  DATA_SCALAR(fixed_sd);//Number of data points
  //-----------------------

  //lOAD PARAMETERS--------
  PARAMETER_VECTOR(xs); //Regression coefficients
  PARAMETER_VECTOR(xt);
  PARAMETER_VECTOR(xst);//Regression coefficients
  PARAMETER(log_tau_s);
  PARAMETER(log_tau_st);//Parameter in the
  PARAMETER(log_tau_t);
  PARAMETER_VECTOR(beta_fixed);//Parameter in the AR1 covariance function
  //The spatio-temporal latent field
  //Parameter for unexplained variation
  //------------------------

  //Transform variables-----
  Type nll = 0;

  Type tau_s = exp(log_tau_s);
  Type tau_t = exp(log_tau_t);
  Type tau_st = exp(log_tau_st);

  nll-=dgamma(tau_s,shape,scale);
  nll-=dgamma(tau_st,shape,scale);
  nll-=dgamma(tau_t,shape,scale);

  //Add priors:
  //------------------------

  //Construct sparce precision matrix for latent field---
  SparseMatrix<Type> Q_s_tau =Q_s*tau_s;
  SparseMatrix<Type> Q_t_tau =Q_t*tau_t;
  SparseMatrix<Type> Q_st_tau =Q_st*tau_st;


  //------------------------

  //Calculates nll-------------------------------
  vector<Type> p_xs = Z_s*xs;
  vector<Type> p_xt = Z_t*xt;
  vector<Type> p_xst = Z_st*xst;
  vector<Type> p_x=X*beta_fixed;
  for(int k=0; (k<beta_fixed.size());k++){
    nll-=dnorm(beta_fixed(k),Type(0),fixed_sd,true);
  }

  //Add prior
    nll -= Type(0.5)*Type(m_s)*log_tau_s -Type(0.5)*log_tau_s*GMRF(Q_s).Quadform(xs);
    nll -= Type(0.5)*Type(m_st)*log_tau_st -Type(0.5)*log_tau_st*GMRF(Q_st).Quadform(xst);
    nll -= Type(0.5)*Type(m_t)*log_tau_t -Type(0.5)*log_tau_t*GMRF(Q_t).Quadform(xt);

  //  nll = SEPARABLE(AR1(rho), GMRF(Q_s))(x); //Alternative approach with use of seperability, includes normalizing constant here because TMB seems to be both faster and more stable when including it.


  vector<Type> mu(n_data);
    for(int j =0; (j<n_data); ++j){

        mu(j) = p_x(j)+p_xs(j)+p_xst(j)+p_xt(j);
        nll -= dpois(y(j), exp(mu(j)),true);
  }
  //---------------------------------------------
  ADREPORT(mu);
  //Report what we want to report----------------
  //---------------------------------------------

  return nll;
}
