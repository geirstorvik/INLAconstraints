#'This function fits the model as done by Rachel Lowe (2021).
#'NB: INLA mode experimental is used.
#'@param formula is the formula in INLA format
#'@param df is the data to be used in the model fitting
#'@param family is the distribution family fitted, e.g Poisson, Gaussian
#'@return
#'@export
#'@name mymodel
mymodel <- function(formula, data = df, family = "nbinomial", config = FALSE,num.threads=num.threads)

{
  model <- inla(formula = formula, data = data, family = family,num.threads=num.threads,control.fixed = list(
                                     prec.intercept =1),
                control.predictor = list(compute = TRUE),control.inla = list(cmin = 0 ),
                verbose = TRUE, inla.mode="experimental")
#rerun(model)
  return(model)
}
