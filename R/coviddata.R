#' coviddata
#'
#'  Test data for Covid-18 extracted form the Norwegian institute of public health
#' The data have been extracted by the script Make_data_covid_test.R. 
#' See the scripts Covid_GCF.R and Covid_SC.R on how to use these data
#' 
#' @docType data
#' 
#' @format A list with two components. 
#' \describe{
#' \item{data}{A data frame with four columns, date, location_code, cases and pop}
#' \item{adj}{A 11x11 matrix giving the adjency matrix between the 11 counties in Norway}
#' }
#' 
#' @keywords datasets
#' 
#' @source <https://github.com/folkehelseinstituttet/surveillance_data/tree/master/covid19>
#' 
"coviddata"