#' SpatioTemporalDataset
#'
#' Simulated data set using the routine simulate_data.R
#' Data are simulated fron the Knorr-Held type IV model with the spatial model defined 
#' through the spatial structure of the 54 regions in Germany and the temporal effect
#' following i RW(2) model.
#' The response is generated from a Poisson model with intercept equal to -1.
#' For each space-time combination, 30 observations are simulated.
#' 
#' @docType data
#' 
#' @format A data.table with 4 coloumns 
#' \describe{
#' \item{main_spatial}{Index for the spatial effect}
#' \item{main_temporal}{Index for the temporal effect}
#' \item{interaction}{Index for the spatio-temporal effect}
#' \item{main_spatial}{The responce variable}
#'  }
#' 
#' @keywords datasets
#' 
#' @source <https://github.com/folkehelseinstituttet/surveillance_data/tree/master/covid19>
#' 
"SpatioTemporalDataset"