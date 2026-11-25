#' Fitted model parameters
#'
#' Parameter sets from the final round of ABC-SMC model fitting, from which
#' DEN.spatial() samples one parameter set per simulation with probability
#' proportional to totaldevs.
#'
#' @format A data.frame with one row per parameter set and columns hometime,
#' mospdeath, betaEnv_mean, betaEnv_cor, Mov_model_type, pdetec, timedevs_R,
#' spacedevs_R and totaldevs (sampling weight).
#' @source {contact author for more detail on fitting}
"finalWeights"
