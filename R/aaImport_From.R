
#' @importFrom dplyr %>%  all_of arrange bind_rows filter group_by left_join mutate pull select summarize
#' summarise
#' @importFrom ggplot2 aes element_blank element_text expand_limits facet_wrap facet_grid geom_boxplot geom_line
#' ggplot ggplotGrob geom_bar geom_hline geom_rect geom_ribbon geom_point guides guide_legend labs scale_x_date theme
#' theme_bw theme_classic xlim ylim xlab ylab
#' @importFrom graphics abline arrows axis axTicks barplot boxplot contour hist
#' @importFrom stats dnorm
NULL



#' Construct an empty `Data` object
#'
#'
#' @export
Make_Data <- function(silent=FALSE) {
  if (!silent) message('Creating an empty `Data` object')
  data <- list()
  data$Metadata <- data.frame(Field=c('Fishery',
                                      'Location',
                                      'Species',
                                      'Common_Name',
                                      'Author',
                                      'Date'),
                              Value=NA)
  data$Ages <- numeric()
  data$Weight_Age_Mean <- numeric()
  data$Weight_Age_SD <- numeric()
  data$Maturity_at_Age <- numeric()
  data$M_at_Age <- numeric()
  data$Post_Spawning_Mortality <- numeric()
  data$CAW <- numeric()
  data$CAW_ESS <- numeric()
  data$Weight_Bins <- numeric()
  data$Weight_Mids <- numeric()
  data$Year <- numeric()
  data$Month <- numeric()
  data$Effort_Mean <- numeric()
  data$Effort_SD <- numeric()
  data$Index_Mean <- numeric()
  data$Index_SD <- numeric()

  data$use_Frwpen <- 1
  data$use_R0rwpen <- 1
  data$model <- 'SLAM'
  class(data) <- 'Data'
  data
}

#' @rdname Make_Data
#' @export
Make <- function(type='Data') {
  if (type=='Data') {
    return(Make_Data())
  }
}

#' Create new Parameters list
#'
#' @return A list of class `Parameters`
#' @export
New_Parameters <- function() {

  parameters <- list()
  parameters$ls50 <- numeric()
  parameters$lsdelta <- numeric()
  parameters$logF_minit <- numeric()
  parameters$logF_ts <- numeric()

  parameters$log_sigmaF_m <- numeric()

  parameters$logR0_m_est <- rep(0, 11)
  parameters$log_sigmaR0 <- numeric()
  parameters$logRec_Devs <- numeric()
  parameters$log_sigmaR  <- numeric()
  class(parameters) <- 'Parameters'
  parameters
}

#' Check a `Data` object
#'
#' @param Data An object of class `Data`
#'
#' @return Nothing
#' @export
Check <- function(Data=NULL, silent=FALSE) {
  if (!inherits(Data, 'Data'))
    stop('Object must be class `Data`')

  # Add checks for missing data etc

  if (all(Data$CAW_ESS<1)) {
    stop('No catch-at-weight Data detected')
  }
  if (!silent)
    message('Data passes all checks')
  invisible(Data)

}

Update <- function(Data) {
  Data$Fit_Effort <- ifelse(sum(is.na(Data$Effort_Mean)) == length(Data$Year) | sum(!is.na(Data$Effort_Mean)) <2,
                            0,1)

  Data$Fit_Index <- ifelse(sum(is.na(Data$Index_Mean)) == length(Data$Year) | sum(!is.na(Data$Index_Mean)) <2,
                           0,1)
  Data$Fit_CAW <- 1
  Data$Fit_CAA <- 0
  maxage <- length(Data$Weight_Age_Mean)-1
  if (is.null(Data$Ages)) Data$Ages <- 0:maxage
  Data
}


check_numeric <- function(sl, Parameters) {
  if(!is.numeric(Parameters[[sl]]))
    stop(paste0('Parameters$', sl), ' is not a numeric value')
}

check_NA <- function(sl, Parameters) {
  if(any(is.na(Parameters[[sl]])))
    stop(paste0('Parameters$', sl), ' is not a numeric value')
}

check_length <- function(sl, len=1, Parameters) {
  if (!(length(Parameters[[sl]])==len))
    stop(paste0('Parameters$', sl), ' is not a length ', len)
}

#' Check Parameters
#'
#' @param Parameters A `Parameters` object
#' @param Data A `Data` object
#' @param silent Logical. Hide the messages?
#'
#' @return Nothing. Returns an error if there is a problem with the Parameters
#' @export
#'
Check_Parameters <- function(Parameters, Data, silent=FALSE) {

  if(!inherits(Parameters, 'Parameters'))
    stop('Object must be class `Parameters`')

  df <- data.frame(Name=names(Parameters),
                   Length=1)

  df$Length[match('logF_ts', df$Name)] <- length(Data$Year)
  df$Length[match('logRec_Devs', df$Name)] <- length(Data$Year)
  df$Length[match('logR0_m_est', df$Name)] <- 11
  df$Length[match('F_init_prior ', df$Name)] <- 2


  for (i in 1:nrow(df)) {
    check_numeric(df$Name[i], Parameters)
    check_NA(df$Name[i], Parameters)
    check_length(df$Name[i], df$Length[i], Parameters)
  }

  if (!silent)
    message('Parameters pass all checks')
}

