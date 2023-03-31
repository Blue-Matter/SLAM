
#' @importFrom dplyr %>%  all_of arrange bind_rows filter group_by left_join mutate pull select summarize
#' summarise
#' @importFrom ggplot2 aes element_blank element_text expand_limits facet_wrap facet_grid geom_boxplot geom_line
#' ggplot ggplotGrob geom_bar geom_rect geom_ribbon geom_point guides guide_legend labs scale_x_date theme
#' theme_bw theme_classic xlim ylim xlab ylab
#' @importFrom graphics abline arrows axis axTicks barplot boxplot contour hist
#' @importFrom stats dnorm
NULL


#' Class Constructor
#'
#' Put some details about it
#'
#' @param type an object
#'
#' @export
Create <- function(type='Data') {
  if (type=='Data') {
    return(Make_Data)
  }
}


Make_Data <- function() {
  data <- list()
  data$Metadata <- data.frame(Field=c('Fishery',
                                      'Location',
                                      'Species',
                                      'Common_Name',
                                      'Author',
                                      'Date'),
                                     Value=NA)
  data$Ages <- numeric()
  data$Weight_Age <- numeric()
  data$Weight_Age_SD <- numeric()
  data$Maturity_at_Age <- numeric()
  data$M_at_Age <- numeric()
  data$Post_Spawning_Mortality <- numeric()
  data$CAW <- numeric()
  data$CAW_ESS
  data$Weight_Bins
  data$Weight_Mids
  data$TS <- data.frame(Year=NA,
                        Month=NA,
                        Effort=NA,
                        Effort_SD=NA,
                        Index=NA,
                        Index_SD=NA)
  data$use_Frwpen <- 1
  data$use_R0rwpen <- 1
  data$model <- 'SLAM'
  class(data) <- 'Data'
  data
}

Make_Parameters <- function() {

}
