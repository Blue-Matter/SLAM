
#' Generate Data and Assessment Reports
#'
#' Generates a HTML report and opens it in the web browser
#'
#' @param x An object of class `Data` or `Assess`
#' @param filename The name of the output file. Defaults to 'Report.html'
#' @param dir The directory to save the report. Defaults to a temporary directory.
#' @param open_file Logical. Open the file in the browser?
#' @param silent Logical. Hide messages?
#' @param ... Other arguments (currently not used)
#'
#' @return Generates a HTML report
#' @export
Report <- function(x,
                   filename = "Report",
                   dir = tempdir(),
                   open_file = TRUE,
                   silent = FALSE, ...) {
  UseMethod('Report')
}

#' @export
#' @rdname Report
Report.Data <- function(x,
                        filename = "Data_Report",
                        dir = tempdir(),
                        open_file = TRUE,
                        silent = FALSE,...) {
  data <- list()
  data$Data <- Check(x, silent=TRUE)

  data$type <- 'Data Report'

  rmd_file <- file.path(system.file(package = "SLAM"), "Report_Data.Rmd")
  rmd <- readLines(rmd_file)

  write(rmd, file = file.path(dir, paste0(filename, ".rmd")))

  if(!silent) message("Creating ", data$type, ": ", file.path(dir, paste0(filename, ".html")))

  out <- rmarkdown::render(file.path(dir, paste0(filename, ".rmd")), "html_document", paste0(filename, ".html"), dir,
                           output_options = list(df_print = "paged"), quiet = TRUE)
  if(!silent) message("Rendering complete.")

  if(open_file) browseURL(out)
  return(invisible(out))
}

#' @export
#' @rdname Report
Report.Assess <- function(x,
                        filename = "Assessment_Report",
                        dir = tempdir(),
                        open_file = TRUE,
                        silent = FALSE,...) {
  data <- x
  data$type <- 'Assessment Report'

  rmd_file <- file.path(system.file(package = "SLAM"), "Report_Assessment.Rmd")
  rmd <- readLines(rmd_file)

  write(rmd, file = file.path(dir, paste0(filename, ".rmd")))

  if(!silent) message("Creating ", data$type, ": ", file.path(dir, paste0(filename, ".html")))

  out <- rmarkdown::render(file.path(dir, paste0(filename, ".rmd")), "html_document", paste0(filename, ".html"), dir,
                           output_options = list(df_print = "paged"), quiet = TRUE)
  if(!silent) message("Rendering complete.")

  if(open_file) browseURL(out)
  return(invisible(out))
}


report_Effort_TS <- function(data) {
  rep <- data$rep
  data <- data$Data
  df <- data.frame(Year=data$Year,
                   Month=data$Month,
                   Effort_Mean=data$Effort_Mean,
                   Effort_SD=data$Effort_SD,
                   Index_Mean=data$Index_Mean,
                   Index_SD=data$Index_SD
  )
  if (is.numeric(df$Month)) {
    df$Month <- month.abb[df$Month]
  }
  df$Month_n <- match(df$Month, month.abb)
  df$Date <- as.Date(paste(df$Year, df$Month_n, 01, sep='-'))

  # Effort
  mu <- log(df$Effort_Mean) -0.5*df$Effort_SD^2
  Lower <- qlnorm(0.1, mu,df$Effort_SD)
  Upper <- qlnorm(0.9, mu,df$Effort_SD)
  dfE <- df %>% select(Year,Month, Effort_Mean, Effort_SD, Date)
  dfE$ymin <- Lower
  dfE$ymax <- Upper

  if (!is.null(rep)) {
    dfE$predE <- rep$StEffort
    mu <- mean(dfE$Effort_Mean, na.rm=T)
    dfE$predE <- dfE$predE/mean(dfE$predE[!is.na(dfE$Effort_Mean)]) * mu
  }

  p1 <- ggplot(dfE, aes(x=Date, y=Effort_Mean, ymin=ymin, ymax=ymax)) +
    geom_ribbon(fill='lightgray') +
    geom_line() +
    labs(x='Date',
         y='Relative Effort',
         title='Index of Effort with 10th and 90th percentiles') +
    expand_limits(y=0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y",
                 limits = c(min(dfE$Date), max = max(dfE$Date)),
                 expand=c(0.01,0))


  if (!is.null(rep)) {
    p1 <- p1 + geom_line(aes(y=predE), color='blue', linetype=2)
  }
  p1
}

report_Index_TS <- function(data) {
  rep <- data$rep
  data <- data$Data
  df <- data.frame(Year=data$Year,
                   Month=data$Month,
                   Effort_Mean=data$Effort_Mean,
                   Effort_SD=data$Effort_SD,
                   Index_Mean=data$Index_Mean,
                   Index_SD=data$Index_SD
  )
  if (is.numeric(df$Month)) {
    df$Month <- month.abb[df$Month]
  }
  df$Month_n <- match(df$Month, month.abb)
  df$Date <- as.Date(paste(df$Year, df$Month_n, 01, sep='-'))

  # Index
  mu <- log(df$Index_Mean) -0.5*df$Index_SD^2
  Lower <- qlnorm(0.1, mu,df$Index_SD)
  Upper <- qlnorm(0.9, mu,df$Index_SD)
  dfI <- df %>% select(Year,Month, Index_Mean, Index_SD, Date)
  dfI$ymin <- Lower
  dfI$ymax <- Upper

  if (!is.null(rep)) {
    dfI$predI <- rep$stpredIndex
    mu <- mean(dfI$Index_Mean, na.rm=T)
    dfI$predE <- dfI$predI/mean(dfI$predI[!is.na(dfI$Index_Mean)]) * mu
  }

  p2 <- ggplot(dfI, aes(x=Date, y=Index_Mean, ymin=ymin, ymax=ymax)) +
    geom_ribbon(fill='lightgray') +
    geom_line() +
    labs(x='Date',
         y='Relative Abundance',
         title='Index of Abundance with 10th and 90th percentiles') +
    expand_limits(y=0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y",
                 limits = c(min(dfI$Date), max = max(dfI$Date)),
                 expand=c(0.01,0))

  if (!is.null(rep)) {
    p2 <- p2 + geom_line(aes(y=predI), color='blue', linetype=2)
  }

  p2
}



report_CAW <- function(data) {
  rep <- data$rep
  data <- data$Data

  df_list <- list()
  for (i in 1:length(data$Year)) {
    df_list[[i]] <- data.frame(Year=data$Year[i],
                               Month=data$Month[i],
                               Weight=data$Weight_Mids,
                               Count=data$CAW[,i])
  }
  df <- do.call('rbind', df_list)

  if (is.numeric(df$Month)) {
    df$Month <- month.abb[df$Month]
  }
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)

  if (!is.null(rep$predCAW)) {
    df$predCAW <- as.vector(rep$predCAW)

    df <- df %>% group_by(Year, Month) %>%
      mutate(n_sample=sum(Count)) %>%
      mutate(predCAW=predCAW*n_sample)
  }

  p <- ggplot(df, aes(x=Weight, y=Count)) +
    facet_grid(Month~Year) +
    geom_bar(stat='identity')  +
    expand_limits(y=0) +
    theme_bw()

  if (!is.null(rep$predCAW)) {
    p <- p + geom_line(aes(y=predCAW), color='blue', linetype=1)
  }

  list(df=df, p=p)

}

report_weight_at_age <- function(data) {

  # Weight-at-Age
  mu <- log(data$Weight_Age_Mean) -0.5*data$Weight_Age_SD^2
  Lower <- qlnorm(0.1, mu,data$Weight_Age_SD ) %>% round(2)
  Upper <- qlnorm(0.9, mu,data$Weight_Age_SD ) %>% round(2)
  df <- data.frame(Age=data$Ages, Mean_Weight=data$Weight_Age_Mean, Lower=Lower, Upper=Upper)

  p <- ggplot(df, aes(x=Age, y=Mean_Weight, ymin=Lower, ymax=Upper)) +
    geom_ribbon(fill='lightgray') +
    geom_line() +
    labs(x='Age (months)', y='Weight',
         title='Mean Weight-at-Age',
         subtitle = ' 10th and 90th percentiles') +
    expand_limits(y=0) +
    theme_bw()

  list(df=df, p=p)
}

report_maturity_at_age <- function(data) {
  df <- data.frame(Age=data$Ages, Maturity=data$Maturity_at_Age)

  p <- ggplot(df, aes(x=Age, y=Maturity)) +
    geom_line(linewidth=1.2) +
    labs(x='Age (months)', y='Probability Mature',
         title='Maturity-at-Age') +
    expand_limits(y=0) +
    theme_bw()
  list(df=df, p=p)
}

report_M_at_age <- function(data) {
  df <- data.frame(Age=data$Ages, M=data$M_at_Age)

  p <- ggplot(df, aes(x=Age, y=M)) +
    geom_line(linewidth=1.2) +
    labs(x='Age (months)', y='Natural Mortality',
         title='Natural Mortality-at-Age') +
    expand_limits(y=0) +
    theme_bw()

  list(df=df, p=p)
}

report_PSM_at_age <- function(data) {
  df <- data.frame(Age=data$Ages, PSM=data$Post_Spawning_Mortality)
  p <- ggplot(df, aes(x=Age, y=PSM)) +
    geom_line(linewidth=1.2) +
    labs(x='Age (months)', y='Post-Spawning Mortality',
         title='Post-Spawning Mortality-at-Age',
         linetype='') +
    theme_bw() +
    expand_limits(y=0)

  list(df=df, p=p)
}

report_seasonal_rec <- function(data) {
  rep <- data$rep
  df <- data.frame(Month=month.abb, Rec=rep$R0_m)
  df$Month <- factor(df$Month, levels=month.abb, ordered = TRUE)
  p1 <- ggplot(df, aes(x=Month, y=Rec, group=1)) +
    expand_limits(y=0) +
    geom_line() +
    theme_bw() +
    labs(y='Relative Recruitment')

  list(df=df, p=p1)
}

assumed_pars <- function(data) {

  df <- data.frame(Parameter=c('h',
                               'sigmaR',
                               'use_F_pen',
                               'sigmaF_m',
                               'use_R0_pen',
                               'sigmaR0'),
                   Value=c(data$Data$h,
                           exp(data$Parameters$log_sigmaR),
                           as.logical(data$Data$use_Frwpen),
                           exp(data$Parameters$log_sigmaF_m),
                           as.logical(data$Data$use_R0rwpen),
                           exp(data$Parameters$log_sigmaR0)),
                   Description=c('Assumed steepness of stock-recruit curve',
                                 'Standard deviation of log-normal recruitment deviations',
                                 'Use the penalty for random walk in F?',
                                 'Standard deviation for random walk penalty for F',
                                 'Use the penalty for random walk in seasonal recruitment?',
                                 'Standard deviation for random walk penalty for seasonal recruitment')
  )
  df
}

nll_vals <- function(data) {
  data.frame(Component=c('Catch-at-Weight',
                        'Index of Effort',
                        'Index of Abundance',
                        'Recruitment deviations',
                        'Random Walk F',
                        'Random Walk Seasonal Recruitment'),
             Value=data$rep$nll_joint)

}
