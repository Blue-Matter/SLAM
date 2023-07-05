
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
                        silent = FALSE,
                        dev=FALSE,
                        ...) {
  data <- list()
  data$Data <- Check(x, silent=TRUE)

  data$type <- 'Data Report'

  if (!dev) {
    rmd_file <- file.path(system.file(package = "SLAM"), "Report_Data.Rmd")
  } else {
    rmd_file <- "inst/Report_Data.Rmd"
  }
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
                        dev=FALSE,
                        silent = FALSE,...) {
  data <- x

  if (!dev) {
    rmd_file <- file.path(system.file(package = "SLAM"), "Report_Assessment.Rmd")
  } else {
    rmd_file <- "inst/Report_Assessment.Rmd"
  }

  rmd <- readLines(rmd_file)

  write(rmd, file = file.path(dir, paste0(filename, ".rmd")))

  if(!silent) message("Creating ", 'Assessment Report', ": ", file.path(dir, paste0(filename, ".html")))

  out <- rmarkdown::render(file.path(dir, paste0(filename, ".rmd")), "html_document", paste0(filename, ".html"), dir,
                           output_options = list(df_print = "paged"), quiet = TRUE)
  if(!silent) message("Rendering complete.")

  if(open_file) browseURL(out)
  return(invisible(out))
}

# ---- Plot Input Parameters ----


#' Plot the at-age schedules
#'
#' @param data An object of class `Data`
#'
#' @return A `ggplot` object
#' @export
#'
#' @examples
#' Data <- Import(Example_Data()[1])
#' plot_at_Age(Data)
plot_at_Age <- function(data) {

  p1 <- plot_weight_at_age(data)$p + theme(axis.text.x=element_blank(),
                                           axis.title.x=element_blank())
  p2 <- plot_M_at_age(data)$p + theme(axis.text.x=element_blank(),
                                      axis.title.x=element_blank())
  p3 <- plot_maturity_at_age(data)$p
  p4 <- plot_PSM_at_age(data)$p
  cowplot::plot_grid(p1,p2,p3,p4, labels=paste0(letters[1:4], ')'))

}

#' @describeIn plot_at_Age Plot mean weight-at-age with 10th and 90th percentiles
#' @export
plot_weight_at_age <- function(data) {

  # Weight-at-Age
  mu <- log(data$Weight_Age_Mean) -0.5*data$Weight_Age_SD^2
  Lower <- qlnorm(0.1, mu,data$Weight_Age_SD ) %>% round(2)
  Upper <- qlnorm(0.9, mu,data$Weight_Age_SD ) %>% round(2)
  df <- data.frame(Age=data$Ages, Mean_Weight=data$Weight_Age_Mean, Lower=Lower, Upper=Upper)

  p <- ggplot(df, aes(x=Age, y=Mean_Weight, ymin=Lower, ymax=Upper)) +
    geom_ribbon(fill='lightgray') +
    geom_line() +
    labs(x='Age (months)', y='Weight') +
    expand_limits(y=0) +
    theme_bw()

  list(df=df, p=p)
}

#' @describeIn plot_at_Age Plot maturity-at-age schedule
#' @export
plot_maturity_at_age <- function(data) {
  df <- data.frame(Age=data$Ages, Maturity=data$Maturity_at_Age)

  p <- ggplot(df, aes(x=Age, y=Maturity)) +
    geom_line(linewidth=1.2) +
    labs(x='Age (months)', y='Probability Mature') +
    expand_limits(y=0) +
    theme_bw()
  list(df=df, p=p)
}

#' @describeIn plot_at_Age Plot natural mortality-at-age schedule
#' @export
plot_M_at_age <- function(data) {
  df <- data.frame(Age=data$Ages, M=data$M_at_Age)

  p <- ggplot(df, aes(x=Age, y=M)) +
    geom_line(linewidth=1.2) +
    labs(x='Age (months)', y='Natural Mortality') +
    expand_limits(y=0) +
    theme_bw()

  list(df=df, p=p)
}
#' @describeIn plot_at_Age Plot post-spawning mortality-at-age schedule
#' @export
plot_PSM_at_age <- function(data) {
  df <- data.frame(Age=data$Ages, PSM=data$Post_Spawning_Mortality)
  p <- ggplot(df, aes(x=Age, y=PSM)) +
    geom_line(linewidth=1.2) +
    labs(x='Age (months)', y='Post-Spawning Mortality',
         linetype='') +
    theme_bw() +
    expand_limits(y=0)

  list(df=df, p=p)
}










# ---- Plot Data (and fits) ----

#' Plot the weight composition data, and model fits to the data
#'
#' @param data An object of class `Data` or class `Assess`
#'
#' @return A `list` with a `data.frame` and  a `ggplot` object
#' @export
#'
#' @examples
#' Data <- Import(Example_Data()[1])
#' plot_CAW(Data)
plot_CAW <- function(data) {

  if (inherits(data, 'Assess')) {
    rep <- data$rep
    data <- data$Data

  }

  df_list <- list()
  for (i in 1:length(data$Year)) {
    df_list[[i]] <- data.frame(Year=data$Year[i],
                               Month=data$Month[i],
                               Weight=data$Weight_Mids,
                               Count=data$CAW[i,])
  }

  df <- do.call('rbind', df_list)

  df2 <- data.frame(Year=data$Year,
                    Month=data$Month,
                    Nsamp=round(data$CAW_ESS,0))



  df <- left_join(df, df2,by = join_by(Year, Month))

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
    geom_text(aes(x=-Inf, y=Inf, label=Nsamp, hjust = 0, vjust = 1)) +
    labs(x='Weight Class', y='Frequency') +
    theme_bw()

  if (!is.null(rep$predCAW)) {
    p <- p + geom_line(aes(y=predCAW), color='blue', linetype=1)
  }

  list(df=df, p=p)

}

#' Plot indices of relative biomass, effort, and catch
#'
#' @param data An object of class `Data` or class `Assess`
#'
#' @return A `ggplot object`
#' @export
#'
#' @examples
#' Data <- Import(Example_Data()[2])
#' plot_Indices(Data)
plot_Indices <- function(data) {
  p1 <- plot_Biomass(data) + theme(axis.title.x=element_blank(),
                                   axis.text.x=element_blank())
  p2 <- plot_Effort(data) + theme(axis.title.x=element_blank(),
                                  axis.text.x=element_blank())
  p3 <- plot_Catch(data)
  suppressWarnings(cowplot::plot_grid(p1, p2,p3, labels=paste(letters[1:3], ')'), nrow=3))
}

#' @describeIn plot_Indices Plot Index of Biomass time-series
#' @export
plot_Biomass <- function(data) {
  plot_TS(data, 'Index')
}

#' @describeIn plot_Indices Plot Index of Effort time-series
#' @export
plot_Effort <- function(data) {
  plot_TS(data, 'Effort')
}

#' @describeIn plot_Indices Plot Catch of Effort time-series
#' @export
plot_Catch <- function(data) {
  plot_TS(data, 'Catch')
}




plot_TS <- function(data, type='Catch') {

  if (inherits(data, 'Assess')) {
    rep <- data$rep
    data <- data$Data
  }

  nms <- paste(type, c('Mean', 'SD'), sep='_')
  df <- data.frame(Year=data$Year,
                   Month=data$Month)
  df$Mean <- data[[nms[[1]]]]
  df$SD <- data[[nms[[2]]]]

  if (is.null(df$Mean)) {
    return(NULL)
  }
  if (is.numeric(df$Month)) {
    df$Month <- month.abb[df$Month]
  }
  df$Month_n <- match(df$Month, month.abb)
  df$Date <- as.Date(paste(df$Year, df$Month_n, 01, sep='-'))

  # 10th and 90th percentiles
  mu <- log(df$Mean) -0.5*df$SD^2
  Lower <- qlnorm(0.1, mu,df$SD)
  Upper <- qlnorm(0.9, mu,df$SD)
  df$ymin <- Lower
  df$ymax <- Upper

  if (!is.null(rep)) {
    if (type=='Catch') {
      df$pred <- rep$predCB
      mu <- mean(df$Mean, na.rm=T)
      df$pred <- df$pred/mean(df$pred[!is.na(df$Mean)]) * mu
    }
    if (type=='Effort') {
      df$pred <- rep$StEffort
      mu <- mean(df$Mean, na.rm=T)
      df$pred <- df$pred/mean(df$pred[!is.na(df$Mean)]) * mu
    }
    if (type=='Index') {
      df$pred <- rep$stpredIndex
      mu <- mean(df$Mean, na.rm=T)
      df$pred <- df$pred/mean(df$pred[!is.na(df$Mean)]) * mu
    }

  }

  if (type=='Index') type <- 'Vulnerable Biomass'
  p1 <- ggplot(df, aes(x=Date, y=Mean, ymin=ymin, ymax=ymax)) +
    geom_ribbon(fill='lightgray') +
    geom_line() +
    labs(x='Date',
         y=paste('Index of ', type)) +
    expand_limits(y=0) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y",
                 limits = c(min(df$Date), max = max(dfE$Date)),
                 expand=c(0.01,0))


  if (!is.null(rep)) {
    p1 <- p1 + geom_line(aes(y=pred), color='blue', linetype=2)
  }
  p1

}


#  ---- Plot Model Predictions ----

#' Plot the predicted seasonal recruitment pattern
#'
#' @param assess An object of class `Assess`
#'
#' @return A list with a data.frame and a ggplot object
#' @export
#'
#' @examples
#' MyData <- Import(Example_Data()[2])
#' MyAssess <- Assess(MyData)
#' plot_seasonal_rec(MyAssess)
plot_seasonal_rec <- function(assess) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  df <- get_seasonal_recruitment(assess)
  p1 <- ggplot(df, aes(x=Month_Name, y=Recruitment, group=1)) +
    expand_limits(y=0) +
    geom_line() +
    theme_bw() +
    labs(y='Relative Recruitment')

  list(df=df, p=p1)
}

#' @describeIn plot_seasonal_rec get data.frame of predicted seasonal recruitment
#' @export
get_seasonal_recruitment <- function(assess) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  months <- assess$Data$Month[1:12]
  df <- data.frame(Month=months, Recruitment=assess$rep$R0_m) %>%
    arrange(Month)
  df$Month_Name <- month.abb[df$Month]
  df$Month_Name <- factor(df$Month_Name, levels=month.abb, ordered = TRUE)
  df

}


#' Plot the predicted selectivity-at-age schedule
#'
#' @param assess An object of class `Assess`
#'
#' @return A list with a data.frame and a ggplot object
#' @export
#'
#' @examples
#' MyData <- Import(Example_Data()[2])
#' MyAssess <- Assess(MyData)
#' plot_selectivity(MyAssess)
plot_selectivity <- function(assess) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  df <- data.frame(Age=assess$Data$Ages, Selectivity=assess$rep$selA)

  p1 <- ggplot(df, aes(x=Age, y=Selectivity)) +
    geom_line() +
    theme_bw()

  list(p=p1, df=df)
}


#' Plot the predicted spawning potential ratio (SPR)
#'
#' @param assess An object of class `Assess`
#'
#' @return A list with a data.frame and a ggplot object
#' @export
#'
#' @examples
#' MyData <- Import(Example_Data()[2])
#' MyAssess <- Assess(MyData)
#' plot_SPR(MyAssess)
plot_SPR <- function(assess) {

  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  df <- make_df(assess)
  ggplot(df, aes(x=Date, y=SPR)) +
    geom_line() +
    labs(x='Date',
         y='Spawning Potential Ratio (SPR)') +
    expand_limits(y=c(0,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y",
                 limits = c(min(df$Date), max = max(df$Date)),
                 expand=c(0.01,0))

}

#' Plot the predicted fishing mortality (F)
#'
#' @param assess An object of class `Assess`
#'
#' @return A list with a data.frame and a ggplot object
#' @export
#'
#' @examples
#' MyData <- Import(Example_Data()[2])
#' MyAssess <- Assess(MyData)
#' plot_F(MyAssess)
plot_F <- function(assess) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')


  df <- make_df(assess)
  ggplot(df, aes(x=Date, y=F)) +
    geom_line() +
    labs(x='Date',
         y='Fishing Mortality') +
    expand_limits(y=c(0)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y",
                 limits = c(min(df$Date), max = max(df$Date)),
                 expand=c(0.01,0))
}



report_rec_devs <- function(data) {
  df <- make_df(data)
  ggplot(df, aes(x=Date, y=log_rec_devs )) +
    geom_line() +
    labs(x='Date',
         y='Log recruitment deviations') +
    expand_limits(y=c(0)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_x_date(date_breaks = '1 month', date_labels = "%b-%Y",
                 limits = c(min(df$Date), max = max(df$Date)),
                 expand=c(0.01,0)) +
    geom_hline(yintercept = 0, linetype=2)
}



#' Report the negative log-likelihoods
#'
#' @param assess An object of class `Assess`
#'
#' @return A data.frame
#' @export
#'
#' @examples
#' MyData <- Import(Example_Data()[2])
#' MyAssess <- Assess(MyData)
#' report_NLL(MyAssess)
report_NLL <- function(assess) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  df <- data.frame(Component=c('Catch-at-Weight',
                               'Index of Effort',
                               'Index of Abundance',
                               'Index of Catch',
                               'Recruitment deviations',
                               'Random Walk F',
                               'Random Walk Seasonal Recruitment',
                               'Total'),
             Value=c(assess$rep$nll_joint, sum(assess$rep$nll_joint)))

  df$Value <- round(df$Value, 3)
  df

}



make_df <- function(data) {

  df <- data.frame(Year=data$Data$Year,
                   Month=data$Data$Month,
                   SPR=data$rep$SPR,
                   F=data$rep$F_m,
                   log_rec_devs=data$rep$logRec_Devs
  )

  if (is.numeric(df$Month)) {
    df$Month <- month.abb[df$Month]
  }
  df$Month_n <- match(df$Month, month.abb)
  df$Date <- as.Date(paste(df$Year, df$Month_n, 01, sep='-'))
  df
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


## ---- Reference Points -----

## Calculate and Plot F and SPR relative to Reference Points

#' @describeIn Calc_RP calculate average monthly F
#' @export
calc_average_monthly_F <- function(assess) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  months <- assess$Data$Month
  df <- data.frame(Month=months, F=assess$rep$F_m) %>%
    group_by(Month) %>% summarise(Mean=mean(F))
  df$Month_Name <- month.abb[df$Month]
  df$Month_Name <- factor(df$Month_Name, levels=month.abb, ordered = TRUE)
  df
}

#' @describeIn Calc_RP calculate average monthly SPR
#' @export
calc_average_monthly_SPR <- function(assess) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  months <- assess$Data$Month
  df <- data.frame(Month=months, SPR=assess$rep$SPR) %>%
    group_by(Month) %>% summarise(Mean=mean(SPR))
  df$Month_Name <- month.abb[df$Month]
  df$Month_Name <- factor(df$Month_Name, levels=month.abb, ordered = TRUE)
  df
}


#' Calculate SPR and F reference points from HARA model
#'
#' @param assess An object of class `Assess`
#' @param utilpow numeric vector of values for the HARA utility exponent
#'
#' @return A data.frame
#' @export

Calc_RP <- function(assess, utilpow=c(0.4, 0.9)) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')

  average_monthly_F <- calc_average_monthly_F(assess)
  average_monthly_F$Variable <- "F"

  average_monthly_SPR <- calc_average_monthly_SPR(assess)
  average_monthly_SPR$Variable <- "SPR"

  season_R0 <- get_seasonal_recruitment(assess)

  F_list <- list()
  SPR_list <- list()
  for (i in seq_along(utilpow)) {
    tt <- calculate_optimal_fishing(R0_m=season_R0$Recruitment,
                                              steepness=assess$Data$h,
                                              Weight_Age_Mean=assess$Data$Weight_Age_Mean,
                                              Maturity_at_Age=assess$Data$Maturity_at_Age,
                                              M_at_Age=assess$Data$M_at_Age,
                                              Post_Spawning_Mortality=assess$Data$Post_Spawning_Mortality,
                                              sel_at_age=assess$rep$selA,
                                              opt_type=1, utilpow=utilpow[i])


    average_monthly_F$Optimal <- tt$F_m
    average_monthly_F$util_pow <- utilpow[i]
    F_list[[i]] <- average_monthly_F %>%
      tidyr::pivot_longer(., cols=c(2, 5), names_to = 'Name', values_to = 'Value')


    average_monthly_SPR$Optimal <- tt$SPR
    average_monthly_SPR$util_pow <- utilpow[i]
    SPR_list[[i]] <- average_monthly_SPR %>%
      tidyr::pivot_longer(., cols=c(2, 5), names_to = 'Name', values_to = 'Value')

  }


  df <- bind_rows(do.call('rbind', F_list), do.call('rbind', SPR_list))
  df$util_pow <- factor(df$util_pow)
  df
}

#' @describeIn Calc_RP plot the mean F and mean SPR and the optimal values from the HARA model
#' @export
plot_RP <- function(assess, utilpow=c(0.4, 0.9)) {
  if(!inherits(assess, 'Assess'))
    stop('Must be object of class `Assess`')
 df <-  Calc_RP(assess, utilpow)
 ggplot(df, aes(x=Month_Name, y=Value, color=Name, group=Name)) +
   facet_grid(util_pow~Variable) +
   geom_line() +
   expand_limits(y=c(0,1)) +
   theme_bw() +
   labs(x='')
}





