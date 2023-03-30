
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
#' @return
#' @export
#'
Report <- function(x,
                   filename = "Report",
                   dir = tempdir(),
                   open_file = TRUE,
                   silent = FALSE, ...) {
  UseMethod('Report', x)
}

#' @describeIn Report Generate a Data report
Report.Data <- function(x,
                        filename = "Report",
                        dir = tempdir(),
                        open_file = TRUE,
                        silent = FALSE,...) {
  data <- x
  data$type <- 'Data Report'

  rmd_file <- file.path(system.file(package = "SLAM"), "Report.Rmd")
  rmd <- readLines(rmd_file)

  write(rmd, file = file.path(dir, paste0(filename, ".rmd")))

  if(!silent) message("Rendering markdown file to HTML: ", file.path(dir, paste0(filename, ".html")))

  out <- rmarkdown::render(file.path(dir, paste0(filename, ".rmd")), "html_document", paste0(filename, ".html"), dir,
                           output_options = list(df_print = "paged"), quiet = TRUE)
  if(!silent) message("Rendering complete.")

  if(open_file) browseURL(out)
  return(invisible(out))
}



report_TS <- function(data) {
  df <- data.frame(Year=data$Year,
                   Month=data$Month,
                   Effort=data$Effort,
                   Effort_SD=data$Effort_SD,
                   Index=data$CPUE,
                   Index_SD=data$CPUE_SD
  )
  df$Month_n <- match(df$Month, month.abb)
  df$Date <- as.Date(paste(df$Year, df$Month_n, 01, sep='-'))

  # Effort
  mu <- log(df$Effort) -0.5*df$Effort_SD^2
  Lower <- qlnorm(0.1, mu,df$Effort_SD)
  Upper <- qlnorm(0.9, mu,df$Effort_SD)
  dfE <- df %>% select(Year,Month, Effort, Effort_SD, Date)
  dfE$ymin <- Lower
  dfE$ymax <- Upper

  p1 <- ggplot(dfE, aes(x=Date, y=Effort, ymin=ymin, ymax=ymax)) +
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
                 expand=c(0.01,0)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank())

  # Index
  mu <- log(df$Index) -0.5*df$Index_SD^2
  Lower <- qlnorm(0.1, mu,df$Index_SD)
  Upper <- qlnorm(0.9, mu,df$Index_SD)
  dfI <- df %>% select(Year,Month, Index, Index_SD, Date)
  dfI$ymin <- Lower
  dfI$ymax <- Upper

  p2 <- ggplot(dfI, aes(x=Date, y=Index, ymin=ymin, ymax=ymax)) +
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

  p <- cowplot::plot_grid(p1,p2, nrow=2, labels=c('a)', 'b)'),
                          rel_heights = c(0.85,1))
  list(df=df, p=p)
}



report_CAW <- function(data) {
  df_list <- list()
  for (i in 1:data$n_month) {
    df_list[[i]] <- data.frame(Year=data$Year[i],
                               Month=data$Month[i],
                               Weight=data$Weight_Mids,
                               Count=data$CAW[,i])
  }
  df <- do.call('rbind', df_list)
  df$Month <- factor(df$Month, ordered = TRUE, levels=month.abb)
  p <- ggplot(df, aes(x=Weight, y=Count)) +
    facet_grid(Month~Year) +
    geom_bar(stat='identity')  +
    expand_limits(y=0) +
    theme_bw()
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

