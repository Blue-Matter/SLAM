% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Report.R
\name{plot_selectivity}
\alias{plot_selectivity}
\title{Plot the predicted selectivity-at-age schedule}
\usage{
plot_selectivity(assess)
}
\arguments{
\item{assess}{An object of class \code{Assess}}
}
\value{
A list with a data.frame and a ggplot object
}
\description{
Plot the predicted selectivity-at-age schedule
}
\examples{
MyData <- Import(Example_Data()[2])
MyAssess <- Assess(MyData)
plot_selectivity(MyAssess)
}
