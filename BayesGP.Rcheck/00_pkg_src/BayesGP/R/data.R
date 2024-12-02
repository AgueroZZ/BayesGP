#' The COVID-19 daily death data in Canada.
#'
#'  A subset of the the COVID-19 daily death data collected between 2020 to 2022. The data is obtained from
#' COVID-19 Data Repository by the Center for Systems Science and Engineering (CSSE) at
#' Johns Hopkins University (Dong et al., 2020).
#'
#' @format 'covid_canada'
#' A data frame with 787 rows and 5 columns:
#' \describe{
#'   \item{Date}{The date of the measurement.}
#'   \item{new_deaths}{The number of new deaths at that date.}
#'   \item{t}{The converted numerical value of 'Date'.}
#'   \item{weekdays 1-6}{Coded as 1 if the date is the corresponding weekday, -1 else if the date is on Sunday, and 0 otherwise.}
#'   \item{index}{The index of that observation.}
#' }
#' @source Dong, E., H. Du, and L. Gardner (2020). An interactive web-based dashboard to track
#' covid-19 in real time. \emph{The Lancet infectious diseases 20 (5)}, 533–534.
"covid_canada"


#' A simulated dataset from the case-crossover model.
#'
#'
#' @format 'ccData'
#' A data frame with 3596 rows and 6 columns.
"ccData"


#' The monthly all-cause mortality for male with age less than 40 in Pennsylvania.
#'
#' @format 'PEN_death'
#' A data frame with 299 rows of observations.
"PEN_death"
