#' Austrian Weather Data 2021
#'
#' This data is a subset of the GeoSphere Austria monthly weather data of 2021
#' averaged using the median. Stations with missing values are removed.
#'
#'
#' @format
#' A data frame with 183 rows and 10 columns:
#' \describe{
#'   \item{name}{Unique name of the weather station in German.}
#'   \item{lon, lat}{Longitude and latitude of the weather station.}
#'   \item{alt}{Altitude of the weather station (meter).}
#'   \item{p}{Average air pressure (hPa).}
#'   \item{s}{Monthly sum of sunshine duration (hours).}
#'   \item{vv}{Wind velocity (meter/second).}
#'   \item{t}{Air temperature in 2 meters above the ground in (°C).}
#'   \item{rsum}{Average daily sum of precipitation (mm).}
#'   \item{rel}{Relative air humidity (percent).}
#' }
#' @source
#' The original data was downloaded here (December 2022): \url{https://data.hub.zamg.ac.at/dataset/klima-v1-1m}.
#'
#' @references
#' Data Source: GeoSphere Austria - \url{https://data.hub.zamg.ac.at}.
#'
#' @examples
#' data(weatherAUT2021)
#' summary(weatherAUT2021)
"weatherAUT2021"
