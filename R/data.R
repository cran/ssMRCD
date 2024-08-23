#' Austrian Weather Data 2021
#'
#' This data is a subset of the GeoSphere Austria monthly weather data of 2021
#' averaged using the median. Stations with missing values are removed.
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
#' The original data was downloaded here (December 2022): \url{https://data.hub.geosphere.at/dataset/klima-v1-1m}.
#'
#' @references
#' Data Source: GeoSphere Austria - \url{https://data.hub.geosphere.at}.
#'
#' @examples
#' data(weatherAUT2021)
#' summary(weatherAUT2021)
"weatherAUT2021"




#' Vienna Weather Time Series (1960-2023)
#'
#' This data is a subset of the GeoSphere Austria daily weather data of the time 1960-2023
#' for the weather station Hohe Warte in Vienna.
#'
#' @format
#' A data frame with 23372 rows and 18 columns including 13 weather measurements:
#' \describe{
#'   \item{time}{Time of measurement in date format.}
#'   \item{cloud_cover}{Daily mean of cloud coverage, values between 1 and 100.}
#'   \item{global_radiation}{Daily sum of global radiation (J/cm^2).}
#'   \item{vapor_pressure}{Daily mean of vapour pressuer (hPa).}
#'   \item{max_wind_speed}{Maximal wind speed (m/s).}
#'   \item{air_pressure}{Daily mean of air pressure (hPa).}
#'   \item{relative_humidity}{Daily mean of relative humidity (percent).}
#'   \item{precipitation}{Daily sum of precepitation (mm).}
#'   \item{sight}{Sight distance at 1pm (m).}
#'   \item{sunshine_duration}{Daily sum of sunshine duration (h).}
#'   \item{temperature_max}{Daily maximum of temperature at 2m air height (°C).}
#'   \item{temperature_min}{Daily minimum of temperature at 2m air height (°C).}
#'   \item{temperature_mean}{Daily mean of temperature at 2m air height (°C).}
#'   \item{wind_velocity}{Daily mean of wind speed (m/s).}
#'   \item{year}{Year of measurement.}
#'   \item{month}{Month of measurement.}
#'   \item{day}{Day of the year of measurement.}
#'   \item{season}{Season of measuremen (1 = winter, 2 = spring, 3 = summer, 4 = fall).}
#' }
#'
#
#' @source
#' The original data was downloaded here (April 2024): \url{https://data.hub.geosphere.at/dataset/klima-v2-1d}.
#'
#' @references
#' Data Source: GeoSphere Austria - \url{https://data.hub.geosphere.at}.
#'
#' @examples
#' data(weatherHoheWarte)
#' summary(weatherHoheWarte)
"weatherHoheWarte"
