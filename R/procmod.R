#' @title Informative Procrustean Matrix Correlation
#' @name procmod
#' @description Estimates corrected Procrustean correlation between
#'              matrices for removing overfitting effect.
#' @details
#' The functions in the ProcMod package aims to estimate and to test correlation
#' between matrices, correcting for the spurious correlations because of the
#' over-fitting effect.
#'
#' The ProcMod package is developed on the metabarcoding.org gitlab
#' (https://git.metabarcoding.org/lecasofts/ProcMod).
#' The gitlab of metabarcoding.org provides up-to-date information and
#' forums for bug reports.
#'
#' @author Christelle Gonindard-Melodelima
#' @author Eric Coissac
#'
#' @docType package
#' @importFrom Rdpack reprompt
NULL

#' DNA metabarcoding Australia South-North Gradient
#'
#' This data set of five \code{data.frame}
#' is a simplified version of a full data set
#' describing biodiversity changes along a South-North
#' gradient on the Australian East Coast, from Sidney to
#' North Cap using a DNA metabarcoding approach.
#' The gradient is constituted of 21 locations.
#'
#' \describe{
#' \item{bacteria}{is a 21 x 2150 \code{data.frame} describing bacterial
#'   community at each one of the 21 locations.
#'   Each number is the relative frequency of a molecular operational
#'   taxonomy unit (MOTU) at a site after data cleaning  and
#'   averaging of 135 pontual measures.}
#'
#' \item{bacteria}{is a 21 x 1393 \code{data.frame} describing eukariote
#'   community at each one of the 21 locations.
#'   Each number is the relative frequency of a molecular operational
#'   taxonomy unit (MOTU) at a site after data cleaning  and
#'   averaging of 135 pontual measures.}
#'
#' \item{climat}{is a 21 x 6 \code{data.frame} describing climatic conditions
#'   at each site using worldclim descriptors (https://www.worldclim.org).
#'   \describe{
#'   \item{Aspect}{}
#'   \item{TempSeasonality}{}
#'   \item{MaxMonTemp}{Max Temperature of Warmest Month}
#'   \item{MeanMonTempRange}{}
#'   \item{AnnMeanTemp}{}
#'   \item{Isothemality}{Mean Diurnal Range / Temperature Annual Range, with
#'     \describe{
#'        \item{Mean Diurnal Range}{Mean of monthly (max temp - min temp)}
#'        \item{Temperature Annual Range}{Max Temperature of Warmest Month - Min Temperature of Coldest Month}
#'      }
#'   }}}
#'
#'  \item{soil}{s a 21 x 6 \code{data.frame} describing soil chemistery
#'   at each site.
#'   Each variable is reduced and centered
#'   \describe{
#'   \item{KLg}{Logarithm of the potassium concentration}
#'   \item{pH}{Soil Ph}
#'   \item{AlLg}{Logarithm of the aluminium concentration}
#'   \item{FeLg}{Logarithm of the iron concentration}
#'   \item{PLg}{Logarithm of the phosphorus concentration}
#'   \item{SLg}{Logarithm of the sulphur concentration}
#'   \item{CaLg}{Logarithm of the calcium concentration}
#'   \item{MgLg}{Logarithm of the magnesium concentration}
#'   \item{MnLg}{Logarithm of the manganese concentration}
#'   \item{CNratio}{carbon / nitrogen concentration ratio}
#'   \item{CLg}{Logarithm of the carbon concentration}
#'   \item{NLg}{Logarithm of the nitrogen concentration}
#'   }}
#'
#'  \item{geography}{}
#' }
#'
#' @docType data
#' @usage data(eukaryotes)
#' @keywords datasets
#' @format five data.frame of 21 rows
#' @rdname australia
#' @author Christelle Gonindard-Melodelima
#' @author Eric Coissac
"eukaryotes"

#' @docType data
#' @usage data(bacteria)
#' @rdname australia
"bacteria"

#' @docType data
#' @usage data(climat)
#' @rdname australia
"climat"

#' @docType data
#' @usage data(soil)
#' @rdname australia
"soil"

#' @docType data
#' @usage data(geography)
#' @rdname australia
"geography"
