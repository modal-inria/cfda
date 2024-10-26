#' Family life states from the Swiss Household Panel biographical survey
#'
#' 2000 16 year-long family life sequences built from the retrospective biographical survey carried out by the Swiss
#' Household Panel (SHP) in 2002. Data from \code{TraMineR} package.
#'
#' @name biofam2
#' @aliases biofam2
#' @docType data
#' @keywords data
#'
#' @usage data(biofam2)
#'
#' @format A data.frame containing three columns:
#' \itemize{
#'   \item \emph{id} id of individuals (2000 different ids)
#'   \item \emph{time} age in years where a change occurs
#'   \item \emph{state} new state.
#' }
#'
#'
#' @details
#' The biofam2 dataset derives from the biofam dataset from \code{TraMineR} package.
#' The biofam2 format is adapted to \code{cfda} functions.
#' The biofam data set was constructed by Müller et al. (2007) from the data of the retrospective biographical survey
#' carried out by the Swiss Household Panel (SHP) in 2002.
#' The data set contains sequences of family life states from age 15 to 30 (sequence length is 16).
#' The sequences are a sample of 2000 sequences of those created from the SHP biographical survey.
#' It includes only individuals who were at least 30 years old at the time of the survey.
#' The biofam data set describes family life courses of 2000 individuals born between 1909 and 1972.
#'
#' The eight states are defined from the combination of five basic states, namely Living with parents (Parent),
#' Left home (Left), Married (Marr), Having Children (Child), Divorced:
#' "Parent", "Left", "Married", "Left+Marr", "Child", "Left+Child", "Left+Marr+Child", "Divorced"
#'
#' @source Swiss Household Panel https://forscenter.ch/projects/swiss-household-panel/
#'
#' @references Müller, N. S., M. Studer, G. Ritschard (2007). Classification de parcours de vie à l'aide de l'optimal matching.
#' In XIVe Rencontre de la Société francophone de classification (SFC 2007), Paris, 5 - 7 septembre 2007, pp. 157–160.
#'
#' @examples
#' data(biofam2)
#' head(biofam2)
#'
#' plotData(biofam2)
#'
#' \donttest{
#' # It is recommended to increase the number of cores to reduce computation time
#' set.seed(42)
#' basis <- create.bspline.basis(c(15, 30), nbasis = 4, norder = 4)
#' fmca <- compute_optimal_encoding(biofam2, basis, nCores = 2)
#'
#' plot(fmca, harm = 1)
#' plot(fmca, harm = 2)
#' plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE)
#' plotComponent(fmca, comp = c(1, 2), addNames = FALSE)
#' }
#' @family datasets
NULL


#' Care trajectories
#'
#' Care trajectories of patients diagnosed with a serious and chronic condition
#'
#' @name care
#' @aliases care
#' @docType data
#' @keywords data
#'
#' @usage data(care)
#'
#' @format A data.frame containing three columns:
#' \itemize{
#'   \item \emph{id} id of individuals (2929 different ids)
#'   \item \emph{time} number of months since the diagnosis
#'   \item \emph{state} new state.
#' }
#'
#'
#' @details
#' In this study, patients were followed from the time they were diagnosed with a serious and chronic condition
#' and their care trajectories were tracked monthly from the time of diagnosis.
#' The status variable contains the care status of each individual for each month of follow-up.
#' Trajectories have different lengths.
#'
#' The four states are:
#' \itemize{
#'  \item{D: diagnosed, but not in care}
#'  \item{C: in care, but not on treatment}
#'  \item{T: on treatment, but infection not suppressed}
#'  \item{S: on treatment and suppressed infection}
#' }
#'
#' @source https://larmarange.github.io/analyse-R/data/care_trajectories.RData
#' https://larmarange.github.io/analyse-R/trajectoires-de-soins.html
#'
#'
#' @examples
#' data(care)
#' head(care)
#'
#' plotData(care)
#'
#' # Individuals has not the same length. In order to compute the encoding,
#' # we keep individuals with at least 18 months of history and work
#' # with the 18 first months.
#' duration <- compute_duration(care)
#' idToKeep <- as.numeric(names(duration[duration >= 18]))
#' care2 <- cut_data(care[care$id %in% idToKeep, ], 18)
#' head(care2)
#' \donttest{
#' # It is recommended to increase the number of cores to reduce computation time
#' set.seed(42)
#' basis <- create.bspline.basis(c(0, 18), nbasis = 10, norder = 4)
#' fmca <- compute_optimal_encoding(care2, basis, nCores = 2)
#'
#' plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE)
#' plot(fmca)
#' plot(fmca, addCI = TRUE)
#' plotComponent(fmca, addNames = FALSE)
#' }
#' @family datasets
NULL


#' Flours dataset
#'
#' Resistance of dough during the kneading process
#'
#' @name flours
#' @aliases flours
#' @docType data
#' @keywords data
#'
#' @usage data(flours)
#'
#' @format \code{flours} is a list of 3 elements:
#' \itemize{
#'   \item \code{data} A matrix of size 241*115 containing the resistance of dough (measured every 2s) during the kneading
#' process. One dough batch = 1 column
#'   \item \code{quality} Quality of cookies baked with the associated dough (1=Good, 2=Medium, 3=Bad)
#'   \item \code{time} time values
#' }
#'
#' @examples
#' data(flours)
#'
#' matplot(flours$time, flours$data, col = flours$quality, type = "l", lty = 1)
#'
#' # convert to categorical data
#' flours_cfd <- convertToCfd(flours$data,
#'     breaks = c(-Inf, 150, 300, 450, 600, Inf),
#'     times = flours$time
#' )
#'
#' plotData(flours_cfd, group = flours$quality)
#'
#'
#' # convert to categorical data after projecting in a basis of functions
#' basis <- create.bspline.basis(c(0, 480), nbasis = 10)
#' flours_fd <- Data2fd(flours$time, flours$data, basis)
#' plot(flours_fd)
#'
#' flours_cfd2 <- convertToCfd(flours_fd, breaks = c(-Inf, 150, 300, 450, 600, Inf))
#'
#' plotData(flours_cfd2, group = flours$quality)
#' @family datasets
#'
NULL
