% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotEncoding.R
\name{plot.fmca}
\alias{plot.fmca}
\title{Plot the optimal encoding}
\usage{
\method{plot}{fmca}(
  x,
  harm = 1,
  states = NULL,
  addCI = FALSE,
  coeff = 1.96,
  col = NULL,
  nx = 128,
  ...
)
}
\arguments{
\item{x}{output of \code{\link{compute_optimal_encoding}} function}

\item{harm}{harmonic to use for the encoding}

\item{states}{states to plot (default = NULL, it plots all states)}

\item{addCI}{if TRUE, plot confidence interval (only when \code{computeCI = TRUE} in \link{compute_optimal_encoding})}

\item{coeff}{the confidence interval is computed with +- coeff * the standard deviation}

\item{col}{a vector containing color for each state}

\item{nx}{number of time points used to plot}

\item{...}{not used}
}
\value{
a \code{ggplot} object that can be modified using \code{ggplot2} package.
}
\description{
Plot the optimal encoding
}
\details{
The encoding for the harmonic \code{h} is \eqn{a_{x}^{(h)} \approx \sum_{i=1}^m \alpha_{x,i}^{(h)}\phi_i}.
}
\examples{
# Simulate the Jukes-Cantor model of nucleotide replacement
K <- 4
Tmax <- 6
PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
lambda_PJK <- c(1, 1, 1, 1)
d_JK <- generate_Markov(n = 10, K = K, P = PJK, lambda = lambda_PJK, Tmax = Tmax)
d_JK2 <- cut_data(d_JK, Tmax)

# create basis object
m <- 6
b <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)
\donttest{
# compute encoding
encoding <- compute_optimal_encoding(d_JK2, b, computeCI = FALSE, nCores = 1)

# plot the encoding produced by the first harmonic
plot(encoding)


# modify the plot using ggplot2
library(ggplot2)
plot(encoding, harm = 2, col = c("red", "blue", "darkgreen", "yellow")) +
  labs(title = "Optimal encoding")
}

}
\seealso{
Other encoding functions: 
\code{\link{compute_optimal_encoding}()},
\code{\link{get_encoding}()},
\code{\link{plotComponent}()},
\code{\link{plotEigenvalues}()},
\code{\link{predict.fmca}()},
\code{\link{print.fmca}()},
\code{\link{summary.fmca}()}
}
\author{
Quentin Grimonprez
}
\concept{encoding functions}
