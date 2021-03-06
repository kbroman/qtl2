% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit1.R
\name{fit1}
\alias{fit1}
\title{Fit single-QTL model at a single position}
\usage{
fit1(
  genoprobs,
  pheno,
  kinship = NULL,
  addcovar = NULL,
  nullcovar = NULL,
  intcovar = NULL,
  weights = NULL,
  contrasts = NULL,
  model = c("normal", "binary"),
  zerosum = TRUE,
  se = TRUE,
  hsq = NULL,
  reml = TRUE,
  blup = FALSE,
  ...
)
}
\arguments{
\item{genoprobs}{A matrix of genotype probabilities, individuals x genotypes.
If NULL, we create a single intercept column, matching the individual IDs in \code{pheno}.}

\item{pheno}{A numeric vector of phenotype values (just one phenotype, not a matrix of them)}

\item{kinship}{Optional kinship matrix.}

\item{addcovar}{An optional numeric matrix of additive covariates.}

\item{nullcovar}{An optional numeric matrix of additional additive
covariates that are used under the null hypothesis (of no QTL)
but not under the alternative (with a QTL). This is needed for
the X chromosome, where we might need sex as a additive
covariate under the null hypothesis, but we wouldn't want to
include it under the alternative as it would be collinear with
the QTL effects.}

\item{intcovar}{An optional numeric matrix of interactive covariates.}

\item{weights}{An optional numeric vector of positive weights for the
individuals. As with the other inputs, it must have \code{names}
for individual identifiers.}

\item{contrasts}{An optional numeric matrix of genotype contrasts, size
genotypes x genotypes. For an intercross, you might use
\code{cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0))} to get
mean, additive effect, and dominance effect. The default is the
identity matrix.}

\item{model}{Indicates whether to use a normal model (least
squares) or binary model (logistic regression) for the phenotype.
If \code{model="binary"}, the phenotypes must have values in
\eqn{[0, 1]}.}

\item{zerosum}{If TRUE, force the genotype or allele coefficients
sum to 0 by subtracting their mean and add another column with
the mean. Ignored if \code{contrasts} is provided.}

\item{se}{If TRUE, calculate the standard errors.}

\item{hsq}{(Optional) residual heritability; used only if
\code{kinship} provided.}

\item{reml}{If \code{kinship} provided: if \code{reml=TRUE}, use
REML; otherwise maximum likelihood.}

\item{blup}{If TRUE, fit a model with QTL effects being random, as in \code{\link[=scan1blup]{scan1blup()}}.}

\item{...}{Additional control parameters; see Details;}
}
\value{
A list containing
\itemize{
\item \code{coef} - Vector of estimated coefficients.
\item \code{SE} - Vector of estimated standard errors (included if \code{se=TRUE}).
\item \code{lod} - The overall lod score.
\item \code{ind_lod} - Vector of individual contributions to the LOD score (not provided if \code{kinship} is used).
\item \code{fitted}  - Fitted values.
\item \code{resid}  - Residuals.
If \code{blup==TRUE}, only \code{coef} and \code{SE} are included at present.
}
}
\description{
Fit a single-QTL model at a single putative QTL position and get detailed results
about estimated coefficients and individuals contributions to the LOD score.
}
\details{
For each of the inputs, the row names are used as
individual identifiers, to align individuals.

If \code{kinship} is absent, Haley-Knott regression is performed.
If \code{kinship} is provided, a linear mixed model is used, with a
polygenic effect estimated under the null hypothesis of no (major)
QTL, and then taken as fixed as known in the genome scan.

If \code{contrasts} is provided, the genotype probability matrix,
\eqn{P}, is post-multiplied by the contrasts matrix, \eqn{A}, prior
to fitting the model. So we use \eqn{P \cdot A}{P A} as the \eqn{X}
matrix in the model. One might view the rows of
\ifelse{html}{\out{<em>A</em><sup>-1</sup>}}{\eqn{A^{-1}}}
as the set of contrasts, as the estimated effects are the estimated
genotype effects pre-multiplied by
\ifelse{html}{\out{<em>A</em><sup>-1</sup>}}{\eqn{A^{-1}}}.

The \code{...} argument can contain several additional control
parameters; suspended for simplicity (or confusion, depending on
your point of view). \code{tol} is used as a tolerance value for linear
regression by QR decomposition (in determining whether columns are
linearly dependent on others and should be omitted); default
\code{1e-12}. \code{maxit} is the maximum number of iterations for
converence of the iterative algorithm used when \code{model=binary}.
\code{bintol} is used as a tolerance for converence for the iterative
algorithm used when \code{model=binary}. \code{eta_max} is the maximum value
for the "linear predictor" in the case \code{model="binary"} (a bit of a
technicality to avoid fitted values exactly at 0 or 1).
}
\examples{
# read data
iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
\dontshow{iron <- iron[,"7"] # subset to chr 7}

# insert pseudomarkers into map
map <- insert_pseudomarkers(iron$gmap, step=5)

# calculate genotype probabilities
probs <- calc_genoprob(iron, map, error_prob=0.002)

# grab phenotypes and covariates; ensure that covariates have names attribute
pheno <- iron$pheno[,1]
covar <- match(iron$covar$sex, c("f", "m")) # make numeric
names(covar) <- rownames(iron$covar)

# scan chromosome 7 to find peak
out <- scan1(probs[,"7"], pheno, addcovar=covar)

# find peak position
max_pos <- max(out, map)

# genoprobs at max position
pr_max <- pull_genoprobpos(probs, map, max_pos$chr, max_pos$pos)

# fit QTL model just at that position
out_fit1 <- fit1(pr_max, pheno, addcovar=covar)

}
\references{
Haley CS, Knott SA (1992) A simple
regression method for mapping quantitative trait loci in line
crosses using flanking markers.  Heredity 69:315--324.

Kang HM, Zaitlen NA, Wade CM, Kirby A, Heckerman D, Daly MJ, Eskin
E (2008) Efficient control of population structure in model
organism association mapping. Genetics 178:1709--1723.
}
\seealso{
\code{\link[=pull_genoprobpos]{pull_genoprobpos()}}, \code{\link[=find_marker]{find_marker()}}
}
