\name{cof_selection}

\alias{cof_selection}

\alias{plot.fullsib_cofactors }

\title{Cofactors Selection procedure for CIM analysis}

\description{
Performs the selection of markers that will be used as cofactors in CIM
analysis. The procedure is fully based on \code{step} function from
\pkg{stats} package
}

\usage{
cof_selection(fullsib, pheno.col=1, addcovar = NULL, k=2, n.cofactor=10,
              stoppage.df=NULL, trace=0, thres.effect=1, selection=1)

\method{plot}{fullsib_cofactors}(x, horizontal = FALSE, grid = FALSE,
cex.mrk = 1, cex.grp = 0.75,\ldots)
}

\arguments{
  \item{fullsib}{an object from class \emph{fullsib}.}

  \item{pheno.col}{Column number in the phenotype matrix (present in
    \emph{fullsib} object) which should be used as the phenotype.}

  \item{addcovar}{Additive covariates. If it is used, one must indicate the
    design matrix for those source of variation. It should be noted that
    additive covariates is included in the model as fixed effects, under
    ordinary linear regression.}
  
  \item{k}{value to be used as penalty in the model selection. If
    \code{k=2} is used corresponds to the Akaike Information Criterion;
    for \code{k=log(n)} is the Bayesian information criterion (BIC) 
    or sometimes called Schwarz Information criterion.}
  
  \item{n.cofactor}{maximum number of cofactors that will be present on
    the model. Default is to consider a model with maximum of 10
    cofactors. See details.}
  
  \item{stoppage.df}{controls the maximum number of markers to be
    included in the model, using the model's degree of freedom as a
    limit to stop the process. Default is NULL (not considering this
    argument),  but \eqn{2\times\sqrt(n)} (being \eqn{n} is number of
    individuals phenotyped) can be used as a limit to control the
    maximum number of cofactors.} 
  
  \item{trace}{If zero is considered none information about the
    selection process is showed, if non-negative integer is considered,
    one can see additional information during the running process. This
    argument is the same from \code{step}, once the model selection is
    based on it.}
  
  \item{thres.effect}{Threshold value to remove non significative effect
    for selected cofactors. Default is \code{thres.effect = 1}, none
    cofactors effect is removed from the analysis.}
  
  \item{selection}{integer indicating which markers should be used for
    this step. It may assumes values of \sQuote{0}, \sQuote{1} and
    \sQuote{2}. If \sQuote{0}: only markers that are not included in the 
    map are used on the selection process. If \sQuote{1}: the
    selection process just consider markers that are on the linkage
    groups. And, if \sQuote{2}: all markers are used for selection
    model.}

  \item{x}{An object from class \emph{fullsib_cofactors}.}
  
  \item{horizontal}{if \code{TRUE}, indicates that the map should be plotted
    horizontally. Default is \code{FALSE}.}
  
  \item{grid}{if \code{TRUE}, displays a grid in the background. Default is
    \code{FALSE}.}
  
  \item{cex.mrk}{the magnification to be used for markers.}
  
  \item{cex.grp}{the magnification to be used for group axis
    annotation.}

  \item{\ldots}{Further arguments, passed to other methods. Currently
  ignored.}


}

\details{
  This function prepares the information contained in the object of
  class \emph{fullsib} for using on a multiple linear regression. By
  default the independent variable are molecular markers that are
  present on the linkages groups (LG), optionally one can also consider
  the markers that were not placed on the map. Genetic Predictors are
  obtained by the multipoint probability at the exactly point that the
  markers is placed on the LG and multiplied by the estimable
  contrasts.

  The genetic predictors are regressed with the phenotype using the
  \code{step} function. This strategy was done to get advantages on
  model selection tools avaliable on \R. In this case, the function
  performs stepwise procedure using information criteria, such 
  as AIC and BIC (defined by the user). If one would like to use P-value
  (or another criteria) to select cofactors, one can do by oneself and
  add these cofactors using \code{cof_definition} function.

  As the selection process is based on \code{step}, one can follow the
  selection just using \code{trace} argument (same present in
  \code{step}) higher than zero.

  Remembering that this is an intermediate step for CIM analysis, and
  sometimes many markers can be selected as cofactors reducing
  statistical power. To solve this situation, it is considered some
  arguments that can control the number of included markers, for example
  \sQuote{n.cofactor} controls the maximum number of included markers
  established by the user. Other option is \sQuote{stoppage.df}, which
  the user can define the maximum number of degree of freedom that model
  can have. In both case, the selection criteria is used to include
  cofactors, but when model reaches one of these limits FIRST the
  process is stopped.

  An other way to control the selection of markers is to verify the non
  significative effects. For this, a linear model is fitted \code{lm}
  function and \code{summary.lm} is used to obtain the p-values of the
  cofactors effects. So \sQuote{thres.effect} is a threshold probability
  used to determine which effect are significative. If 1 is used none
  effect is considered non significative and removed from the analysis.

  The method plot was designed for the user evaluated the cofactors
  saturation and dispersion along the genome. Helping with the decision
  of window size dimension and best selection options.
}

\value{
  An object of class \emph{fullsib_cofactors} returned, which has the
  same structure of an object of class \emph{fullsib} with the inclusion
  of an extra component (\sQuote{cofactor}) that is a list with the
  components: \code{names.cof}, \code{matrix.cof} and \code{trait.cof}.

  \code{names.cof} is a data frame showing which are the selected
  markers and their linkage groups. \code{NA} value is used for markers
  that are not placed on linkage groups.

  \code{matrix.cof} is a matrix with contains all the cofactors
  effect. During the CIM analysis columns are dropped in function of
  window size.

  \code{trait.cof} is the indication for which pheno.col was considered
  for cofactor selection step.
}

\references{
  Burnham, K.P., Anderson, D.R. (2004) Understanding AIC and BIC in
  Model Selection. \emph{Sociological Methods & Research}, 33: 261-304

  Sakamoto, Y., Ishiguro, M., and Kitagawa, G. (1986), Akaike
  \emph{Information Criterion Statistics}, Tokyo: KTK Scientific Publishers.
}  
  
  

\author{Rodrigo Gazaffi, \email{rgazaffi@gmail.com}
  \code{plot_fullsib_cofactors } was based on \code{draw.map} that was
  written by Marcelo Mollinari \email{mmollina@esalq.usp.br} and modified
  by Rodrigo Gazaffi, allowing cofactors inclusion.
}

\seealso{
  \code{\link[fullsibQTL]{create_fullsib}},

  \code{\link[fullsibQTL]{cof_definition}},
  
  \code{step} from \pkg{stats} package
  
 
}
\examples{

  data(example_QTLfullsib)

  fullsib <- create_fullsib(example_QTLfullsib,
                            list(LG1_final, LG2_final, LG3_final, LG4_final),
                            step=0,map.function="kosambi",condIndex=3.5)



  ###############################################
  ## cofactor selection using BIC (n.ind = 300)

  ### just using markers that are placed on linkage groups (default)
  cofs.fs <- cof_selection(fullsib, pheno.col=1, k = log(300), selection =1)
  cofs.fs
  plot(cofs.fs)


  ###just markers that are not placed in LG
  ###Attention: just used if there are markers that are not placed on
  ###genetic map
  \dontrun{
  cofs.fs1 <- cof_selection(fullsib, pheno.col=1, k = log(300), selection =0)
  }
  ###all markers are used
  ###Attention: just used if fullsib$unlinked.mkr is not NULL
  \dontrun{
   cofs.fs2 <- cof_selection(fullsib, pheno.col=1, k = log(300), selection=2)
   }



  ###############################################
  ## cofactor selection using AIC criteria

  \dontrun{  
  ### using AIC alowing at maximum of 20 markers
  cofs.fs3 <- cof_selection(fullsib, pheno.col=1, k = 2, n.cofactor =  20)
  plot(cofs.fs3)

  ### using AIC alowing at maximum of 5 markers
  cofs.fs4 <- cof_selection(fullsib, pheno.col=1, k = 2, n.cofactor = 5)
  plot(cofs.fs4)

  ### using AIC alowing the model with maximum of 18 D.f.
  ### here 6 marker, because each cofactor has 3 effects
  cofs.fs5 <- cof_selection(fullsib, pheno.col=1, k = 2, stoppage.df = 18)
  plot(cofs.fs5)

  ### it is possible to see the selection process
  cofs.fs6 <- cof_selection(fullsib, pheno.col=1, k=2, n.cofactor=5, trace=1)
  plot(cofs.fs6)

  ### Same selection, but removing effects that are not signif. at 5\%
  cofs.fs7 <- cof_selection(fullsib, pheno.col=1, k = 2, n.cofactor = 5, thres.effect = 0.05) 
  }

  ##cofactor selection using additive covariate
  \dontrun{
  covar <- matrix(rep(c(1,-1), each=150), ncol=1)
  cofs.fs8 <- cof_selection(fullsib, pheno.col=2, addcovar=covar, k = 2)
  }

} 


\keyword{utilities}
%\keyword{plot}
