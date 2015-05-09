\name{stackerReg}
\alias{stackerReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stackerReg(models, train_mat, pred_mat, cost, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
%%     ~~Describe \code{models} here~~
}
  \item{train_mat}{
%%     ~~Describe \code{train_mat} here~~
}
  \item{pred_mat}{
%%     ~~Describe \code{pred_mat} here~~
}
  \item{cost}{
%%     ~~Describe \code{cost} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (models, train_mat, pred_mat, cost, type) 
{
    l1m <- LiblineaR(data = train_mat[, models], target = training$Class, 
        type = type, cost = cost, epsilon = 0.01, svr_eps = NULL, 
        bias = TRUE, wi = NULL, cross = 0, verbose = FALSE)
    log_ensemble <- predict(l1m, newx = pred_mat[, models], proba = TRUE)
    return(log_ensemble$probabilities[, 2])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
