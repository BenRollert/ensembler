\name{stackerL1}
\alias{stackerL1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stackerL1(dataset, models, traindf, testdf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{models}{
%%     ~~Describe \code{models} here~~
}
  \item{traindf}{
%%     ~~Describe \code{traindf} here~~
}
  \item{testdf}{
%%     ~~Describe \code{testdf} here~~
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
function (dataset, models, traindf, testdf) 
{
    train_mat <- ensembler(dataset = dataset, models = models, 
        newdata = traindf, outcome = 2, type = "prob")
    l1m <- LiblineaR(data = train_mat, target = training$Class, 
        type = 6, cost = 1, epsilon = 0.01, svr_eps = NULL, bias = TRUE, 
        wi = NULL, cross = 0, verbose = FALSE)
    pred_mat <- ensembler(dataset = dataset, models = models, 
        newdata = testdf, outcome = 2, type = "prob")
    log_ensemble <- predict(l1m, newx = pred_mat, proba = TRUE)
    return(log_ensemble$probabilities[, 2])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
