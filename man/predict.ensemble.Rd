\name{predict.ensemble}
\alias{predict.ensemble}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predict.ensemble(model_list, newdata, type, outcome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model_list}{
%%     ~~Describe \code{model_list} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{outcome}{
%%     ~~Describe \code{outcome} here~~
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
function (model_list, newdata, type, outcome) 
{
    pred_fact <- function(newdata, type, outcome) {
        function(model) {
            out <- predict(model, newdata = newdata, type = type)
            out <- out[, outcome]
            out
        }
    }
    predict_filt <- pred_fact(newdata = newdata, type = "prob", 
        outcome = 2)
    out <- sapply(model_list, predict_filt)
    if (class(out) == "numeric") {
        names(out) <- sapply(model_list, function(x) x$method)
    }
    else {
        colnames(out) <- sapply(model_list, function(x) x$method)
    }
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
