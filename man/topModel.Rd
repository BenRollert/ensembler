\name{topModel}
\alias{topModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
topModel(n, ordered_acc, pred_mat, reference)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{ordered_acc}{
%%     ~~Describe \code{ordered_acc} here~~
}
  \item{pred_mat}{
%%     ~~Describe \code{pred_mat} here~~
}
  \item{reference}{
%%     ~~Describe \code{reference} here~~
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
function (n, ordered_acc, pred_mat, reference) 
{
    if (n > 1) {
        top_mean <- apply(pred_mat[, ordered_acc$model[1:n]], 
            1, mean)
        pred <- prediction(top_mean, reference)
        perf <- performance(pred, "tpr", "fpr")
        auc.tmp <- performance(pred, "auc")
        auc <- as.numeric(auc.tmp@y.values)
        cbind(as.integer(n), auc, paste(ordered_acc$model[1:n], 
            collapse = ", "))
    }
    else {
        top_mean <- pred_mat[, ordered_acc$model[1]]
        pred <- prediction(top_mean, reference)
        perf <- performance(pred, "tpr", "fpr")
        auc.tmp <- performance(pred, "auc")
        auc <- as.numeric(auc.tmp@y.values)
        cbind(as.integer(n), auc, paste(ordered_acc$model[1]))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
