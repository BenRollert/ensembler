\name{topComb}
\alias{topComb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
topComb(ordered_acc, pred_mat, reference)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (ordered_acc, pred_mat, reference) 
{
    top_ens <- sapply(1:length(models), topModel, ordered_acc = ordered_acc, 
        pred_mat = pred_mat, reference = reference)
    top_m <- data.table(t(top_ens))
    setnames(top_m, c("number_models", "AUC", "model_names"))
    top_m
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
