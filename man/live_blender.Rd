\name{live_blender}
\alias{live_blender}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
live_blender(dataset, models, Cl.thickness, Cell.size, Cell.shape, Marg.adhesion, Epith.c.size, Bare.nuclei, Bl.cromatin, Normal.nucleoli, Mitoses)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{models}{
%%     ~~Describe \code{models} here~~
}
  \item{Cl.thickness}{
%%     ~~Describe \code{Cl.thickness} here~~
}
  \item{Cell.size}{
%%     ~~Describe \code{Cell.size} here~~
}
  \item{Cell.shape}{
%%     ~~Describe \code{Cell.shape} here~~
}
  \item{Marg.adhesion}{
%%     ~~Describe \code{Marg.adhesion} here~~
}
  \item{Epith.c.size}{
%%     ~~Describe \code{Epith.c.size} here~~
}
  \item{Bare.nuclei}{
%%     ~~Describe \code{Bare.nuclei} here~~
}
  \item{Bl.cromatin}{
%%     ~~Describe \code{Bl.cromatin} here~~
}
  \item{Normal.nucleoli}{
%%     ~~Describe \code{Normal.nucleoli} here~~
}
  \item{Mitoses}{
%%     ~~Describe \code{Mitoses} here~~
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
function (dataset, models, Cl.thickness, Cell.size, Cell.shape, 
    Marg.adhesion, Epith.c.size, Bare.nuclei, Bl.cromatin, Normal.nucleoli, 
    Mitoses) 
{
    df <- data.frame(Cl.thickness = Cl.thickness, Cell.size = Cell.size, 
        Cell.shape = Cell.shape, Marg.adhesion = Marg.adhesion, 
        Epith.c.size = Epith.c.size, Bare.nuclei = Bare.nuclei, 
        Bl.cromatin = Bl.cromatin, Normal.nucleoli = Normal.nucleoli, 
        Mitoses = Mitoses)
    df <- data.frame(lapply(df, as.factor))
    model_list <- modelList(dataset = dataset, models = models)
    pred_vec <- predict(model_list, newdata = df, type = "prob", 
        outcome = 2)
    mean(pred_vec)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
