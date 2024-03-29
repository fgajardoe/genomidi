% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Genomidi.R
\name{playSymphony}
\alias{playSymphony}
\title{playSymphony}
\usage{
playSymphony(symphony, play = FALSE)
}
\arguments{
\item{symphony}{An already composed object of class Symphony (See the composeSymphony function)}

\item{play}{Boolean. Enable/Disable the automatic playing of the symphony (Default: FALSE. Currently not working)}
}
\value{
A text (.asc) and a MIDI (.mid) file.
}
\description{
This function produces the files required for listening the genomic symphony.
}
\examples{
playSymphony(symphony)
}
