\name{wcc-package}
\alias{wcc-package}
\alias{wcc-package}
\docType{package}
\title{
wcc
}
\description{
Functions for Windowed Cross Correlation.
}
\details{
Calculates Windowed Cross Correlation for estimating association between pairs of nonstationary time series.
    Provides support for surrogate analysis for nonparametric test of significance.
    Calculates aggregate statistics over a range of parameter values.
    Plots the results as wcc plots and heat maps.
}
\author{
Steven Boker,
Minquan Xu,
Sareena Chadha,
Christopher Welker,
Jingyun Wu,
Pascal Deboeck

Maintainer: Steven Boker <smb3u@virginia.edu>
}
\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.

Moulder, R., Boker, S., Ramseyer, F., & Tschacher, W. (2018). Determining synchrony between behavioral time series: An application of surrogate data generation for establishing falsifiable null-hypotheses. Psychological Methods. 23:4 pp 757–773
}
\keyword{package}

%% Uncomment below to imitate parts of library(help = wcc)
%\section{The \file{DESCRIPTION} File}{\packageDESCRIPTION{wcc}}
%\section{Documentation Index}{\packageIndices{wcc}}

\seealso{
See \link{wccCalc} for the core function that calculates a windowed cross correlation matrix.\cr
See \link{wccPeakPick} for the core function that estimates time lag of maximum association from a windowed cross correlation matrix.\cr
See \link{wccAggregate} to call \code{wccCalc} and \code{wccPeakPick} and then calculate aggregate statistics from a windowed cross correlation matrix.  \cr
See \link{wccPlot} to plot a section of a windowed cross correlation matrix and peak picking line from a file created by \link{wccAggregate}.\cr
See \link{wccSurrogateDyads} to generate a distribution of the aggregate statistics conforming to the null hypothesis that the pairing of the dyad does not matter.\cr
See \link{wccFindDyadParam} to explore a range of parameter values for \code{wccCalc} that are compared between surrogate and real dyads.\cr
See \link{wccHeatmap} to visualize results of selected aggregated statistics for combinations of parameters calculated by \link{wccFindDyadParam}.\cr
}
\examples{
\donttest{
#Create a windowed cross correlation plot
tSeries1 <- sin(c(1:1000)/10) + rnorm(1000, mean=0, sd=.5)
tSeries2 <- sin(1+c(1:1000)/10) + rnorm(1000, mean=0, sd=.5)
wccPlot(inSeries1=tSeries1, inSeries2=tSeries2, startwindow=1, endwindow=500, 
    wMax=100, tMax=100, wInc=1, tInc=1, Lsize=8, pspan=.25, type="Max", 
    samplespersecond=10)

# Create two arrays of timeseries with dyadic dependence
array1 <- matrix(NA, nrow=10, ncol=500)
array2 <- matrix(NA, nrow=10, ncol=500)
for(i in 1:10) {
    array1[i,] <- sin(c(1:500)/runif(1, min=5, max=20)) + rnorm(500, mean=0, sd=.5)
    array2[i,] <- array1[i,] + rnorm(500, mean=0, sd=.5)
}
# Select parameters to explore
wMaxVec <- c(50,100)
tMaxVec <- c(25,50)
wccFDPout <- wccFindDyadParam(inArray1=array1, inArray2=array2, wMaxvector=wMaxVec,
    tMaxvector=tMaxVec, nSurrogates=30, samplespersecond=1, method="r")
# Plot a heatmap
wccHeatmap(xparam=wccFDPout$wMax, yparam=wccFDPout$tMax, aggstat=wccFDPout$maxMean,
    xlabel="Window Max", ylabel="Max Lag")
}
}
