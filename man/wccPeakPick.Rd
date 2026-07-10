%
%   Copyright 2001-2026 by the individuals mentioned in the source code history
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{wccPeakPick}
\alias{wccPeakPick}

\title{wccPeakPick}

\description{
This function finds the maximum or minimum correlation peak nearest to zero lag for each row in a windowed crosscorrelation matrix.
}

\usage{
wccPeakPick(tAllCor=NA, Lsize=8, pspan=.25, type="Max")
}

\arguments{
\item{tAllCor}{A numeric matrix returned by the \code{wccCalc} function.}
\item{Lsize}{A numeric value for the required width of the peak (default=8).  For instance, \code{Lsize} = 8 means that a peak must have 4 values descending on either side in order to be called a peak.  Thus for elapsed time T and lag offset L if a peak is found at \code{tAllCor[T,L]} then \code{tAllCor[T,(L-1):(L-Lsize)] < tAllCor[T,L]} and \code{tAllCor[T(L+1):(L+Lsize)] < tAllCor[T,L]}. An \code{Lsize} that is too small may find spurious peaks that are just noise, whereas an \code{Lsize} that is too large may miss actual peaks.}
\item{pspan}{A numeric value used by the \link{loess} function to indicate what proportion of data to use in order to do the \link{loess} smooth  (default=.25).}
\item{type}{A character string (default=\sQuote{Max}).  \sQuote{Max} indicates that the maximum correlation closest to zero lag should be picked whereas \sQuote{Min} indicates that the minimum correlation closest to zero lag should be picked.}
}

\details{
This function operates on the result of \code{wccCalc}, a windowed crosscorrelation matrix in order to find the peak correlation closest to a lag of zero.  The lag and the value of the peak correlation are returned for each row of the windowed crosscorrelation matrix, i.e., for the elapsed time of each window calculated by \code{wccCalc}.  The \code{wccCalc} windows are sampled at time lagged intervals from one another so that the association between the time series can be estimated even though the time lag of maximum association between the two time series may be itself time varying.  For instance, when two individuals converse, their body motions may be associated such that the leader and follower take turns as the two individuals swap roles as speaker and listener.  

This same kind of nonstationary association has been observed between many physiological time series.  At any particular elapsed time, the maximum correlation closest to a lag of zero (whether a positive or negative lag) is taken to be the lagged offset between the two time series.  For series that have oscillatory structure (such as head nodding in conversation)  the maximum correlation closest to a lag of zero can be interpreted as the phase lag between the series.  The variability in the lag is one measure of the degree of nonstationarity present in the association between two time series.
}
    
\value{
Returns a list with two vectors: \code{maxIndex} is the lag of the peak for each window and \code{maxValue} is the calculated correlation at that peak.  If \code{wccPeakPick} is called with type = \sQuote{Min} then the returned vectors are named \code{minIndex} and \code{minValue}.  If no peak is found within the row tAllCor[T,] then maxIndex[T] and maxValue[T] are returned as NAs
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.}

\seealso{
\code{\link{wccCalc}} for the definition of the windowed crosscorrelation matrix. \code{\link{loess}} for the definition of the smoothing function.
}

\examples{
#Calculate windowed cross correlation between two time series
tSeries1 <- sin(c(1:1000)/10) + rnorm(1000, mean=0, sd=.5)
tSeries2 <- sin(c(1:1000)/15) + rnorm(1000, mean=0, sd=.5)
wccCalcOut <- wccCalc(inSeries1=tSeries1, inSeries2=tSeries2, wMax=50, tMax=50, wInc=1,
     tInc=1, method="r")

#Calculate the peak picking vectors
wccPeakPick(wccCalcOut, Lsize=8, pspan=.25, type="Max")


}
