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

\name{wccHeatmap}
\alias{wccHeatmap}

\title{wccHeatmap}

\description{
This function plots a heatmap of a chosen aggregate statistic from the results of \code{wccFindDyadParam}
}

\usage{
wccHeatmap(xparam=NA, yparam=NA, aggstat=NA, xlabel=NA, ylabel=NA, pdffile=NA)
}

\arguments{
\item{xparam}{One of the parameter columns from the dataframe returned by \code{wccFindDyadParam}.  This parameter will be used as the x axis of the heatmap plot.  }
\item{yparam}{One of the parameter columns from the dataframe returned by \code{wccFindDyadParam}.  This parameter will be used as the y axis of the heatmap plot. }
\item{aggstat}{One of the aggregate statistics columns from the dataframe returned by wccFindDyadParam.  \code{aggstat} will be the value that determines the color in each cell of the heatmap.}
\item{xlabel}{A character string that will be used as the label for the x axis.}
\item{ylabel}{A character string that will be used as the label for the y axis.}
\item{pdffile}{A character string with the filename of a pdf file to create.  If pdffile=\code{NA} then the plot is sent to the default device.}
}

\details{
This function produces a heat map plot of selected aggregated statistics returned by \code{wccFindDyadParam} for combinations of selected parameter values.  First run \code{wccFindDyadParam} using two or more parameter vectors for a grid search.  Next select two colums of parameters and one aggregated statistic to plot and send it to wccHeatmap.  If the optional parameter pdffile exists, then a publication quality pdf file will be saved.  Otherwise the heatmap will be displayed on the default device.
}
    
\value{
Returns nothing.
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.}

\seealso{
\code{\link{wccFindDyadParam}} for the definition of the dataframe columns used by \code{wccHeatmap}. 
}

\examples{
\donttest{
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
# Plot the heatmap
wccHeatmap(xparam=wccFDPout$wMax, yparam=wccFDPout$tMax, aggstat=wccFDPout$maxMean,
    xlabel="Window Max", ylabel="Max Lag")
}
}
