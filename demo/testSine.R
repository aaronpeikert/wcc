#
#   Copyright 2001-2026 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# ---------------------------------------------------------------------
# Program: testSine.R
#  Author: Steven Boker
#    Date: Mon Mar 30 13:29:49 EDT 2026
#
# This program tests wccAggregate.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Mon Mar 30 13:29:54 EDT 2026
#      Created testSine.R
#
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Variables 
# ---------------------------------------------------------------------
#
# ---------------------------------------------------------------------
startTime <- Sys.time()

# ----------------------------------
# Read libraries and set options.

set.seed(260224)

library(wcc)

totalDyads <- 20
seriesLength <- 1000

totalSurrogates <- 50
samplespersecond <- 10

WCCwIncvector <- c(1)
WCCwMaxvector <- c(50,100)
WCCtMaxvector <- c(50,100)
WCCtIncvector <- c(1)
PPLsizevector <- c(8)
PPpspanvector <- c(.25)
testPPtype <- "Max"



# ----------------------------------
# Simulate sine time series.

series1matrix <- matrix(NA, nrow=totalDyads, ncol=seriesLength)
series2matrix <- matrix(NA, nrow=totalDyads, ncol=seriesLength)

for(i in 1:totalDyads) {
    t1matrix <- sin(c(1:seriesLength)/runif(1, min=5, max=20)) + rnorm(seriesLength, mean=0, sd=.5)
    t2matrix <- sin(c(1:seriesLength)/runif(1, min=5, max=20)) + rnorm(seriesLength, mean=0, sd=.5)
    lagoffset <- floor(runif(1, min=5, max=20))
    predictor <- c(1:seriesLength)/runif(1, min=20, max=30)
    series1matrix[i,1:(lagoffset-1)] <- rnorm((lagoffset-1), mean=0, sd=.5)
    series2matrix[i,1:(lagoffset-1)] <- rnorm((lagoffset-1), mean=0, sd=.5)
    tsel1 <- 1:(seriesLength-(lagoffset-1))
    tsel2 <- lagoffset:seriesLength
    series1matrix[i,tsel2] <- t1matrix[tsel2] * sin(predictor[tsel1]) + t2matrix[tsel1] * cos(predictor[tsel1]) + rnorm(seriesLength-(lagoffset-1), mean=0, sd=.1)
    series2matrix[i,tsel2] <- t2matrix[tsel2] * cos(predictor[tsel1]) + t1matrix[tsel1] * sin(predictor[tsel1]) + rnorm(seriesLength-(lagoffset-1), mean=0, sd=.1)
}

# ----------------------------------
# Test wccAggregate and plot with several timeseries.

WCCwInc <- 1
WCCwMax <- 100
WCCtMax <- 100
WCCtInc <- 1
PPLsize <- 8
PPpspan <- .4
PPtype <- "Max"

for (i in 1:10) {
    startwindow <- 50
    endwindow <- 550
    # uncomment then next line and dev.off() to produce publication quality plots.
#    pdf(paste0("testWCCsinesB", i,".pdf"), height=7,width=9)    
    wccPlot(inSeries1=series1matrix[i,], inSeries2=series2matrix[i,], startwindow=startwindow, endwindow=endwindow, wMax=WCCwMax, tMax=WCCtMax, wInc=WCCwInc, tInc=WCCtInc, 
        Lsize=PPLsize, pspan=PPpspan, type=PPtype, samplespersecond=samplespersecond, windcross=FALSE)
#    dev.off()
}




sinOut20 <- wccFindDyadParam(inArray1=series1matrix, inArray2=series2matrix, wIncvector=WCCwIncvector, wMaxvector=WCCwMaxvector, tMaxvector=WCCtMaxvector, tIncvector=WCCtIncvector,
                             Lsizevector=PPLsizevector, pspanvector=PPpspanvector, type=testPPtype, nSurrogates=totalSurrogates, samplespersecond=samplespersecond, windcross=TRUE) 

wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$maxMeanQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)
wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$maxVarQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)
wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$totalMeanQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)
wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$totalVarQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)
wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$zeroLagMeanQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)
wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$zeroLagVarQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)
wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$lagMeanQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)
wccHeatmap(xparam=sinOut20$wMax, yparam=sinOut20$tMax, aggstat=sinOut20$lagVarQdiff, xlabel="Window Size", ylabel="Maximum Lag", pdffile=NA)


Sys.time()-startTime

# ----------------------------------
# Quit.

