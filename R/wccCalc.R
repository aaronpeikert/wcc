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
# Program: wccCalc.R
#  Author: Steven Boker
#    Date: Wed Jan 28 08:52:17 EST 2026
#
# This function runs a windowed cross correlation and returns a matrix whose columns
# are the lags and whose rows are the elapsed time of each window.
#
# This either calls windcross or a version of windcross in R.
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Wed Jan 28 08:52:19 EST 2026
#      Created wccCalc.R
#
# ---------------------------------------------------------------------

# ----------------------------------
# Calculate WCC.  

wccCalc <- function(inSeries1, inSeries2, wMax=50, tMax=50, wInc=1, tInc=1, windcross=TRUE) {
    if (!is.numeric(inSeries1) | !is.numeric(inSeries2) | !is.vector(inSeries1) | !is.vector(inSeries2) | length(inSeries1) != length(inSeries2)) {
        stop(paste0("Warning: inSeries1 and inSeries2 must be numeric vectors of equal length."))
    }
    if (!is.numeric(wMax) | wMax < 5  ) {
        stop(paste0("Warning: wMax must be a numeric greater than or equal to 5."))
    }
    if ( !is.numeric(wInc) |  wInc < 1 ) {
        stop(paste0("Warning: wInc must be a numeric greater than or equal to 1."))
    }
    
    maxRowLen <- length(inSeries1)
    tStart <- wMax + tMax
	nRow <- floor((maxRowLen - tStart) / wInc)
	nCol <- 2*(floor(tMax / tInc)) + 1
    centerCol <- (floor(tMax / tInc)) + 1
    if (nRow < 0) {
        stop(paste0("Warning: bad choice for wMax and/or wInc parameters. The result matrix has ", nRow, " rows."))
    }
    if (nCol < 0) {
        stop(paste0("Warning: bad choice for tMax and/or tInc  parameters. The result matrix has ", nCol, " columns."))
    }
    if (windcross==TRUE) {
        inSeries1[is.na(inSeries1)] <- -999999.0
        inSeries2[is.na(inSeries2)] <- -999999.0
        tData <- .Call("windcross", as.numeric(inSeries1), as.numeric(inSeries2), as.numeric(wMax), as.numeric(tMax), as.numeric(wInc), as.numeric(tInc), PACKAGE = "wcc")
        tData[tData < -99] <- NA
        return(tData)
    }
    else {
        tData <- matrix(NA, nrow=nRow, ncol=nCol)
        elapsedSeq <-  seq(from=tStart, to=maxRowLen, by=wInc)
        for (tRow in 1:nRow) {
            elapsedIndex <- elapsedSeq[tRow]
            lag0Sel <- seq(from=elapsedIndex, to=(1+ elapsedIndex-wMax), by=-1)
            lag0win1 <- inSeries1[lag0Sel]
            lag0win2 <- inSeries2[lag0Sel]
            tData[tRow,centerCol] <- cor(lag0win1,lag0win2, use="pairwise.complete.obs")
            for(i in 1:(floor(tMax / tInc))) {
                tData[tRow,(centerCol+i)] <- cor(lag0win1, inSeries2[lag0Sel-i], use="pairwise.complete.obs")
                tData[tRow,(centerCol-i)] <- cor(lag0win2, inSeries1[lag0Sel-i], use="pairwise.complete.obs")
            }
        }
        return(tData)
    }
}
