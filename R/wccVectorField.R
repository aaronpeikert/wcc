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
# Program: wccVectorField.R
#  Author: Steven Boker
#    Date: Wed Apr 15 13:07:19 EDT 2026
#
# This program calculates a best forward slope for every cell in a wccCalc matrix
#
#   You must previously have run wccCalc to create a wccCalc matrix
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Wed Apr 15 13:07:20 EDT 2026
#      Created wccVectorField.R
#
# ---------------------------------------------------------------------

iwccMaxDeriv <- function(tVec=NA, type="Zero") {  # Internal function to return the max, min or derivative nearest zero among the 8 directions
    if (type=="Min")    
        return(order(tVec)[1])
    if (type=="Zero")    
        return(order(abs(tVec))[1])
    return(order(tVec)[8])
}

# ----------------------------------
# Calculate the derivative with the maximum, minimum or value nearest zero
#  Return the value of the derivative and the angle in radians where this derivative points.

wccVectorFieldCalc <- function(tAllCor=NA, type="Zero", kernelWidth=4) { 
    nRows <- dim(tAllCor)[1] - (kernelWidth*2)
    nCols <- dim(tAllCor)[2]
    theSlopes <- matrix(NA, nrow=nRows, ncol=nCols-(kernelWidth*2))
    theAngles <- matrix(NA, nrow=nRows, ncol=nCols-(kernelWidth*2))
    W1 <- wccGLLAWMatrix(embed=(kernelWidth*2)+1, tau=1, deltaT=1, order=2)
    tAngles <- matrix(c(0,pi/2,atan(1),atan(-1),atan(1/2),atan(-1/2),atan(2),atan(-2)), nrow=nCols-2, ncol=8, byrow=TRUE)
    for (tRow in 1:nRows) {
        #  Embedding in the direction of elapsed time (angle=0)
        embeddedMatrix <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        tSel <- (kernelWidth+1):(nCols-kernelWidth)
        embeddedMatrix[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        stepMax <- 100
        stepCount <- stepMax
        stepIndex <- 0
        for (i in 1:kernelWidth) {
            embeddedMatrix[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth-i), tSel-stepIndex]
            embeddedMatrix[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth+i), tSel+stepIndex]
            stepCount <- stepCount-1
            if (stepCount == 0) {
                stepIndex <- stepIndex+1
                stepCount <- stepMax
            }
        }
        tSlope1 <- (embeddedMatrix %*% W1)[,2] # derivatives in the direction of elapsed time (angle=0)

        #  Embedding in the direction of elapsed time (angle=pi/2)
        embeddedMatrix <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        tSel <- (kernelWidth+1):(nCols-kernelWidth)
        embeddedMatrix[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        stepMax <- 1
        stepCount <- stepMax
        stepIndex <- 1
        for (i in 1:kernelWidth) {
            embeddedMatrix[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth), tSel-stepIndex]
            embeddedMatrix[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth), tSel+stepIndex]
            stepCount <- stepCount-1
            if (stepCount == 0) {
                stepIndex <- stepIndex+1
                stepCount <- stepMax
            }
        }
        tSlope2 <- (embeddedMatrix %*% W1)[,2] # derivatives in the direction of elapsed time (angle=pi/2)

        #  Embedding in the direction of positive  and negative lag (angle=pi/4 and -pi/4)
        embeddedMatrix <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        embeddedMatrix2 <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        tSel <- (kernelWidth+1):(nCols-kernelWidth)
        embeddedMatrix[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        embeddedMatrix2[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        stepMax <- 1
        stepCount <- stepMax
        stepIndex <- 1
        for (i in 1:kernelWidth) {
            embeddedMatrix[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth-i), tSel-stepIndex]
            embeddedMatrix[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth+i), tSel+stepIndex]
            embeddedMatrix2[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth-i), tSel-stepIndex]
            embeddedMatrix2[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth+i), tSel+stepIndex]
            stepCount <- stepCount-1
            if (stepCount == 0) {
                stepIndex <- stepIndex+1
                stepCount <- stepMax
            }
        }
        tSlope3 <- (embeddedMatrix %*% W1)[,2] # derivatives in the direction of positive lag (angle=pi/4)
        tSlope4 <- (embeddedMatrix2 %*% W1)[,2] # derivatives in the direction of positive lag (angle=-pi/4)

        #  Embedding in the direction of positive  and negative lag (angle=-atan(1/2) and atan(-1/2) )
        embeddedMatrix <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        embeddedMatrix2 <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        tSel <- (kernelWidth+1):(nCols-kernelWidth)
        embeddedMatrix[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        embeddedMatrix2[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        stepMax <- 2
        stepCount <- 1
        stepIndex <- 0
        for (i in 1:kernelWidth) {
            embeddedMatrix[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth-i), tSel-stepIndex]
            embeddedMatrix[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth+i), tSel+stepIndex]
            embeddedMatrix2[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth-i), tSel-stepIndex]
            embeddedMatrix2[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth+i), tSel+stepIndex]
            stepCount <- stepCount-1
            if (stepCount == 0) {
                stepIndex <- stepIndex+1
                stepCount <- stepMax
            }
        }
        tSlope5 <- (embeddedMatrix %*% W1)[,2] # derivatives in the direction of positive lag (angle=atan(1/2))
        tSlope6 <- (embeddedMatrix2 %*% W1)[,2] # derivatives in the direction of positive lag (angle=atan(-1/2))

        #  Embedding in the direction of positive  and negative lag (angle=-atan(2) and atan(-2) )
        embeddedMatrix <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        embeddedMatrix2 <- matrix(NA, nrow=(nCols-(kernelWidth*2)), ncol=(kernelWidth*2)+1)
        tSel <- (kernelWidth+1):(nCols-kernelWidth)
        embeddedMatrix[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        embeddedMatrix2[,kernelWidth+1] <- tAllCor[tRow+kernelWidth, tSel]
        stepMax <- 2
        stepCount <- 1
        stepIndex <- 0
        for (i in 1:kernelWidth) {
            embeddedMatrix[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth+stepIndex), tSel-i]
            embeddedMatrix[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth-stepIndex), tSel+i]
            embeddedMatrix2[,(kernelWidth+1)+i] <- tAllCor[tRow+(kernelWidth-stepIndex), tSel-stepIndex]
            embeddedMatrix2[,(kernelWidth+1)-i] <- tAllCor[tRow+(kernelWidth+stepIndex), tSel+stepIndex]
            stepCount <- stepCount-1
            if (stepCount == 0) {
                stepIndex <- stepIndex+1
                stepCount <- stepMax
            }
        }
        tSlope7 <- (embeddedMatrix %*% W1)[,2] # derivatives in the direction of positive lag (angle=atan(2))
        tSlope8 <- (embeddedMatrix2 %*% W1)[,2] # derivatives in the direction of positive lag (angle=atan(-2))

        t0 <- cbind(tSlope1,tSlope2,tSlope3,tSlope4,tSlope5,tSlope6,tSlope7,tSlope8)
        t1 <- apply(t0,1,iwccMaxDeriv,type=type)
        theSlopes[tRow,] <- t0[cbind(1:nrow(t0),t1)]
        theAngles[tRow,] <- tAngles[cbind(1:nrow(t0),t1)]
    }
    
    return(list(slope=theSlopes,angle=theAngles))
}

# ----------------------------------
# Plot a vector field or topopgraphic map optionally over the top of a WCC plot

wccVectorFieldPlot <- function(VFframe=NA, startwindow=1, endwindow=200, wMax=50, tMax=50, wInc=1, tInc=1, kernelWidth=4, samplespersecond=1, type="Zero") {
    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar)) # code line i + 1
    minx <- wMax + tMax + 2*kernelWidth + 1 + (startwindow * wInc / samplespersecond) 
    maxx <- wMax + tMax + 2*kernelWidth + 1 + (endwindow * wInc / samplespersecond)
    plot(c(minx, maxx), c(-tMax,tMax), type='n', xlab="Elapsed Time (seconds)", ylab="Lag (seconds)")
    for(i in seq(1,dim(VFframe$angle)[1],by=2)) {
        for (j in seq(1,dim(VFframe$angle)[2],by=2)) {
            x <- minx + (i * wInc)
            y <- tMax - (j * tInc)
            lines(c(x,(x+2*cos(VFframe$angle[i,j]))),c(y,(y-2*sin(VFframe$angle[i,j]))),type='l',lwd=.5)
        }
    }
}

# ----------------------------------
# 

