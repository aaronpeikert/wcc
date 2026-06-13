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
# Program: peakPick.R
#  Author: Steven Boker, Minquan Xu
#    Date: Wed Jan 28 08:55:03 EST 2026
#
# This function takes a pair of timeseries runs wcc and peak picking and returns 
#   a wccAgg dataframe with aggregated statistics.  Optionally it saves the wcc matrix
#   and peakpicking vector to a csv file.
#
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- 2001
#      Created peakPick.R
#
# ---------------------------------------------------------------------

# input data structure(Splus object matrix):
#    row: time
#    column: cross correlations
# output: a list of local peak indice and values

# Parameters:
# ---------------------------------------------------------------------
# tAllCor: Splus object matrix, this matrix is created from the output
#           file of windcross program using Splus function such as
#           "scan". ex. tAllCor <- matrix(scan("windcross.dat"), ncol=n, byrow=T), 
#           where n is the number of columns in windcross.dat file. 
# Lsize: local search region, the value should be larger than 0 and 
#        less than 1/2 length of one row
# graphs: number of graphs to draw for one input data object, the
#         value should be larger than 0 and less than number of rows
#         of one input data object
# pspan: see loess() for span
# type: local maximum or local minimum, valid values are: "Min" and "Max" 
# tFileName: root characters for .pdf file
#---------------------------------------------------------------------

wccPeakPick <- function(tAllCor=NA, Lsize=8, pspan=.25, type="Max") { 
    graphs <- 0
    tFileName <- "peak"
    #----------------check for validness of parameters -------------------
    colLen <- length(tAllCor[1,]) # col length --- number of columns 
    rowLen <- length(tAllCor[,1]) # row length --- number of rows
    tLsize <- floor((1/2)*colLen) # maximum local search region

    if(Lsize<1 || Lsize>tLsize) { # Lsize too small or too large
          errorStr<- paste("Lsize should be >0 and <= ", tLsize, sep="")
          stop(errorStr) # print error message and stop the program
    }
    if(graphs<0||graphs>rowLen) { # num of graphics to print is too small or large
           errorStr <- paste("graphs should be >=0 and <= ", rowLen, sep="")
           stop(errorStr) # print error message and stop the program
    }
    if(pspan<0 || pspan>1) { # invalid pspan value
           stop("pspan should be >0 and <1\n") # print error message and stop program
    }
    if(type!="Max"&&type!="Min"&&type!="max"&&type!="min"){ # only two types
           stop("valid types are: max|Max or Min|min \n") # print error message and stop
    }
    #-----------Initilization------------------------------------------------
    drawgraph <- 0 # graphics drawn
    colLen <- length(tAllCor[1,]) # col length  
    rowLen <- length(tAllCor[,1]) # row length 
    xSequence <- seq(-(colLen-1), (colLen-1), by=1) # X axis for each graph
    mx<- rep(NA, (2*colLen-1)) #vector for keeping temp peak value for a row
                               #data points will be 2*colLen-1 after smooth
    tIndex <- rep(NA, rowLen) #vector of peak index---one peak index for each row
    tValue <- rep(NA, rowLen) #vector of peak value---one peak value for each row 

    #------------- type == max or Max ---------------------------------------
    if(type=="Max"||type=="max") { #compute local maximum
          for(rowNo in c(1: rowLen)) { #acess each row
                #eliminate missing value
                miss <- is.na(tAllCor[rowNo, ]) #transfer a row to be T--"NA"  and F not "NA"
                #initialize the position of NA in a row 
                missposition <- 0 
                for(mIndex in c(1:colLen)) { #evaluate an entire row
                      missposition <- mIndex #the position of NA
                      if(miss[missposition]) break #find one
                      missposition <- missposition+1 #increase count
                      if(missposition==colLen+1) break #No NA in this row
                }
                #cat("missposition=", missposition, "\n")
                #if has missing value 
                if(missposition <= colLen) next #skip a row with NA
           
                else { # no missing value
                      drawgraph <- drawgraph+1 #number of graph to draw
                      tCor <- tAllCor[rowNo, ] #number of columns
                      #smooth
                      t1 <- predict(loess(tCor~c(1:colLen), degree=2, span=pspan, ))
                      #data points is set to n
                      t2 <- spline(c(1:colLen), t1, n=(2*colLen-1))$y

                      # show calculate progress
                      # cat("row=", rowNo, "\n")
                      #----------- process a row --find max value and max index------ 
                      windowWidth <- 0  # searched region
                      lookAhead <- 0  # look ahead data points
                      for(j in 1:(colLen-1)) {  # search from 1 to colLen -1
                            windowWidth <- windowWidth+1 # increase search ed region
                            # select the search region, the center of search region 
                            # is in the middle of t2, notice that t2 has 2*colLen-1
                            # data points.
                            tSelect <- (colLen - windowWidth):(colLen+windowWidth)
                            mx[j] <- max(t2[tSelect], na.rm=T) # store temp max value
                            if(j==1) mmx <- mx[j] # mmx is final local max value
                            #remember current max
                            else { # if j != 1
                                  if(mx[j]>mmx) { # new temp max value , note only one                                                     # max value in tSelect
                                        lookAhead <- 0 # set stop search criterion to 0,
                                                 # the criterion is that if we find
                                                 # Lsize data points less than current
                                                 # max, then stop and the current max 
                                                 # value is the local maximum we wanted
                                         mmx <- mx[j] # update new max value
                                  }
                            else if(mx[j]<=mmx) { # if other values are less than
                                  # current maximum
                                  # increase the count---how many neighbor data                                            # point are less than current maximum 
                                  lookAhead <- lookAhead+1 
                                  if(lookAhead>=Lsize) break # meet criterion
                            }
                            }#else
                      }#for j--max value and index for each row
                      #use match function to find the index
                      Index <- match(mmx, t2[tSelect])+tSelect[1]-1
                      # tSelect[1] is the first index of the selected window

                      # relative position to the middle point
                      position <- Index -colLen
                      
                      #according to the local maximium definition 
                      if(position >(colLen - Lsize - 1) || position < (-(colLen - Lsize -1))) { # fail
                            tIndex[rowNo] <- NA
                            tValue[rowNo] <- NA
                      }
                      else { # found a local maximum
                            tIndex[rowNo] <- position
                            tValue[rowNo] <-mmx
                      }

                      #draw plots      
                      if(drawgraph <= graphs) {
                           # define graphic file name
                           tepsfile <- paste(tFileName, "Max", rowNo, ".eps", sep="") 
                           # title of the graph
                           tmain <- paste("max Index", tFileName, "r", rowNo,"w", Lsize, sep="")
                           # write to postscript format
                           postscript(tepsfile, height=6.4, horizontal=F)
           
                           # draw borders and their labels
                           plot(c(-(colLen-1), (colLen-1)), c(-1,1), xlab="Lag", ylab="Cross Correlation", main=tmain, type="n")
           
                           # draw the curve
                           lines(xSequence, t2, type="l")
           
                           # draw the local maximum
                           lines(c(position,position), c(-1,1), type="l", lty=8)
           
                           # draw the axies
                           lines(c(0,0), c(-1,1), type="l", lty=4)
                           lines(c(-(colLen-1), (colLen-1)), c(0,0), type="l", lty=4)
                           dev.off() # term off other device in order run drawing procedure
                      } # if drawgraph
                } #else no missing value
          } # for rowNo-- process each row 
    #end of process a row
    return(list(maxIndex=tIndex, maxValue=tValue))
    }#if type=max
    
    #-------------------------------------------------------------------------------
    # type == Min or min
    #-------------------------------------------------------------------------------
    else if(type=="Min"||type=="min") {
          for(rowNo in c(1: rowLen)) { 
                #eliminate missing value
                miss <- is.na(tAllCor[rowNo, ])
                missposition <- 0
                for(mIndex in c(1:colLen)) {
                      missposition <- mIndex #the position of NA
                      if(miss[missposition]) break #find one
                      missposition <- missposition+1
                      if(missposition==colLen+1) break
                }
                #cat("missposition=", missposition, "\n")
                #if has missing value 
                if(missposition <= colLen) next #skip a row with NA
           
                else { # no missing value
                      drawgraph <- drawgraph+1
                      tCor <- tAllCor[rowNo, ]
                      t1 <- predict(loess(tCor~c(1:colLen), degree=2, span=pspan, ))
                      t2 <- spline(c(1:colLen), t1, n=(2*colLen-1))$y
                      # show calculate progress
                      #cat("row=", rowNo, "\n")

                      #process a row --find max value and max index 
                      #-------------------------------------------------------------
                      windowWidth <- 0
                      lookAhead <- 0
                      for(j in 1:(colLen-1)) { 
                          windowWidth <- windowWidth+1
                          tSelect <- (colLen - windowWidth):(colLen+windowWidth)
                          mx[j] <- min(t2[tSelect], na.rm=T)
                          if(j==1) mmx <- mx[j]
                          #remember current max
                          else {
                          if(mx[j]<mmx) { #only one value in tSelect
                                lookAhead <- 0
                                mmx <- mx[j]
                          }
                          else if(mx[j]>=mmx) {
                                lookAhead <- lookAhead+1
                                if(lookAhead>=Lsize) break
                          }
                          }#else
                      }#for j--max value and index for each row
                      #use match function to find the index
                      Index <- match(mmx, t2[tSelect])+tSelect[1]-1
                      position <- Index -colLen
 
                      if(position >(colLen - Lsize -1) || position < (-(colLen -Lsize -1))) { # fail
                            tIndex[rowNo] <- NA
                            tValue[rowNo] <- NA
                      }
                      else { 
                            tIndex[rowNo] <- position
                            tValue[rowNo] <-mmx
                      }

                      #draw first 10 plots      
                      if(drawgraph <= graphs) {
                            tepsfile <- paste(tFileName, "Min", rowNo, ".eps", sep="") 
                            tmain <- paste( "min Index", tFileName, "r", rowNo,"w", Lsize, sep="")
                            postscript(tepsfile, height=6.4, horizontal=F)
                            plot(c(-(colLen-1), (colLen-1)), c(-1,1), xlab="Lag", ylab="Cross Correlation", main=tmain, type="n")

                            lines(xSequence, t2, type="l")
                            lines(c(position,position), c(-1,1), type="l", lty=8)
                            lines(c(0,0), c(-1,1), type="l", lty=4)
                            lines(c(-(colLen-1), (colLen-1)), c(0,0), type="l", lty=4)
                            dev.off()
                      } # if drawgraph 
                } #else no missing value
          } # rowNo-- process each row 
     #end of process a row
     #-----------------------------------------------------------------------

    return(list(minIndex=tIndex, minValue=tValue))
    }#if type=mix
}

