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

#---------------------------------------------------------
# wccGLLAWMatrix -- Calculates a GLLA linear transformation matrix to 
#                create approximate derivatives
#
# Input:  embed -- Embedding dimension 
#           tau -- Time delay in samples
#        deltaT -- Interobservation interval in time units
#         order -- Highest order of derivatives (2, 3, or more)

wccGLLAWMatrix <- function(embed=NA, tau=NA, deltaT=1, order=2) {
    if (!is.numeric(embed) | !is.numeric(tau) | !is.numeric(deltaT) | !is.numeric(order)) {
        stop(paste0("Warning: embed, tau, and deltaT must be specified and must be numerics."))
    }
    if (embed < order+1) {
        stop(paste0("Warning: embed must be greater than or equal to order+1."))
    }
    if (tau < 1) {
        stop(paste0("Warning: tau must be greater than or equal 1."))
    }
    if (deltaT <= 0) {
        stop(paste0("Warning: deltaT must be greater than 0."))
    }
    if (order < 1 ) {
        stop(paste0("Warning: order must be greater than 0."))
    }
    L <- rep(1,embed)
    for(i in 1:order) {
        L <- cbind(L,(((c(1:embed)-mean(1:embed))*tau*deltaT)^i)/factorial(i)) 
    }
    return(L%*%solve(t(L)%*%L))
}


#---------------------------------------------------------
# wccGLLAEmbed -- Creates a time-delay embedding of a variable 
#              given a vector and an optional grouping variable
#              Requires equal interval occasion data ordered by occasion.
#              If multiple individuals, use the ID vector as "groupby"
#
# Input:      x -- vector to embed
#         embed -- Embedding dimension (2 creates an N by 2 embedded matrix)
#           tau -- rows by which to shift x to create each time delay column 
#       groupby -- grouping vector
#         label -- variable label for the columns
#      idColumn -- if TRUE, return ID values in column 1
#                  if FALSE, return the embedding columns only.
#
# Returns:  An embedded matrix where column 1 has the ID values, and the
#           remaining columns are time delay embedded according to the arguments.

wccGLLAEmbed <- function(x, embed=4, tau=1, groupby=NA, label="x", idColumn=TRUE) {
    
    if (!is.numeric(embed) | !is.numeric(tau) ) {
        stop(paste0("Warning: embed and tau must be whole numbers."))
    }
    if (embed < 2) {
        stop(paste0("Warning: embed must be a whole number greater than or equal to 2."))
    }
    if (!is.vector(x)) {
        stop(paste0("Warning: x must be a vector."))
    }
    if (tau < 1) {
        stop(paste0("Warning: tau must be greater than or equal 1."))
    }
    minLen <- (tau + 1 + ((embed - 2) * tau))
    if (!is.vector(groupby) | length(groupby[!is.na(groupby[])])<1) {
        groupby <- rep(1,length(x))
    }
    x <- x[!is.na(groupby[])]
    groupby <- groupby[!is.na(groupby[])]
    if (length(x) < minLen)
        stop(paste0("Warning: The length of x is shorter than the minimum length required by the specified embed and tau."))
    if (length(groupby) != length(x))
        stop(paste0("Warning: The length of x must be equal to the length of groupby."))
    embeddedMatrix <- matrix(NA, length(x) + (embed*tau), embed+1)
    colNames <- c("ID", paste(label, "0", sep=""))
    for (j in 2:embed) {
        colNames <- c(colNames, paste(label, (j-1)*tau, sep=""))
    }
    dimnames(embeddedMatrix) <- list(NULL, colNames)
    tRow <- 1
    for (i in unique(groupby)) {
        tx <- x[groupby==i]
        if (length(tx) < minLen)
            next
        tLen <- length(tx) - minLen
        embeddedMatrix[tRow:(tRow+tLen), 1] <- i
        for (j in 1:embed) {
            k <- 1 + ((j-1)*tau)
            embeddedMatrix[tRow:(tRow+tLen), j+1] <- tx[k:(k+tLen)]
        }
        tRow <- tRow + tLen + 1
    }
    if (idColumn==TRUE) {
        return(embeddedMatrix[1:(tRow-1),])
    }
    return(embeddedMatrix[1:(tRow-1), 2:(embed+1)])
}

#---------------------------------------------------------
# gllaPeriod -- Calculates the period of an oscillation from eta and zeta 

# wccGLLAPeriod <- function(eta, zeta) {2*pi/(sqrt(-(eta+((zeta^2)/4))))}

