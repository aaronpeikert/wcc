# ---------------------------------------------------------------------
# Program: testWccPeakPick.R
#
# Verifies that the matrix-smoother rewrite of wccPeakPick reproduces
# the original per-row loess + spline implementation exactly (up to
# floating point), on WCC grids from the testSine simulation.
# ---------------------------------------------------------------------

library(wcc)

set.seed(260224)

# The original implementation, embedded verbatim as the reference
# (plotting code removed; it was dead, graphs was hard-coded to 0).
wccPeakPickLegacy <- function(tAllCor=NA, Lsize=8, pspan=.25, type="Max") {
    colLen <- length(tAllCor[1,])
    rowLen <- length(tAllCor[,1])
    mx <- rep(NA, (2*colLen-1))
    tIndex <- rep(NA, rowLen)
    tValue <- rep(NA, rowLen)
    findMin <- (type=="Min" || type=="min")
    for (rowNo in c(1:rowLen)) {
        if (any(is.na(tAllCor[rowNo, ]))) next
        tCor <- tAllCor[rowNo, ]
        t1 <- predict(loess(tCor~c(1:colLen), degree=2, span=pspan))
        t2 <- spline(c(1:colLen), t1, n=(2*colLen-1))$y
        windowWidth <- 0
        lookAhead <- 0
        for (j in 1:(colLen-1)) {
            windowWidth <- windowWidth + 1
            tSelect <- (colLen - windowWidth):(colLen + windowWidth)
            mx[j] <- if (findMin) min(t2[tSelect], na.rm=TRUE) else max(t2[tSelect], na.rm=TRUE)
            if (j == 1) mmx <- mx[j]
            else {
                better <- if (findMin) (mx[j] < mmx) else (mx[j] > mmx)
                if (better) {
                    lookAhead <- 0
                    mmx <- mx[j]
                }
                else {
                    lookAhead <- lookAhead + 1
                    if (lookAhead >= Lsize) break
                }
            }
        }
        Index <- match(mmx, t2[tSelect]) + tSelect[1] - 1
        position <- Index - colLen
        if (position > (colLen - Lsize - 1) || position < (-(colLen - Lsize - 1))) {
            tIndex[rowNo] <- NA
            tValue[rowNo] <- NA
        }
        else {
            tIndex[rowNo] <- position
            tValue[rowNo] <- mmx
        }
    }
    if (findMin) return(list(minIndex=tIndex, minValue=tValue))
    return(list(maxIndex=tIndex, maxValue=tValue))
}

makeSeries <- function(n) {
    sin(c(1:n)/runif(1, min=5, max=20)) + rnorm(n, mean=0, sd=.5)
}

sameOrBothNA <- function(a, b, tol=1e-10) {
    all(is.na(a) == is.na(b)) && all(abs(a - b) <= tol, na.rm=TRUE)
}

n <- 1000
for (d in 1:2) {
    s1 <- makeSeries(n)
    s2 <- makeSeries(n)
    for (wMax in c(50, 100)) {
        for (tMax in c(50, 100)) {
            grid <- wccCalc(s1, s2, wMax=wMax, tMax=tMax)
            grid[is.na(grid)] <- 0
            for (type in c("Max", "Min")) {
                for (pspan in c(.25, .4)) {
                    new <- wccPeakPick(grid, Lsize=8, pspan=pspan, type=type)
                    old <- wccPeakPickLegacy(grid, Lsize=8, pspan=pspan, type=type)
                    stopifnot(sameOrBothNA(new[[1]], old[[1]], tol=0),
                              sameOrBothNA(new[[2]], old[[2]]))
                }
            }
            cat(sprintf("dyad %d wMax=%3d tMax=%3d: new == legacy\n", d, wMax, tMax))
        }
    }
}

# Rows containing NA must yield NA index and value.
gridNA <- wccCalc(makeSeries(n), makeSeries(n), wMax=50, tMax=50)
gridNA[is.na(gridNA)] <- 0
gridNA[5, 17] <- NA
res <- wccPeakPick(gridNA, Lsize=8, pspan=.25, type="Max")
stopifnot(is.na(res$maxIndex[5]), is.na(res$maxValue[5]),
          !is.na(res$maxIndex[6]))

# Cache: second build of the same smoother matrix must come from the cache.
key <- paste0(ncol(gridNA), "_", 0.25)
stopifnot(!is.null(wcc:::.wccPeakPickCache[[key]]))
t1 <- system.time(wcc:::wccSmoothMatrix(ncol(gridNA), 0.25))["elapsed"]
stopifnot(t1 < 0.1)

cat("All wccPeakPick tests passed.\n")
