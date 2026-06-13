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
# Program: HeatMaps.R
#  Author: Steven Boker
#    Date: Thu Feb 12 11:07:51 EST 2026
#
# This program plots a heatmap of the proportions returned by wccFindDyadParam
#
# ---------------------------------------------------------------------
# Revision History
#  Steve Boker  -- Thu Feb 12 11:07:58 EST 2026
#      Created HeatMaps.R
#
# ---------------------------------------------------------------------


# ----------------------------------
# Read libraries and set options.

library(pheatmap)
library(gtable)


# ----------------------------------
# function to plot one tested parameter against another for an arbitrary aggregate statistic

wccHeatmap <- function(xparam=NA, yparam=NA, aggstat=NA, xlabel=NA, ylabel=NA, pdffile=NA) {
    if (!is.numeric(xparam) | !is.vector(xparam) | length(xparam) < 4  | !is.numeric(yparam) | !is.vector(yparam) | length(yparam) < 4  ) {
        stop(paste0("Warning: xparam and yparam must be numeric vectors with length greater or equal to 4."))
    }
    if (!is.numeric(aggstat) | !is.vector(aggstat) | length(aggstat) < 4 ) {
        stop(paste0("Warning: aggstat must be a numeric vector with length greater or equal to 4."))
    }
    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar)) # code line i + 1
    
    tXlabel <- sort(unique(xparam))
    tYlabel <- sort(unique(yparam))
    tMatrix <- matrix(NA, nrow=length(tYlabel), ncol=length(tXlabel))
    i <- 1
    for (tx in tXlabel) {
        j <- 1
        for (ty in tYlabel) {
            tMatrix[j,i] <- mean(aggstat[yparam==ty & xparam==tx])
            j <- j+1
        }
        i <- i+1
    }
    dimnames(tMatrix) <- list(tYlabel,tXlabel)
    # ----------------------------------
    # Plot it.

    if (!is.na(pdffile))     pdf(pdffile, height=5,width=5)
    p <- pheatmap(silent=TRUE,height=3,width=3,
      tMatrix,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      scale = "none"
    )

    grid.newpage()
    pushViewport(viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        heights = unit(c(1, 10), "null"),
        widths  = unit(c(1, 10), "null")
      )
    ))
    
    # Y-axis label
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    if(!is.na(ylabel))    grid.text(ylabel, rot = 90)
    else  grid.text("   ", rot = 90)
    popViewport()

    # X-axis label
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    if(!is.na(xlabel))    grid.text(xlabel,x=0.4)
    else  grid.text("   ", x=0.4)
    popViewport()

    # Heatmap
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
    grid.draw(p$gtable)
    popViewport(2)
    if (!is.na(pdffile))     dev.off()
    
}





