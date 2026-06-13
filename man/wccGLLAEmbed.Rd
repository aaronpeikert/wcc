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

\name{wccGLLAEmbed}
\alias{wccGLLAEmbed}

\title{wccGLLAEmbed}

\description{
This function creates a time delay embedded matrix from a time series.
}

\usage{
wccGLLAEmbed(x, embed=4, tau=1, groupby=NA, label="x", idColumn=TRUE)
}

\arguments{
\item{x}{A numeric vector.}
\item{embed}{A numeric wholenumber value indicating the embedding dimension, the number of columns in the target time delay embedded matrix.  embed must be greater than or equal to order + 1}
\item{tau}{An numeric wholenumber value indicating the number of samples by which each column in the time delay embedded matrix is lagged. tau must be greater than or equal to 1.}
\item{groupby}{If \code{NA}, then do not group the rows of the output matrix.  If not missing, it is a numeric wholenumber vector of the same length as x.  Each element in the vector x is treated as belonging to the group specified by the wholenumber in the matching element of groupby.}
\item{label}{A character string specifying the label prefix for each column of the output matrix.  The number of the column is pasted onto this prefix so that columns have unique names.}
\item{idColumn}{A logical specifying whether to rbind the group as specified by groupby to the first column of the time delay embedded matrix.  This column is always named "ID"}
}

\details{
wccGLLAEmbed creates a time delay embedded matrix from a time series.  A time delay embedded matrix has rows corresponding to the samples of elapsed time in the time series and columns that are lagged by tau samples from each other.  For example, for a time series x of length N composed of a single group and with embed = 3 and tau=1, the time delay embedding matrix would have N-2 rows and 3 columns.  Column 1 would contain x[1:(N-2)], column 2 would contain x[2:(N-1)], and column 3 would contain x[3:N].  Thus each row contains a time lagged vector representing not only the state of the system at a moment in time, but also a snapshot of the dynamics of the system at that moment.  The number of columns in a N x D time delay embedding matrix is often referred to as the embedding dimension of the matrix, since each row can be thought of as a point in a D dimensional space.

When multiple groups' or individuals' time series are time delay embedded, one must prevent the data from group 1 to overlap on the same row with data from group 2 and so forth.  Thus, if a group has N elements, the resulting submatrix for that group will have embed columns and N - ((embed - 1) * tau) rows.  If any group has less than (tau + 1 + ((embed - 2) * tau)) elements then that group will not appear in the time delay embedding matrix since there is not sufficient data in the time series to compose 1 row.

Missing values are retained in the output time delay embedded matrix.  This is useful for methods such as latent differential equations where full information maximum likelihood can account for the missingness, but will result in missing valued rows in the output of Generalized Local Linear Approximation as calculated using the wccGLLAWMatrix function. 

Time delay embedding  is a technique that represents a time series in a state space that contains both the instantaneous state of a system and its local dynamics.  Takens (1981) showed that time delay embedding captures all of the dynamics in a system if the embedding dimension is sufficient.  Choosing the appropriate embedding dimension for representing the dynamics of a system is something of an art (Sauer, Yorke, and Casdagli, 1991).  For calculating time series deriviatives using Generalized Local Linear Approximation, the minimum embedding dimension is recommended to be two plus the order of the highest derivative to be calculated.  Higher embedding dimensions than required will perform smoothing on the time series, which may be useful, but will also attenuate the values of the derivatives calculated.

}
    
\value{
Returns a time delay embedding matrix.  Optionally this includes a column of ID numbers representing the group to which each row belongs.
}

\references{
Boker, S. M., Rotondo, J. L., Xu, M., & King, K. (2002). Windowed cross-correlation and peak picking for the analysis of variability in the association between behavioral time series. Psychological methods, 7(3), 338.

Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010) Generalized Local Linear Approximation of Derivatives from Time Series. In Statistical Methods for Modeling Human Dynamics: An Interdisciplinary Dialogue, S.-M. Chow & E. Ferrar (Eds). Boca Raton, FL: Taylor & Francis, pp 179–212.

Sauer, T., Yorke, J., & Casdagli, M. (1991). Embedology. Journal of Statistical Physics, 65(3,4), 95-116.

Takens, F. (1981). Detecting strange attractors in turbulence. In D. Rand & L.-S. Young(Eds.), Dynamical systems and turbulence, warwick 1980: Proceedings of a
symposium held at the University of Warwick (pp. 366-381). Berlin: Springer-Verlag.
}

\seealso{
\link{wccGLLAWMatrix} for the definition of the matrix required to estimate derivatives using Generalized Local Linear Approximation. 
}

\examples{
# Create a time delay embedding matrix from a single time series.
theSeries <- sin((1:20)/10)
length(theSeries)
embeddedMatrix <- wccGLLAEmbed(theSeries, embed=4, tau=1, label="x", idColumn=FALSE)
dim(embeddedMatrix)
# Calculate derivatives for the series
embeddedMatrix \%*\% wccGLLAWMatrix(embed=4, tau=1, deltaT=1, order=2)


theFrame <- data.frame(theSeries=c(sin((1:20)/10)), sin((1:20)/15),theGroups=c(rep(1,20),rep(2,20)))
dim(theFrame)
embeddedMatrix <- wccGLLAEmbed(theFrame$theSeries, embed=4, tau=1, label="x", 
     groupby=theFrame$theGroups, idColumn=TRUE)
dim(embeddedMatrix)
# Calculate derivatives for the series
embeddedMatrix[,2:5] \%*\% wccGLLAWMatrix(embed=4, tau=1, deltaT=1, order=2)

}
