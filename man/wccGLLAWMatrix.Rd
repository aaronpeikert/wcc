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

\name{wccGLLAWMatrix}
\alias{wccGLLAWMatrix}

\title{W Matrix}

\description{
This function returns a matrix, W, which when premultiplied by a time delay embedded matrix will calculate the time derivatives of the time series in the time delay embedded matrix.
}

\usage{
wccGLLAWMatrix(embed=NA, tau=NA, deltaT=1, order=2)
}

\arguments{
\item{embed}{A numeric wholenumber value indicating the embedding dimension, the number of columns in the target time delay embedded matrix.  embed must be greater than or equal to order + 1.  \sQuote{embed} must be equal to that given to \sQuote{wccGLLAEmbed}  to create the target time delay embedded matrix.}
\item{tau}{An numeric wholenumber value indicating the number of samples by which each column in the time delay embedded matrix is lagged. \sQuote{tau}  must be greater than or equal to 1.   \sQuote{tau} must be equal to that given to \sQuote{wccGLLAEmbed}  to create the target time delay embedded matrix.}
\item{deltaT}{A numeric value indicating the elapsed time between samples in the time series that was time delay embedded.  deltaT must be greater than 0.}
\item{order}{A numeric wholenumber value indicating the highest derivative to be estimated when the W matrix is premultiplied by the time delay embedding matrix.}
}

\details{
The W matrix is defined by the Generalized Local Linear Approximation method for calculating derivatives of time series.  First time delay embed the time series using \sQuote{wccGLLAEmbed}  .  Then call \sQuote{wccGLLAWMatrix}   to return the appropriate W matrix for that embedding.  Finally postmultiply the time delay embedding matrix by the W matrix, which creates an output matrix whose rows match the rows of the time delay embedded matrix and whose columns are the calculated derivatives of the time series.
}
    
\value{
Returns the W matrix that can postmultiply a time delay embedded matrix to calculate derivatives.
}

\references{
Boker, S. M., Deboeck, P. R., Edler, C., & Keel, P. K. (2010) Generalized Local Linear Approximation of Derivatives from Time Series. In Statistical Methods for Modeling Human Dynamics: An Interdisciplinary Dialogue, S.-M. Chow & E. Ferrar (Eds). Boca Raton, FL: Taylor & Francis, pp 179–212.}

\seealso{
\link{wccGLLAEmbed} for creating a time delay embedding matrix. 
}

\examples{
# Create a time delay embedding matrix from a single time series.
theSeries <- sin((1:20)/10)
length(theSeries)
embeddedMatrix <- wccGLLAEmbed(theSeries, embed=4, tau=1, label="x", idColumn=FALSE)
dim(embeddedMatrix)
# Calculate derivatives for the series
embeddedMatrix \%*\% wccGLLAWMatrix(embed=4, tau=1, deltaT=1, order=2)

}
