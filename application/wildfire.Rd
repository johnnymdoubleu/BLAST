\encoding{UTF-8}
\name{Wildfire}
\alias{Wildfre}
\title{Portugal Wildfire Data}
\description{
  The \code{wildfire} data consist of daily burnt area (in hectares) for each wildfire and Canadian Forest Fire Weather Index System (van Wagner et al., 1974) indices of Portugal between 1980 and 2019.
}

\usage{wildfire}
\format{
  The \code{wildfire} dataframe contains the following columns:
  \itemize{
    \item{\code{Burnt_Area}}{:
	    daily burnt area in hectares.
    }
    \item{\code{DSR}}{: 
      Daily Severity Rating (DSR), a numeric rating of the difficulty of controlling fires. 
    }
    \item{\code{FWI}}{: 
      Fire Weather Index (FWI), rating of fire intensity. 
    }
    \item{\code{BUI}}{: 
      Buildup Index (BUI), rating of the total amount of fuel available for combustion. 
    }
    \item{\code{ISI}}{: 
      Initial Spread Index (ISI), rating of the expected rate of fire spread. 
    }
    \item{\code{FFMC}}{: 
      Fine Fuel Moisture Code (FFMC), rating of the moisture content of litter and other cured fine fuels. 
    }
    \item{\code{DMC}}{: 
      Duff Moisture Code (DMC), rating of the average moisture content of loosely compacted organic layers of moderate depth. 
    }    
    \item{\code{DC}}{: 
      Drought Code (DC), a rating of the average moisture content of deep, compact organic layers. 
    }
    \item{\code{day}, \code{month}, \code{year}}{:
      time stamp of each daily datapoints.
    }
  }
}
\source{Instituto Dom Luiz}
\references{
  Lee, M. W., de Carvalho, M., Paulin, D., Pereira, S., Trigo, R. and Da Camara, C. (2024)
  A Bayesian Lasso for Tail Index Regression.
  \emph{-}, \bold{-}, --.
  
  van Wagner, C. E. and others. (1974)
  Structure of the Canadian Forest Fire Weather Index.
  vol. \bold{1333}, \emph{Environment Canada, Forestry Service Ottawa}.  
}

\examples{
## preview of the data
data("wildfire")
head(wildfire, 10)

## visualizing the data
}
\keyword{wildfire}
