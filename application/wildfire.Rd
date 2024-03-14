\encoding{UTF-8}
\name{Wildfire}
\alias{Wildfire}
\title{Portugal Wildfire Data}
\description{
  The \code{wildfire} data from Portugal consist of daily burnt area (in hectares) for each wildfire and Canadian Forest Fire Weather Index System indices between 1980 to 2019.
}

\usage{wildfire}
\format{
  \code{wildfire} is a data frame with 14,609 occurances (rows) and 11 variables (columns).

  The \code{wildfire} data frame contains the following columns:
  \itemize{
    \item{\code{Burnt_Area}}{:
	    daily burnt area in hectares.
    }
    \item{\code{DSR}}{: 
      Daily Severity Rating (DSR), a numeric rating of the difficulty of controlling fires. 
    }
    \item{\code{FWI}}{: 
      Fire Weather Index (FWI), a numeric rating of fire intensity. 
    }
    \item{\code{BUI}}{: 
      Buildup Index (BUI), a numeric rating of the total amount of fuel available for combustion. 
    }
    \item{\code{ISI}}{: 
      Initial Spread Index (ISI), a numeric rating of the expected rate of fire spread. 
    }
    \item{\code{FFMC}}{: 
      Fine Fuel Moisture Code (FFMC), a numeric rating of the moisture content of litter and other cured fine fuels. 
    }
    \item{\code{DMC}}{: 
      Duff Moisture Code (DMC), a numeric rating of the average moisture content of loosely compacted organic layers of moderate depth. 
    }    
    \item{\code{DC}}{: 
      Drought Code (DC), a rating of the average moisture content of deep, compact organic layers. 
    }
    \item{\code{day}, \code{month}, \code{year}}{:
      timestamp to date for each datapoints.
    }
  }
}
\source{Instituto Dom Luiz}
\references{
  Lee, M. W., de Carvalho, M., Paulin, D., Pereira, S., Trigo, R. and Da Camara, C. (2024)
  A Bayesian Lasso for Tail Index Regression. \emph{Submitted}.
}

\examples{
## preview of the data
data("wildfire")
head(wildfire, 10)
summary(wildfire)

require(ggplot2)
## visualizing the data by month
ggplot(wildfire, aes(x = month, y = Burnt_Area, color = month)) + 
     geom_point(size = 6) + ylab("Burnt Area (ha)")+ xlab("Month")
     theme_minimal() +
     theme(plot.title = element_text(hjust = 0.5, size = 20),
            legend.title = element_blank(),
            legend.text = element_text(size = 20),
            axis.text.x = element_text(hjust = 0.35),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 15))
}
\keyword{wildfire}
