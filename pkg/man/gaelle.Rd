\name{gaelle}

\alias{gaelle}

\docType{data}

\title{Messerli et al metabolomic data}

\description{ Data consists of replicates of 14 mustants to study
metabolite pathways of genetically modified samples of Arbidopsis
thaliana. Each mutant has 4 replicates but one has 3. }


\usage{data(gaelle)}
\format{
The format is:
matrix with 55 observations and measured on 43 variables.
}

\references{
http://www.plantphysiol.org/cgi/content/abstract/143/4/1484
}

\examples{
data(gaelle)
heatmap(gaelle)
}

