\name{gaelle}

\alias{gaelle}

\docType{data}

\title{Messerli et. al. metabolomic data}

\description{The dataset consists of replicates of 14 mustants to study
metabolite pathways of genetically modified samples of \emph{Arabidopsis thaliana}. Values
of 43 metabolites are measured for each sample which are supposed to monitor their
genetic changes. The data involve two mutants defective in starch
biosynthesis, \emph{pgm} and \emph{isa2}; four defective in starch
degradation \emph{sex1}, \emph{sex4}, \emph{mex1}, and
\emph{dpe2}; a mutant for comparison that accumulates starch as a
pleitropic effect, \emph{tpt}; four uncharacterised mutants,
\emph{deg172}, \emph{deg263}, \emph{ke103}, and \emph{sex3};
and three wild type plants, \emph{WsWT}, \emph{RLDWT}, and
\emph{ColWT}. There are four replicates of all samples except the
last for which there are three \url{http://www.plantphysiol.org/cgi/content/abstract/143/4/1484}. 
}


\usage{data(gaelle)}
\format{
Matrix with 55 observations measured on 43 variables.
}

\references{
Messerli, G., Partovi Nia, V., Trevisan, M., Kolbe, A., Schauer, N., Geigen-berger, P., Chen, J., Davison, A. C., Fernie, A. R. and Zeeman, S. C. (2007) Rapid classification of phenotypic mutants of Arabidopsis via metabolite fingerprinting. Plant Physiology 143, 1481-1492.
}

\examples{
data(gaelle)
heatmap(gaelle)
}

