\name{DSBayes}
\alias{DSBayes}
\title{Bayesian subgroup analysis in clinical trials}
\description{Calculate posterior modes and credible intervals of parameters of 
the Dixon-Simon model for subgroup analysis (with binary covariates) in clinical trials.}
\usage{
     DSBayes(obj, thetahat, C, lvector, control=list(), ...)
}
\arguments{
    \item{obj}{
      \code{obj} is the object from a regression model, for example, 
      linear regression or Cox proportional-hazards regression.
      If \code{obj} is specified, then \code{thetahat} and \code{C} 
      should be set to NULL. 
      }

    \item{thetahat}{
      \code{thetahat} is a vector of regression coefficients without the intercept.
      If \code{thetahat} is specified, then \code{C} should be provided 
      as well, and \code{obj} should be set at NULL. 
      } 

    \item{C}{
      \code{C} is a variance covariance matrix of regression.
      If \code{C} is given, then \code{thetahat} should also be provided,
      and \code{obj} should be set at NULL. 
      } 

    \item{lvector}{
      A vector or a matrix that denotes linear combination of the 
      parameters for which posterior estimates are desired.
      Note that, the order of the  \code{lvector} should be as follows:
      the first parameter should always be the treatment indicator, 
      then a set of binary covariates, and then the interactions between the
      treatment with covariates. See *Examples*.
      } 

    \item{control}{A list of control parameters. See *Details*.}

    \item{...}{Additional arguments.}
}


\value{
\describe{
\item{tol}{A relative accuracy for numerical quadrature. Default is \code{tol = 1.e-03}.}

\item{epsilon}{A small positive quantity to ensure proper posterior resulting from
  Jeffreys' prior. Default is \code{epsilon = 0.005}.}

\item{ci}{Level of the credible interval. Default is \code{ci = 0.95}.}

\item{k}{
  A constant value to determine the interval width for searching the Bayesian credible interval,
  from lower to upper for a maximum of the density function. Default value for 
  \code{k} is, k = qnorm((6+ci)/7) = 2.45.
  }

\item{transform}{
  If \code{transform = NULL}, then no transformation is performed.  
  If \code{transform = "logit"}, which is at default, then logit transformation 
  is applied for posterior density function to find the credibile interval,
  \code{logit(x) = log(x/(1-x))}. }                            
} 
} 

\author{Ravi Varadhan <rvaradhan@jhmi.edu> and Wenliang Yao (maintainer) <yaow080@gmail.com>}

 
\examples{

# ex1 - use given thetahat and C matrix, and set "obj=NULL".
# an example from the clinical trial reported by Fisher(1988)

thetahat  	<- c(-1.57,-0.52,-0.39,.68, 1.09, 0.68, 0.91)
names(thetahat) <- c("trt","Sex","Age","Stage","trt*sex","trt*age","trt*stage")
p <- length(thetahat)

C <- matrix(NA, p, p)
C[upper.tri(C, diag=TRUE)] <- c( .1502, .0141, .0505, .0198, .0042, .0506,
 .0389, -.0038, .0041, .0538, -.0361, -.0505, -.0042, .0039, .1037, -.0445,
 -.0042, -.0507, -.0041, -.0046, .1066,-.1209, .0037, -.0041, -.0536, -.0025,
 .0120, .1474)
C[lower.tri(C)] <- t(C)[lower.tri(t(C))]

# define lvector
trt 	<- rep(1,8)
cov 	<- as.data.frame(matrix(rep(0,24), ncol=3))
lmatrix<-as.matrix(cbind(trt,cov,rep(1:0,each=4),rep(rep(0:1,each=2),2), rep(0:1,4)))
dimnames(lmatrix)[[2]]<-c("trt","Sex","Age","Stage","trt*sex","trt*age","trt*stage")

lvector <- lmatrix[2,]    # for 1 subset  
#> lvector
#      trt       Sex       Age     Stage   trt*sex   trt*age trt*stage 
#        1         0         0         0         1         0         1 
# treatment effect for the subset of Female under 65 at stage C.
# in this case the reference group is Male, under 65 years, at stage B.

#lvector <- lmatrix       # for all 8 subsets

result <- DSBayes(NULL, thetahat, C, lvector)


################################################################################
# ex2 - use "obj" option, and set "thetahat=NULL" and "C=NULL" 
# To run ex2, you need to remove hashmark(#).

#data(simsolvd)
#simsolvd$event <- 1-simsolvd$censor
#obj <- glm(event~trt*(age+beat+lvef+cardratio+sodium),
#                 family = "binomial", data = simsolvd)
#
#para    <- as.data.frame(matrix(rep(rep(0,5),5), ncol=5))
#lmatrix <- as.matrix(cbind(rep(1,5),para[1:5,],diag(1,5)))
#dimnames(lmatrix)[[2]] <- c("trt","age","beat","lvef","cardratio","sodium",
#"trt*age","trt*beat","trt*lvef","trt*cardratio","trt*sodium")
  
#lvector   <- lmatrix[2,] 	    # for 1 subset 
#out <- DSBayes(obj, NULL, NULL, lvector)

}
\references{

   Dixon D. and Simon R. (1991).  Bayesian Subset Analysis. \emph{Biometrics}, 47, 871-881

}