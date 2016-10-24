\name{POST-package}
\alias{POST-package}
 \Rdversion{3.3.1}
 \docType{package}
 \title{
   Projection onto Orthogonal Space Testing for High Dimensional Data
 }
 \description{
   a tool to indentify gene sets that are associated with a endpoint of interest 
   by projecting the original data of a set to an orthogonal subspace.
 }
 \details{
   \tabular{ll}{
     Package: \tab POST\cr
     Type: \tab Package\cr
     Version: \tab 0.99.1\cr
     Date: \tab 2016-10-24\cr
     License: \tab GPL (>=2)\cr
     LazyLoad: \tab yes\cr
   }
   The POST (Projection onto Orthogonal Space Testing) is performed by calling function \emph{POSTglm} for linear model or 
   \emph{POSTcoxph} for Cox proportional hazard model. The genomic data such as gene expression is passed through ExpressionSet; 
   The gene set definition is defined by a gene set collection. 
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}

      Maintainer: Xueyuan Cao \email{xueyuan.cao@stjude.org}
 }
 \references{
  Hotelling H. (1936). Relations between two sets of variables. Biometrika, 28, 321-327
  
 }
 \keyword{ package }
 \examples{
 ## load data
 data(sampExprSet)
 data(sampGeneSet)
 test<-POSTglm(exprSet=sampExprSet,                         
               geneSet=sampGeneSet,                         
               lamda=0.95,                      
               seed=13,                        
               nboots=100,                      
               modd='Group ~ ',   
               family=binomial(link = "logit")) 
 }
