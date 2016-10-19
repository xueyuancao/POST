\name{POSTglm}
\alias{POSTglm}
 \title{POST of Generalized Linear Model}
 \description{Peform POST analysis with linear model}
\usage{
POSTglm (exprSet, geneSet, lamda = 0.95, seed = 13, nboots = 100, 
    modd = "Group ~ ", family = binomial(link = "logit"), ...) 
}
 \arguments{
 \item{exprSet}{an ExpressionSet class contains minimum of exprs (expression matrix) of genomic data 
                such as gene expression and phenoData (AnnotatedDataFrame of end point data). 
                Please refer to Biobase for details on how to create such an ExpressionSet expression set. }
 \item{geneSet}{a GeneSetCollection class with minimum of setName and geneIDs for each GeneSet. 
                Please refer to \emph{GSEABase} for how to create such a GeneSetCollection class.}
 \item{lamda}{a predefined fraction} 
 \item{seed}{seed of random number generator}
 \item{nboots}{number of bootstrap samples, default=100}
 \item{modd}{a character string to define the model}
 \item{family}{a description of the error distribution and link function to be used in the model. 
               For glm, this can be a character string naming a family function, 
               a family function or the result of a call to a family function.}
 \item{...}{For glm: arguments to be used to form the default control argument if it is not supplied directly}                         
 }
 \details{The function performs POST analysis for association of gene sets with phenotype 
          in generalized linear model framework.}
 \value{
       a data frame with 5 columns
       \item{GeneSet}{Gene set id}
       \item{Nprobe}{numbe of probes in the ExpressionSet annotated to the gene set}
       \item{Nproj}{number of projected vectors included in the POST analysis}
       \item{Stat}{statistic of quadratic form}  
       \item{p.value}{p value} 
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
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
