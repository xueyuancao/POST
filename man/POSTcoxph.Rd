\name{POSTcoxph}
\alias{POSTcoxph}
 \title{POST Analysis of Cox Proportional Hazard Model}
 \description{Peform POST analysis with Cox proportional hazard model}
\usage{
POSTcoxph (exprSet, geneSet, lamda = 0.95, nboots = 100, 
           model = "Surv(EFSTIME, EFSCENSOR) ~ strata(arm2)", 
           seed = 13, ...) 
}
 \arguments{
 \item{exprSet}{an ExpressionSet class contains minimum of exprs (expression matrix) of genomic 
                data such as gene expression and phenoData (AnnotatedDataFrame of end point data). 
                Please refer to Biobase for details on how to create such an ExpressionSet expression set. }
 \item{geneSet}{a GeneSetCollection class with minimum of setName and geneIDs for each GeneSet. 
                Please refer to \emph{GSEABase} for how to create such a GeneSetCollection class.}
 \item{lamda}{a predefined fraction of 0 to 1: the fraction of variation retained in the selected orthogonal projections} 
 \item{nboots}{number of bootstrap samples, default=100}
 \item{model}{a character string to define model}
 \item{seed}{seed of random number generator}
 \item{...}{control arguments to be used in Cox proportinal hazard model,
            default control arguments is used if they are not supplied directly}                                                 
 }
 \details{The function performs POST analysis for association of gene set with time to event endpoint 
          in Cox proportional hazard model framework.}
 \value{
       a data frame with 5 columns
       \item{GeneSet}{Gene set id}
       \item{Nprobe}{numbe of probes in the ExpressionSet annotated to the gene set}
       \item{Nproj}{number of projected vectors included in the POST analysis}
       \item{Stat}{statistic of quadratic form}  
       \item{p.value}{p value of generalized Chi-square approximation} 
 }
 \author{Xueyuan Cao \email{Xueyuan.cao@stjude.org}, Stanley Pounds \email{stanley.pounds@stjude.org}}
 \examples{
 ## load data
  data(sampExprSet)
  data(sampGeneSet)
 test<-POSTcoxph(exprSet=sampExprSet,                         
               geneSet=sampGeneSet,                         
               lamda=0.95,                      
               seed=13,                        
               nboots=100,                      
               model='Surv(time, censor) ~ ') 
 }
