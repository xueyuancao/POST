POSTcoxph<-function(exprSet,                       
                  geneSet,                         
                  lamda=0.95,                      
                  nboots=100,                      
                  modd='Surv(EFSTIME, EFSCENSOR) ~ strata(arm2)',  
                  seed=13)                                                                                         
{
  sttime<-Sys.time()
  str.opt<-options()$stringsAsFactors
  if (!str.opt)
  {
    options(stringsAsFactors=TRUE)
  }
  
  #extract data
  expr<-exprs(exprSet)
  clin<-pData(exprSet)
  
  vars<-modd
  splc<-c("~", ",",  "+", "(", ")")
  for (i in 1:length(splc))
  {
    vars<-unlist(strsplit(vars, splc[i], fixed =TRUE)) 
  }
  vars<-gsub(" ", "", vars, fixed=TRUE)
  vars<-vars[!is.element(toupper(vars), c("","SURV", "FACTOR","AS.FACTOR","STRATA"))]
  
  if (sum(!is.element(vars, dimnames(clin)[[2]]))!=0)
  { print("Model Specification Error:"); break()} 
  var2<-unlist(strsplit(modd, "~"))
  var2<-gsub(" ", "", var2, fixed=TRUE)
  var2<-var2[!is.element(var2, '')]
  if ( length(var2)==2 ) newmod<-paste(var2[1], "~ prb + ", var2[2], sep="")
  if ( length(var2)==1 ) newmod<-paste(var2[1], "~ prb", sep="")
  nsub<-nrow(clin)
  res<-NULL
  
  set.seed(seed)
  perms<-matrix(NA, nboots, nsub)
  for (i in 1:nboots) perms[i,]<-sample(nsub,replace=TRUE)
    
  for (i in 1:length(geneSet)){
    thisSet<-geneSet[[i]]
    this<-geneIds(thisSet)
    this.pot<-expr[is.element(dimnames(expr)[[1]], this), ]
    if (nrow(this.pot)>1)
    {
       egd<-eigen(cov(t(this.pot)))
       val<-egd$values
       vec<-egd$vectors
       if (is.null(lamda)) nf<-length(val)
          else if (is.numeric(lamda)) 
                   nf<-min(c(1:length(val))[cumsum(val)/sum(val) >=lamda])
       if (nf >1)
       {          
         AA<-diag(val[1:nf])
         prj.pot<-t(t(this.pot)%*%vec[, 1:nf])
         zstat<-apply(prj.pot, 1, function(prb, mdat){
                tdd<-cbind.data.frame(prb, mdat)
                tt<-coxph(as.formula(newmod), data=tdd)
                return(tt$coefficients[1]/sqrt(diag(tt$var)[1]))
                }, clin[, vars])
         qstat<-t(zstat)%*%AA%*%zstat
             
         boots<-NULL
         for (j in 1:nboots)
         {
            this.boot<-perms[j,]
            boot.t<-apply(prj.pot, 1, function(prb, mdat){
                tdd<-cbind.data.frame(prb, mdat)
                tt<-coxph(as.formula(newmod), data=tdd)
                return(tt$coefficients[1]/sqrt(diag(tt$var)[1]))
                }, clin[this.boot, vars])
          boots<-rbind(boots, boot.t)
         }
         sigma<-cor(boots)
         Asigma<-AA%*%sigma
         npdlam<-nearPD(Asigma)$eigenvalues
         lambda<-npdlam[npdlam>0]          
         gchiappr<-davies(qstat, lambda)$Qq
         this.res<-cbind(Gene=setName(thisSet), Nprobe=nrow(this.pot),  
                         Ndim=nf, Stat=qstat, Pvalue=gchiappr)
         res<-rbind(res, this.res)
         #print(paste(i, nf, sep=": "))
       }
      if (nf ==1)
      {
        tdd<-cbind.data.frame(prb=t(this.pot)%*%vec[, 1], clin[, vars])
        tt<-coxph(as.formula(newmod), data=tdd)
        zstat<-tt$coefficients[1]/sqrt(diag(tt$var)[1])
        pval<-pchisq(zstat, 1, lower.tail=FALSE)
        this.res<-cbind(Gene=setName(thisSet),Nprobe=nrow(this.pot), 
                        Ndim=1,  Stat=zstat, Pvalue=pval)
        res<-rbind(res, this.res) 
        #print(paste(i, nf, sep=": "))
      }    
    }    
    if ( nrow(this.pot)==1)
    {
        tdd<-cbind.data.frame(prb=as.numeric(this.pot), clin[, vars])
        tt<-coxph(as.formula(newmod), data=tdd)
        zstat<-tt$coefficients[1]/sqrt(diag(tt$var)[1])
        pval<-pchisq(zstat, 1, lower.tail=FALSE)
        this.res<-cbind(Gene=setName(thisSet),Nprobe=nrow(this.pot), 
                        Ndim=1,  Stat=zstat, Pvalue=pval)
        res<-rbind(res, this.res)
        #print(paste(i, 1, sep=": "))
    }
  } 
  dimnames(res)[[2]]<-c("GeneSet", "Nprobe", "Nproj", 'Stat', 'p.value')
  options(stringsAsFactors=str.opt)
  endtime<-Sys.time()
  print(paste("POST taking: ", 
              format(round(difftime(endtime, sttime, units ="mins"),2)), sep="")) 
  return(res)
}
