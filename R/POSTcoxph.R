POSTcoxph<-function(exprSet,                       
                  geneSet,                         
                  lamda=0.95,                      
                  nboots=100,                      
                  model='Surv(EFSTIME, EFSCENSOR) ~ strata(arm2)',  
                  seed=13,...)                                                                                         
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
  
  spec.str<-c("SURV", "FACTOR","AS.FACTOR","STRATA")
  vars<-unlist(strsplit(model, "[,() ~]+"))
  vars<-vars[!is.element(toupper(vars), spec.str)] 
  if (sum(!is.element(vars, dimnames(clin)[[2]]))!=0)
  { print("Model Specification Error:"); break()}
   
  var2<-unlist(strsplit(model, "~"))
  var2<-gsub(" ", "", var2, fixed=TRUE)
  if (var2[2]!="" ) newmodel<-paste(var2[1], "~ prb + ", var2[2], sep="")
  if (var2[2]=="" ) newmodel<-paste(var2[1], "~ prb", sep="")
  nsub<-nrow(clin)
  res<-NULL
  
  set.seed(seed)
  perms<-matrix(NA, nboots, nsub)
  for (i in seq_len(nboots)) perms[i,]<-sample(nsub,replace=TRUE)
  
  nSet<-length(geneSet)  
  for (i in seq_len(nSet))
  {
    thisSet<-geneSet[[i]]
    thisIds<-geneIds(thisSet)
    this.pot<-expr[is.element(dimnames(expr)[[1]], thisIds), ]
    if (nrow(this.pot)>1)
    {
       egd<-eigen(cov(t(this.pot)))
       val<-egd$values
       vec<-egd$vectors
       if (is.null(lamda)) nf<-length(val)
       if (is.numeric(lamda)) 
           nf<-min(seq_len(length(val))[cumsum(val)/sum(val) >=lamda])
       if (nf >1)
       {          
         AA<-diag(val[seq_len(nf)])
         prj.pot<-t(t(this.pot)%*%vec[, seq_len(nf)])
         zstat<-apply(prj.pot, 1, function(prb, mdat){
                tdd<-cbind.data.frame(prb, mdat)
                tt<-coxph(as.formula(newmodel), data=tdd ,...)
                return(tt$coefficients[1]/sqrt(diag(tt$var)[1]))
                }, clin[, vars])
         qstat<-t(zstat)%*%AA%*%zstat
             
         boots<-matrix(NA, nboots, nf)
         for (j in seq_len(nboots))
         {
            this.boot<-perms[j,]
            boot.res<-apply(prj.pot, 1, function(prb, mdat){
                tdd<-cbind.data.frame(prb, mdat)
                tt<-coxph(as.formula(newmodel), data=tdd,...)
                return(tt$coefficients[1]/sqrt(diag(tt$var)[1]))
                }, clin[this.boot, vars])
            boots[j, ]<-boot.res
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
        tt<-coxph(as.formula(newmodel), data=tdd ,...)
        zstat<-tt$coefficients[1]/sqrt(diag(tt$var)[1])
        pval<-pchisq(zstat, 1, lower.tail=FALSE)
        this.res<-cbind(Gene=setName(thisSet),Nprobe=nrow(this.pot), 
                        Ndim=1,  Stat=zstat, Pvalue=pval)
        res<-rbind(res, this.res) 
      }    
    }    
    if ( nrow(this.pot)==1)
    {
        tdd<-cbind.data.frame(prb=as.numeric(this.pot), clin[, vars])
        tt<-coxph(as.formula(newmodel), data=tdd,...)
        zstat<-tt$coefficients[1]/sqrt(diag(tt$var)[1])
        pval<-pchisq(zstat, 1, lower.tail=FALSE)
        this.res<-cbind(Gene=setName(thisSet),Nprobe=nrow(this.pot), 
                        Ndim=1,  Stat=zstat, Pvalue=pval)
        res<-rbind(res, this.res)
    }
  } 
  dimnames(res)[[2]]<-c("GeneSet", "Nprobe", "Nproj", 'Stat', 'p.value')
  options(stringsAsFactors=str.opt)
  endtime<-Sys.time()
  print(paste("POST taking: ", 
              format(round(difftime(endtime, sttime, units ="mins"),2)), sep="")) 
  return(res)
}
