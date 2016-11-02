POSTglm<-function(exprSet,                         
                  geneSet,                         
                  lamda=0.95,                      
                  seed=13,                         
                  nboots=100,                      
                  model='Group ~ ',                 
                  family=binomial(link = "logit"), 
                  ...)                                                    
{
  sttime<-Sys.time()
  str.opt<-options()$stringsAsFactors
  if (!str.opt)
  {
    options(stringsAsFactors=TRUE)
  }  
  expr<-exprs(exprSet)
  clin<-pData(exprSet)
  
  spec.str<-c("FACTOR","AS.FACTOR","STRATA")
  vars<-unlist(strsplit(model, "[,() ~]+"))
  vars<-vars[!is.element(toupper(vars), spec.str)] 
  if (sum(!is.element(vars, dimnames(clin)[[2]]))!=0)
     {print("Model Specification Error:"); break()}
      
  var2<-unlist(strsplit(model, "~"))
  var2<-gsub(" ", "", var2, fixed=TRUE)
  if (var2[2]!="" ) newmodel<-paste(var2[1], "~ prb + ", var2[2], sep="")
  if (var2[2]=="" &length(vars)>1) newmodel<-paste(var2[1], "~ prb", sep="")
  if (var2[2]=="" &length(vars)==1) newmodel<-"ph ~ prb"

  nsub<-nrow(clin)
  res<-NULL
   
  set.seed(seed)
  perms<-matrix(NA, nboots, nsub)
  for (i in seq_len(nboots)) perms[i,]<-sample(nsub,replace=TRUE)
   
  nSet<-length(geneSet)  
  for (i in seq_len(nSet))
  {
    thisSet<-geneSet[[i]]
    this<-geneIds(thisSet)
    this.pot<-expr[is.element(dimnames(expr)[[1]], this), ]
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
              if (length(vars)==1)  tdd<-cbind.data.frame(prb, ph= mdat)
                else tdd<-cbind.data.frame(prb, mdat)             
              tres<-summary(glm(as.formula(newmodel), family=family, data=tdd,...))
              return(tres$coefficients[2, 3])
              }, clin[, vars])
         qstat<-t(zstat)%*%AA%*%zstat
         
         # bootstrap to estimate the correlation structure, under null.     
         boots<-matrix(NA, nboots, nf)
         for (j in seq_len(nboots))
         {
            this.boot<-perms[j,]
            boot.res<-apply(prj.pot, 1, function(prb, mdat){
              if (length(vars)==1)  tdd<-cbind.data.frame(prb, ph= mdat)
                  else tdd<-cbind.data.frame(prb, mdat)
              tres<-summary(glm(as.formula(newmodel), family=family, data=tdd,...))
              return(tres$coefficients[2, 3])
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
       }
      if (nf ==1)
      {
        tdd<-cbind.data.frame(prb=t(this.pot)%*%vec[, 1], clin[, vars])
        tres<-summary(glm(as.formula(newmodel), family=family, data=tdd,...))
        zstat<-tres$coefficients[2, 3]^2
        pval<-pchisq(zstat, 1, lower.tail=FALSE)
        this.res<-cbind(Gene=setName(thisSet),Nprobe=nrow(this.pot), 
                        Ndim=1,  Stat=zstat, Pvalue=pval)
        res<-rbind(res, this.res) 
      }    
    }    
    if (nrow(this.pot)==1)
    {
        tdd<-cbind.data.frame(prb=as.numeric(this.pot), clin[, vars])
        tres<-summary(glm(as.formula(newmodel), family=family, data=tdd,...))
        zstat<-tres$coefficients[2, 3]^2
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
