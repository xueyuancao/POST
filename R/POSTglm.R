POSTglm<-function(exprSet,                         
                  geneSet,                         
                  lamda=0.95,                      
                  seed=13,                         
                  nboots=100,                      
                  modd='Group ~ ',                 
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
  
  vars<-modd
  splc<-c("~", "+", "(", ")")
  for (i in 1:length(splc))
  {
    vars<-unlist(strsplit(vars, splc[i], fixed =TRUE)) 
  }
  vars<-gsub(" ", "", vars, fixed=TRUE)
  vars<-vars[!is.element(toupper(vars), c("", "FACTOR","AS.FACTOR","STRATA"))]
  if (sum(!is.element(vars, dimnames(clin)[[2]]))!=0)
     {print("Model Specification Error:"); break()} 
  var2<-unlist(strsplit(modd, "~"))
  var2<-gsub(" ", "", var2, fixed=TRUE)
  var2<-var2[!is.element(var2, '')]
  
  #Define new model form to include an extra term for projected data
  if (length(var2)==2 ) newmod<-paste(var2[1], " ~ prb + ", var2[2], sep="")
  if (length(var2)==1 &length(vars)>1) newmod<-paste(var2[1], " ~ prb", sep="")
  if (length(var2)==1 &length(vars)==1) newmod<-"ph ~ prb"
  nsub<-nrow(clin)
  res<-NULL
   
  set.seed(seed)
  perms<-matrix(NA, nboots, nsub)
  for (i in 1:nboots) perms[i,]<-sample(nsub,replace=TRUE)
  
  for (i in 1:length(geneSet))
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
          else if (is.numeric(lamda)) 
               nf<-min(c(1:length(val))[cumsum(val)/sum(val) >=lamda])
       if (nf >1)
       {          
         AA<-diag(val[1:nf])
         prj.pot<-t(t(this.pot)%*%vec[, 1:nf])
         zstat<-apply(prj.pot, 1, function(prb, mdat){
              if (length(vars)==1)  tdd<-cbind.data.frame(prb, ph= mdat)
                else tdd<-cbind.data.frame(prb, mdat)             
              #tres<-summary(glm(as.formula(newmod), family=family, data=tdd,...))
              tres<-summary(glm(as.formula(newmod), family=family, data=tdd))
              return(tres$coefficients[2, 3])
              }, clin[, vars])
         qstat<-t(zstat)%*%AA%*%zstat
         
         # bootstrap to estimate the correlation structure, under null.     
         boots<-NULL
         for (j in 1:nboots)
         {
            this.boot<-perms[j,]
            boot.t<-apply(prj.pot, 1, function(prb, mdat){
              if (length(vars)==1)  tdd<-cbind.data.frame(prb, ph= mdat)
                  else tdd<-cbind.data.frame(prb, mdat)
              tres<-summary(glm(as.formula(newmod), family=family, data=tdd,...))
              return(tres$coefficients[2, 3])
              }, clin[this.boot, vars])
          boots<-rbind(boots, boot.t)
         }
         sigma<-cor(boots)
         Asigma<-AA%*%sigma
         npdlam<-nearPD(Asigma)$eigenvalues
         lambda<-npdlam[npdlam>0] 
         #Asigma.decomp<-eigen(Asigma)
         #lambda<-Asigma.decomp$values[Asigma.decomp$values>0]           
         gchiappr<-davies(qstat, lambda)$Qq
         this.res<-cbind(Gene=setName(thisSet), Nprobe=nrow(this.pot),  
                         Ndim=nf, Stat=qstat, Pvalue=gchiappr)
         res<-rbind(res, this.res)
         #print(paste(i, nf, sep=": "))
       }
      if (nf ==1)
      {
        tdd<-cbind.data.frame(prb=t(this.pot)%*%vec[, 1], clin[, vars])
        tres<-summary(glm(as.formula(newmod), family=family, data=tdd,...))
        zstat<-tres$coefficients[2, 3]^2
        pval<-pchisq(zstat, 1, lower.tail=FALSE)
        this.res<-cbind(Gene=setName(thisSet),Nprobe=nrow(this.pot), 
                        Ndim=1,  Stat=zstat, Pvalue=pval)
        res<-rbind(res, this.res) 
        #print(paste(i, nf, sep=": "))
      }    
    }    
    if (nrow(this.pot)==1)
    {
        tdd<-cbind.data.frame(prb=as.numeric(this.pot), clin[, vars])
        tres<-summary(glm(as.formula(newmod), family=family, data=tdd,...))
        zstat<-tres$coefficients[2, 3]^2
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
