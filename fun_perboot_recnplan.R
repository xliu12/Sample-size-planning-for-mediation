library(parallel)
mccores=2
# perboot.powerdist -------------
perboot.powerdist=function(nrawdata,nemp,
                           a_true,b_true,cp_true,alpha,nboot,
                           Npilot,Nplan){
  
  perboot.power.1pilot=function(inrawdata){
    nemp=nemp; a_true=a_true; b_true=b_true; cp_true=cp_true; alpha=0.05; nboot=nboot; N=Npilot; Nplan=Nplan;
    Var.e_M=1-a_true^2
    Var.e_Y=1-b_true^2-cp_true^2-2*a_true*b_true*cp_true
    powerdist=NULL
    #for(r in 1:(1000/nrawdata) ){
    
    x=rnorm(N)
    m=a_true*x+rnorm(N,0,sqrt(Var.e_M))
    y=b_true*m+cp_true*x+rnorm(N,0,sqrt(Var.e_Y))
    
    txx=t(x)%*%x
    tmx=t(m)%*%x
    tmm=t(m)%*%m
    
    ahat=tmx/txx
    bhat=(t(c(txx)*m-c(tmx)*x)%*%y)/(txx*tmm-tmx*tmx) 
    cphat=(t(c(tmm)*x-c(tmx)*m)%*%y)/(txx*tmm-tmx*tmx) 
    Var.e_M.hat=t(m-c(ahat)*x)%*%(m-c(ahat)*x)/(N-1)
    Var.e_Y.hat=t(y-c(cphat)*x-c(bhat)*m)%*%(y-c(cphat)*x-c(bhat)*m)/(N-2)  
    
    xrep=rnorm(nemp*Nplan)
    mrep=c(ahat)*xrep+rnorm(nemp*Nplan,0,sqrt(Var.e_M.hat))
    yrep=c(cphat)*xrep+c(bhat)*mrep+rnorm(nemp*Nplan,0,sqrt(Var.e_Y.hat))
    xrep=matrix(xrep,nrow = Nplan)
    mrep=matrix(mrep,nrow = Nplan)
    yrep=matrix(yrep,nrow = Nplan)
    xrep_mrep_yrep=rbind(xrep,mrep,yrep)
    
    ind.boot=function(xmy){
      indboot=NULL
      for(b in 1:nboot){
        bootsp=sample(1:(length(xmy)/3), length(xmy)/3, replace=T)
        
        x=xmy[1:(length(xmy)/3)][bootsp]
        m=xmy[(length(xmy)/3+1):(length(xmy)*2/3)][bootsp]
        y=xmy[(length(xmy)*2/3+1):length(xmy)][bootsp]
        
        txx=t(x)%*%x
        tmx=t(m)%*%x
        tmm=t(m)%*%m
        
        ahat=tmx/txx
        bhat=(t(c(txx)*m-c(tmx)*x)%*%y)/(txx*tmm-tmx*tmx) 
        
        indboot=c(indboot,ahat*bhat)
      }
      #ifsig=prod(quantile(indboot,c(alpha/2,1-alpha/2))>0)
      ifsig=( prod(quantile(indboot,c(alpha/2,1-alpha/2)))>0 )
      return(ifsig)
    }
    
    ind_b_sig=apply(xrep_mrep_yrep,2,ind.boot)
    
    power.emp=mean(ind_b_sig)
    
    return(power.emp)
    
  }
  
  
  #powerdist.list = lapply(1:nrawdata, perboot.power.1pilot)
  powerdist.list = mclapply(1:nrawdata, perboot.power.1pilot, mc.preschedule = F, mc.cores = mccores)
  powerdist = unlist(powerdist.list)
  return(powerdist)
}

# perboot.recnplan.assurance --------------

perboot.recnplan.assurance=function(Npilot,Nplan.start,nstep,
                                  power.desired=0.8,
                                  assurance.desired,
                                  #meanpower.desired=NULL,
                                  #medianpower.desired=NULL,
                                  a_true,b_true,cp_true,alpha,
                                  nrawdata,nemp,nboot){
  
  powerdist.start=perboot.powerdist(nrawdata,nemp,
                                  a_true,b_true,cp_true,alpha,nboot,
                                  Npilot,Nplan=Nplan.start)
  if(mean(powerdist.start>=power.desired)>assurance.desired){
    
    Nplan=Nplan.start-nstep
    repeat{
      powerdist=perboot.powerdist(nrawdata,nemp,
                                a_true,b_true,cp_true,alpha,nboot,
                                Npilot,Nplan)
      print(c( Nplan, mean(powerdist>=power.desired) ))
      if(mean(powerdist>=power.desired)<=assurance.desired) break
      Nplan=Nplan-nstep
    }
    recnplan=Nplan+nstep
    
  } 
  
  else{
    
    Nplan=Nplan.start+nstep
    repeat{
      powerdist=perboot.powerdist(nrawdata,nemp,
                                a_true,b_true,cp_true,alpha,nboot,
                                Npilot,Nplan)
      print(c( Nplan, mean(powerdist>=power.desired) ))
      if(mean(powerdist>=power.desired)>=assurance.desired) break
      Nplan=Nplan+nstep
    }
    recnplan=Nplan
    
  }
  
  return(recnplan)
  
}

# perboot.recnplan.mean -----------------

perboot.recnplan.mean=function(Npilot,Nplan.start,nstep,
                             power.desired=0.8,
                             #assurance.desired,
                             meanpower.desired=0.8,
                             #medianpower.desired=NULL,
                             a_true,b_true,cp_true,alpha,
                             nrawdata,nemp,nboot){
  
  powerdist.start=perboot.powerdist(nrawdata,nemp,
                                  a_true,b_true,cp_true,alpha,nboot,
                                  Npilot,Nplan=Nplan.start)
  if(mean(powerdist.start)>meanpower.desired){
    
    Nplan=Nplan.start-nstep
    repeat{
      powerdist=perboot.powerdist(nrawdata,nemp,
                                a_true,b_true,cp_true,alpha,nboot,
                                Npilot,Nplan)
      print(c( Nplan, mean(powerdist) ))
      if(mean(powerdist)<=meanpower.desired) break
      Nplan=Nplan-nstep
    }
    recnplan=Nplan+nstep
    
  } 
  
  else{
    
    Nplan=Nplan.start+nstep
    repeat{
      powerdist=perboot.powerdist(nrawdata,nemp,
                                a_true,b_true,cp_true,alpha,nboot,
                                Npilot,Nplan)
      print(c( Nplan, mean(powerdist) ))
      if(mean(powerdist)>=meanpower.desired) break
      Nplan=Nplan+nstep
    }
    recnplan=Nplan
    
  }
  
  return(recnplan)
  
}
































