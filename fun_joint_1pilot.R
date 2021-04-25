
# alpha=.05
# nrawdata=1
# nemp=1000
# a_true=.14
# b_true=.14
# cp_true=0
# N=5000
# Nplan=551
# 
# 
# system.time(joint.power.1pilot(nemp,
#                           a_true,b_true,cp_true,alpha,
#                           N,Nplan))

joint.power.1pilot=function(nemp,
                              a_true,b_true,cp_true,alpha,
                              N,Nplan){
  Var.e_M=1-a_true^2
  Var.e_Y=1-b_true^2-cp_true^2-2*a_true*b_true*cp_true
  powerdist=NULL

  
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
  
  ind.joint=function(xmy){
    x=xmy[1:(length(xmy)/3)]
    m=xmy[(length(xmy)/3+1):(length(xmy)*2/3)]
    y=xmy[(length(xmy)*2/3+1):length(xmy)]
    
    txx=t(x)%*%x
    tmx=t(m)%*%x
    tmm=t(m)%*%m
    tmy=t(m)%*%y
    txy=t(x)%*%y
    
    ahat=tmx/txx
    bhat=(txx*tmy-tmx*txy)/(txx*tmm-tmx*tmx) 
    cphat=(tmm*txy-tmx*tmy)/(txx*tmm-tmx*tmx) 
    Var.e_M.hat=t(m-c(ahat)*x)%*%(m-c(ahat)*x)/(length(x)-1)
    Var.e_Y.hat=t(y-c(cphat)*x-c(bhat)*m)%*%(y-c(cphat)*x-c(bhat)*m)/(length(x)-2)  
    
    Var.ahat.hat=Var.e_M.hat/txx
    Var.bhat.hat=Var.e_Y.hat*txx/(txx*tmm-tmx*tmx) 
    
    ifsig=(ahat^2/Var.ahat.hat>qf(alpha,1,length(x)-1,lower.tail = F))*
      (bhat^2/Var.bhat.hat>qf(alpha,1,length(x)-2,lower.tail = F))
    
    return(ifsig)
  }
  
  ind_b_sig=apply(xrep_mrep_yrep,2,ind.joint)
  
  power.emp=mean(ind_b_sig)
  
  return(power.emp)
  
}



joint.powerdist=function(nrawdata,nemp,
                         a_true,b_true,cp_true,alpha,
                         N,Nplan){
  powerdist=NULL
  for(raw in 1:nrawdata){
    temp=joint.power.1pilot(nemp,
                            a_true,b_true,cp_true,alpha,
                            N,Nplan)
    powerdist=c(powerdist,temp)
  }
  return(powerdist)
}
 


