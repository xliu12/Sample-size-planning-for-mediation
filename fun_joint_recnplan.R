

joint.recnplan.assurance=function(Npilot,Nplan.start,nstep,
                        power.desired=0.8,
                        assurance.desired,
                        #meanpower.desired=NULL,
                        #medianpower.desired=NULL,
                        a_true,b_true,cp_true,alpha,
                        nrawdata,nemp){
  source("fun_joint_1pilot.R")
  
  powerdist.start=joint.powerdist(nrawdata,nemp,
                  a_true,b_true,cp_true,alpha,
                  Npilot,Nplan=Nplan.start)
  if(mean(powerdist.start>=power.desired)>assurance.desired){
    
    Nplan=Nplan.start-nstep
    repeat{
      powerdist=joint.powerdist(nrawdata,nemp,
                                      a_true,b_true,cp_true,alpha,
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
      powerdist=joint.powerdist(nrawdata,nemp,
                                      a_true,b_true,cp_true,alpha,
                                      Npilot,Nplan)
      print(c( Nplan, mean(powerdist>=power.desired) ))
      if(mean(powerdist>=power.desired)>=assurance.desired) break
      Nplan=Nplan+nstep
    }
    recnplan=Nplan
    
  }
  
  return(recnplan)
  
}

##############################################################################

joint.recnplan.mean=function(Npilot,Nplan.start,nstep,
                                  power.desired=0.8,
                                  #assurance.desired,
                                  meanpower.desired=0.8,
                                  #medianpower.desired=NULL,
                                  a_true,b_true,cp_true,alpha,
                                  nrawdata,nemp){
  source("fun_joint_1pilot.R")
  
  powerdist.start=joint.powerdist(nrawdata,nemp,
                                  a_true,b_true,cp_true,alpha,
                                  Npilot,Nplan=Nplan.start)
  if(mean(powerdist.start)>meanpower.desired){
    
    Nplan=Nplan.start-nstep
    repeat{
      powerdist=joint.powerdist(nrawdata,nemp,
                                      a_true,b_true,cp_true,alpha,
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
      powerdist=joint.powerdist(nrawdata,nemp,
                                      a_true,b_true,cp_true,alpha,
                                      Npilot,Nplan)
      print(c( Nplan, mean(powerdist) ))
      if(mean(powerdist)>=meanpower.desired) break
      Nplan=Nplan+nstep
    }
    recnplan=Nplan
    
  }
  
  return(recnplan)
  
}

##############################################################################
n.sigasigb=function(a_true,b_true,cp_true,alpha,power.desired){
  z_alpha=qnorm(1-alpha/2)
  z_beta=qnorm(power.desired)
  n_siga=( (1-a_true^2)/a_true^2 )*( z_alpha+z_beta )^2
  n_sigb=( (1-b_true^2-cp_true^2-2*a_true*b_true*cp_true)/(b_true^2*(1-a_true^2)) )*( z_alpha+z_beta )^2
  n_sigasigb=max(n_siga,n_sigb)
  return(n_sigasigb)
}






























