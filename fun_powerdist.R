
source('fun_powerdist.joint.R')
source('fun_powerdist.perboot.R')

powerdist<- function(
  std.ahat=.39, std.bhat=.39, std.cphat=0,
  method="joint", # perboot uses the percentile bootstrap test; it will take longer than the joint significance test 
  Npilot=250,
  Nplan=68,
  alpha=0.05, 
  K = 1000, # for uncertainty in effect size estimates -- K draws from the distribution of a b cp   
  nemp=1000, # for empirical power
  nboot=1000 # for percentile bootstrap
){
  
  if (method=='joint') {
    powerdist_vec = joint.powerdist(
      Npilot=Npilot,
      Nplan=Nplan,
      a_true=std.ahat, b_true=std.bhat, cp_true=std.cphat,
      alpha=alpha,
      nrawdata=K, nemp=nemp
    )
  }
  
  if (method=='perboot') {
    powerdist_vec = perboot.powerdist(
      Npilot=Npilot,
      Nplan=Nplan,
      a_true=std.ahat, b_true=std.bhat, cp_true=std.cphat,
      alpha=alpha,
      nrawdata=K, nemp=nemp, nboot = nboot
    )
  }
  
  return( powerdist_vec )
}



