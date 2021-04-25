# reccommended sample size
rm(list = ls())
source("fun_joint_recnplan.R"); source("fun_perboot_recnplan.R")

recnplan <- function(
  std.ahat=.39, std.bhat=.39, std.cphat=0,
  method="joint",
  goal="meanpower",
  Npilot=250,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000){
  
  a_true=std.ahat; b_true=std.bhat; cp_true=std.cphat
  if(power.desired!=meanpower.desired) power.desired=meanpower.desired
  if(is.null(Nplan.start)) Nplan.start=ceiling(n.sigasigb(a_true,b_true,cp_true,alpha,power.desired))
  if(method=="joint" & goal=="assurance"){
    rec_nplan=joint.recnplan.assurance(
      Npilot=Npilot,
      Nplan.start=Nplan.start,
      nstep=nstep,
      power.desired=power.desired,
      assurance.desired=assurance.desired,
      a_true=std.ahat, b_true=std.bhat, cp_true=std.cphat,
      alpha=alpha,
      nrawdata=nrawdata, nemp=nemp)
  }
  if(method=="joint" & goal=="meanpower"){
    rec_nplan=joint.recnplan.mean(
      Npilot=Npilot,
      Nplan.start=Nplan.start,
      nstep=nstep,
      power.desired=power.desired,
      meanpower.desired=meanpower.desired,
      a_true=std.ahat, b_true=std.bhat, cp_true=std.cphat,
      alpha=alpha,
      nrawdata=nrawdata, nemp=nemp)
  }
  if(method=="perboot" & goal=="assurance"){
    rec_nplan=perboot.recnplan.assurance(
      Npilot=Npilot,
      Nplan.start=Nplan.start,
      nstep=nstep,
      power.desired=power.desired,
      assurance.desired=assurance.desired,
      a_true=std.ahat, b_true=std.bhat, cp_true=std.cphat,
      alpha=alpha,
      nrawdata=nrawdata, nemp=nemp, nboot=nboot)
  }
  if(method=="perboot" & goal=="meanpower"){
    rec_nplan=perboot.recnplan.mean(
      Npilot=Npilot,
      Nplan.start=Nplan.start,
      nstep=nstep,
      power.desired=power.desired,
      meanpower.desired=meanpower.desired,
      a_true=std.ahat, b_true=std.bhat, cp_true=std.cphat,
      alpha=alpha,
      nrawdata=nrawdata, nemp=nemp, nboot=nboot)
  }
  
  return(rec_nplan)
  
}


# empirical example 1: n=73 -------
std.ahat1=0.250
std.bhat1=0.258
std.cphat1=0.161
n1=Npilot1=73

joint_recnplan1_80assurance = recnplan(
  std.ahat=std.ahat1, std.bhat=std.bhat1, std.cphat=std.cphat1,
  method="joint",
  goal="assurance",
  Npilot=Npilot1,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

joint_recnplan1_80assurance

joint_recnplan1_80mean = recnplan(
  std.ahat=std.ahat1, std.bhat=std.bhat1, std.cphat=std.cphat1,
  method="joint",
  goal="meanpower",
  Npilot=Npilot1,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

joint_recnplan1_80mean

perboot_recnplan1_80assurance = recnplan(
  std.ahat=std.ahat1, std.bhat=std.bhat1, std.cphat=std.cphat1,
  method="perboot",
  goal="assurance",
  Npilot=Npilot1,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

perboot_recnplan1_80assurance

perboot_recnplan1_80mean = recnplan(
  std.ahat=std.ahat1, std.bhat=std.bhat1, std.cphat=std.cphat1,
  method="perboot",
  goal="meanpower",
  Npilot=Npilot1,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

perboot_recnplan1_80mean

# empirical example 2: n=332 -------
n2=Npilot2=332
std.ahat2=-0.1547
# se^2=MSE/(n*sx^2)
sx2=sqrt(5.2555/(0.2679^2)/332)
# sy^2=MSE/(1-R2)
sm2=sqrt(5.2555/(1-0.0239))
sy2=sqrt(1.6470/(1-0.4189))
std.bhat2=0.4638*sm2/sy2
std.cphat2=-0.1357*sx2/sy2

joint_recnplan2_80assurance = recnplan(
  std.ahat=std.ahat2, std.bhat=std.bhat2, std.cphat=std.cphat2,
  method="joint",
  goal="assurance",
  Npilot=Npilot2,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

joint_recnplan2_80assurance

joint_recnplan2_80mean = recnplan(
  std.ahat=std.ahat2, std.bhat=std.bhat2, std.cphat=std.cphat2,
  method="joint",
  goal="meanpower",
  Npilot=Npilot2,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

joint_recnplan2_80mean

perboot_recnplan2_80assurance = recnplan(
  std.ahat=std.ahat2, std.bhat=std.bhat2, std.cphat=std.cphat2,
  method="perboot",
  goal="assurance",
  Npilot=Npilot2,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

perboot_recnplan2_80assurance

perboot_recnplan2_80mean = recnplan(
  std.ahat=std.ahat2, std.bhat=std.bhat2, std.cphat=std.cphat2,
  method="perboot",
  goal="meanpower",
  Npilot=Npilot2,
  Nplan.start=NULL,
  nstep=1,
  power.desired=0.8,
  assurance.desired=.8,
  meanpower.desired=.8,
  alpha=.05,
  nrawdata=1000, nboot=2000, nemp=1000
)

perboot_recnplan2_80mean


save.image("emp_recnplan.RData")





