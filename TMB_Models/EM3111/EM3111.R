
######################################
#Implementing TMB on supercomputer
######################################

require(TMB)
source("/blue/edvcamp/nfisch/Chapter_4/fit_tmb.R")

run_TMB<-function(N, model, surv, Est_Rsd, REML, PE){
 
 if(substr(model,8,8)=="1"){
  comp<-"rnd"
 }else if (substr(model,8,8) %in% c("2","3")){
  comp<-"corr"
 }
 
  if(substr(model,1,4)=="OM11"){
    OM<-paste0("RF_Dome_",PE,"_40yr_")
    which_dat<-paste0("dat_file_Indexnom_comp",comp,".dat")
  }else if(substr(model,1,4)=="OM12"){
    OM<-paste0("RF_Dome_",PE,"_80yr_")
    which_dat<-paste0("dat_file_Indexnom_comp",comp,".dat")
  }else if(substr(model,1,4)=="OM21"){
    OM<-paste0("RF_Dome_",PE,"_40yr_")
    which_dat<-paste0("dat_file_Indexalt_comp",comp,".dat")
  }else if(substr(model,1,4)=="OM22"){
    OM<-paste0("RF_Dome_",PE,"_80yr_")
    which_dat<-paste0("dat_file_Indexalt_comp",comp,".dat")
  }else if(substr(model,1,4)=="OM31"){
    OM<-paste0("GM_",PE,"_40yr_")
    which_dat<-paste0("dat_file_Indexalt_comp",comp,".dat")
  }else if(substr(model,1,4)=="OM32"){
    OM<-paste0("GM_",PE,"_80yr_")
    which_dat<-paste0("dat_file_Indexalt_comp",comp,".dat")
  }else if(substr(model,1,4)=="OM41"){
    OM<-paste0("GM_",PE,"_40yr_")
    which_dat<-paste0("dat_file_Indexnom_comp",comp,".dat")
  }else if(substr(model,1,4)=="OM42"){
    OM<-paste0("GM_",PE,"_80yr_")
    which_dat<-paste0("dat_file_Indexnom_comp",comp,".dat")
  }
  
setwd(paste0("/blue/edvcamp/nfisch/Chapter_4/",substr(model,9,14)))
#Compile and load model 
if(surv==TRUE){
   TMB_name<-paste0(substr(model,9,14),"_wsurv")
} else if (surv==FALSE){
   TMB_name<-paste0(substr(model,9,14),"_nosurv")
  }
  
 compile(paste0(TMB_name,".cpp"))
 dyn.load(dynlib(TMB_name))
 
 for (i in N){
 
#Reading dat file for specific OM replicate
  if (substr(model,8,8)=="1"){
   x<-scan(paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",OM,substr(model,5,7),".",substr(model,8,8),"_",i,"/",which_dat))
  }else if (substr(model,8,8) %in% c("2","3")){
   x<-scan(paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",OM,substr(model,5,7),".",as.numeric(substr(model,8,8))-1,"_",i,"/",which_dat))
  }
  
  #Extracting arguments to make a dat file for TMB
  dat<-list(fyear=x[1], lyear=x[2], fage=x[3], lage=x[4], years=c(as.integer(x[1]):as.integer(x[2])), ages=c(as.integer(x[3]):as.integer(x[4])),
            M_vec=x[5:25],
            obs_harv=x[26:(25+as.integer(x[2]))],
            obs_fishery_cpue=x[(25+as.integer(x[2])+1):(25+as.integer(x[2])*2)],
            obs_fishery_comp=matrix(x[(25+as.integer(x[2])*2+1):(25+as.integer(x[2])*2+(x[2]*(x[4]+1)))],nrow=as.integer(x[2]), ncol=as.integer(x[4]+1), byrow=T),
            SS_fishery=x[(25+as.integer(x[2])*2+(x[2]*(x[4]+1))+1):((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2]))],
            obs_FIM_CPUE=x[((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])+1):((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2)],
            obs_FIM_comp=matrix(x[((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+1):((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+(x[2]*(x[4]+1)))],nrow=as.integer(x[2]), ncol=as.integer(x[4]+1), byrow=T),
            SS_FIM=x[((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+(x[2]*(x[4]+1))+1):((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+(x[2]*(x[4]+1))+as.integer(x[2]))],
            Fecund_aa=x[((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+(x[2]*(x[4]+1))+as.integer(x[2])+1):((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+(x[2]*(x[4]+1))+as.integer(x[2])+x[4]+1)],
            Wtage=x[((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+(x[2]*(x[4]+1))+as.integer(x[2])+x[4]+2):((25+as.integer(x[2])*2+(x[2]*(x[4]+1)))+as.integer(x[2])*2+(x[2]*(x[4]+1))+as.integer(x[2])+(x[4]+1)*2)],
            test=x[(length(x)-2):length(x)])

  if (substr(model,4,4)=="1"){
  lmax_F<-c(-3.397353, -3.274413, -3.133463, -2.975869, -2.804664, -2.624553, -2.441559, -2.262295, -2.093026, -1.938777, -1.802763, -1.686263, -1.588892,
   -1.509104,-1.444715, -1.393338, -1.352676, -1.320676, -1.295588, -1.275967, -1.260648, -1.248699, -1.239385, -1.232128, -1.226475, -1.222071, -1.218641,
   -1.215970,-1.213890, -1.212270, -1.234272, -1.257415, -1.281760, -1.307367, -1.334304, -1.362638, -1.392442, -1.423792, -1.456769, -1.491458)
  }else if (substr(model,4,4)=="2"){
  lmax_F<-c(-3.452034, -3.397353, -3.338161, -3.274413, -3.206141, -3.133463, -3.056598, -2.975869, -2.891711, -2.804664, -2.715369, -2.624553, -2.533006,
   -2.441559, -2.351049, -2.262295, -2.176059, -2.093026, -2.013777, -1.938777,-1.868367, -1.802763, -1.742064, -1.686263, -1.635262, -1.588892, -1.546928,
   -1.509104, -1.475133, -1.444715,-1.417548, -1.393338, -1.371803, -1.352676, -1.335710, -1.320676, -1.307365, -1.295588, -1.285173, -1.275967,-1.267833,
   -1.260648, -1.254303, -1.248699, -1.243753, -1.239385, -1.235531, -1.232128, -1.229125, -1.226475,-1.224135, -1.222071, -1.220249, -1.218641, -1.217222,
   -1.215970, -1.214865, -1.213890, -1.213029, -1.212270,-1.223131, -1.234272, -1.245697, -1.257415, -1.269433, -1.281760, -1.294402, -1.307367, -1.320665,
   -1.334304, -1.348291, -1.362638, -1.377351, -1.392442, -1.407919, -1.423792, -1.440072, -1.456769, -1.473894, -1.491458)
  }

  #Parameters
  par <- list(log_M_scalar=log(0.099),
              log_q=-12.65,
              log_q_FIM=-14.275,
              log_recruit_devs=rep(0,dat$lyear+dat$lage),
              steepness=0.99,
              log_R0_FLA=17.3204,
              log_sigma_rec=log(0.3),
              log_cv_fishery=log(0.05),
              log_cv_fishery_CPUE=log(0.25),
              log_cv_FIM_CPUE=-1.84,
              FIM_sellogis_k=2,
              FIM_sellogis_midpt=2,
              ltheta=1,
              ltheta_FIM=3,
              B1=3.170388,
              B2=-4.187365,
              B3=1.048499,
              B4=1.028458,
              B5=-4.906849,
              B6=1.116807,
              log_fint=lmax_F,
              log_lsel_sd=-2,
              rho_sela=-3,
              rho_selt=-3,
              log_sel_devs=rep(0,length(dat$fyear:dat$lyear)*length(2:10)))
  
  #Parameter names
  parm_names<-names(MakeADFun(dat, par, DLL=TMB_name)$par)
  
  #Param bounds
  lower_bounds<-c(-5,-20,-20,rep(-10,dat$lyear+dat$lage),0,10,-5,-5,-5,-5,-2,0,-10,-10,-10,-10,-10,-10,-10,-10,rep(-20,dat$lyear),-10,-5,-5,rep(-10,length(dat$fyear:dat$lyear)*length(2:10)))
  upper_bounds<-c( 2,  1,  1,rep( 10,dat$lyear+dat$lage),1,25, 2, 2, 2, 2, 5,20,20, 20, 20, 20, 20, 20, 20, 20,rep(  0,dat$lyear),  1, 5, 5,rep(  1,length(dat$fyear:dat$lyear)*length(2:10)))
  
  if(surv==TRUE & Est_Rsd==TRUE){
    fixed<-list(log_M_scalar=factor(NA),steepness=factor(NA),log_cv_fishery=factor(NA),log_cv_fishery_CPUE=factor(NA))  
    
    if(REML==TRUE){reffects=c("log_recruit_devs","log_sel_devs","log_fint","log_R0_FLA","log_q","B1","B2","B3","B4","B5","B6","FIM_sellogis_k","FIM_sellogis_midpt")
    } else if(REML==FALSE){reffects=c("log_recruit_devs","log_sel_devs")}
    
  } else if(surv==TRUE & Est_Rsd==FALSE){
    fixed<-list(log_M_scalar=factor(NA),steepness=factor(NA),log_sigma_rec=factor(NA),log_cv_fishery=factor(NA),log_cv_fishery_CPUE=factor(NA))  
    
    if(REML==TRUE){reffects=c("log_recruit_devs","log_sel_devs","log_fint","log_R0_FLA","log_q","B1","B2","B3","B4","B5","B6","FIM_sellogis_k","FIM_sellogis_midpt")
    } else if(REML==FALSE){reffects=c("log_recruit_devs","log_sel_devs")}
        
  } else if(surv==FALSE & Est_Rsd==TRUE){
   fixed<-list(log_M_scalar=factor(NA),log_q_FIM=factor(NA),steepness=factor(NA),log_cv_fishery=factor(NA),log_cv_fishery_CPUE=factor(NA),
               log_cv_FIM_CPUE=factor(NA),FIM_sellogis_k=factor(NA),FIM_sellogis_midpt=factor(NA),ltheta_FIM=factor(NA))
                            
    if(REML==TRUE){reffects=c("log_recruit_devs","log_sel_devs","log_fint","log_R0_FLA","log_q","B1","B2","B3","B4","B5","B6")
    } else if(REML==FALSE){reffects=c("log_recruit_devs","log_sel_devs")} 
   
  } else if(surv==FALSE & Est_Rsd==FALSE){
   fixed<-list(log_M_scalar=factor(NA),log_q_FIM=factor(NA),steepness=factor(NA),log_sigma_rec=factor(NA),log_cv_fishery=factor(NA),
               log_cv_fishery_CPUE=factor(NA),log_cv_FIM_CPUE=factor(NA),FIM_sellogis_k=factor(NA),FIM_sellogis_midpt=factor(NA),
               ltheta_FIM=factor(NA))
  
    if(REML==TRUE){reffects=c("log_recruit_devs","log_sel_devs","log_fint","log_R0_FLA","log_q","B1","B2","B3","B4","B5","B6")
    } else if(REML==FALSE){reffects=c("log_recruit_devs","log_sel_devs")} 
  }
  
  l<-lower_bounds[-which(parm_names %in% c(names(fixed),reffects))]
  u<-upper_bounds[-which(parm_names %in% c(names(fixed),reffects))]
  
  SCAA <- MakeADFun(dat, par, DLL=TMB_name, map=fixed, random=reffects);

  counter<-1  
  tryCatch({
   SCAA_fit <- fit_tmb(obj=SCAA, startpar=SCAA$par, lower=l, upper=u, newtonsteps=1, getsd=TRUE,bias.correct=TRUE,getHessian=TRUE)
  }, error=function(e){
  counter<<-0
  SCAA_fit<<-list(NA)
  })
 
   convcounter<-1
   if(is.null(SCAA_fit$hessian)){
   jfactor<-10
   while(is.null(SCAA_fit$hessian) & convcounter < 2){
    if(Est_Rsd==TRUE){
    par <- list(log_M_scalar=log(0.099),log_q=jitter(-12.65,jfactor),log_q_FIM=jitter(-14.275, jfactor),
                log_recruit_devs=jitter(rep(0,dat$lyear+dat$lage),jfactor),steepness=0.99,log_R0_FLA=jitter(17.3204,jfactor),
                log_sigma_rec=jitter(log(0.3),jfactor),log_cv_fishery=log(0.05),log_cv_fishery_CPUE=log(0.25),log_cv_FIM_CPUE=jitter(-1.84,jfactor),
                FIM_sellogis_k=jitter(2,jfactor),FIM_sellogis_midpt=jitter(2,jfactor),ltheta=jitter(1,jfactor),ltheta_FIM=jitter(3,jfactor),
                B1=jitter(3.17,jfactor),B2=jitter(-4.18,jfactor),B3=jitter(1.04,jfactor),B4=jitter(1.02,jfactor),B5=jitter(-4.9,jfactor),
                B6=jitter(1.11,jfactor),log_fint=jitter(lmax_F,jfactor),log_lsel_sd=jitter(-2,jfactor),rho_sela=rnorm(1,-1.5,1), rho_selt=rnorm(1,-1.5,1),
                log_sel_devs=jitter(rep(0,length(dat$fyear:dat$lyear)*length(2:10)),jfactor))
    } else if(Est_Rsd==FALSE){
    par <- list(log_M_scalar=log(0.099),log_q=jitter(-12.65,jfactor),log_q_FIM=jitter(-14.275, jfactor),
                log_recruit_devs=jitter(rep(0,dat$lyear+dat$lage),jfactor),steepness=0.99,log_R0_FLA=jitter(17.3204,jfactor),
                log_sigma_rec=log(0.3),log_cv_fishery=log(0.05),log_cv_fishery_CPUE=log(0.25),log_cv_FIM_CPUE=jitter(-1.84,jfactor),
                FIM_sellogis_k=jitter(2,jfactor), FIM_sellogis_midpt=jitter(2,jfactor),ltheta=jitter(1,jfactor),ltheta_FIM=jitter(3,jfactor),
                B1=jitter(3.17,jfactor),B2=jitter(-4.18,jfactor),B3=jitter(1.04,jfactor),B4=jitter(1.02,jfactor),B5=jitter(-4.9,jfactor),
                B6=jitter(1.11,jfactor),log_fint=jitter(lmax_F,jfactor),log_lsel_sd=jitter(-2,jfactor),rho_sela=rnorm(1,-1.5,1), rho_selt=rnorm(1,-1.5,1),
                log_sel_devs=jitter(rep(0,length(dat$fyear:dat$lyear)*length(2:10)),jfactor))
    }
    SCAA <- MakeADFun(dat, par, DLL=TMB_name, map=fixed, random=reffects)
     tryCatch({
      SCAA_fit <- fit_tmb(obj=SCAA, startpar=SCAA$par, lower=l, upper=u, newtonsteps=1, getsd=TRUE,bias.correct=TRUE,getHessian=TRUE)
     }, error=function(e){
     counter<<-0
     SCAA_fit<<-list(NA)
     })
    convcounter<-sum(convcounter, 1)
   }
  }
  
  if (substr(model,8,8)=="1"){
   dname<-paste0(OM,substr(model,5,7),".",substr(model,8,8),"_",i)
  } else if (substr(model,8,8) %in% c("2","3")){
   dname<-paste0(OM,substr(model,5,7),".",as.numeric(substr(model,8,8))-1,"_",i)
  }
  #Extracting arguments to make a dat file for TMB

  dir.create(paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",dname,"/",model))
  
  if(counter==1){
  summ_sdr<-summary(SCAA_fit$SD)

  repor<-SCAA$report(SCAA$env$last.par.best)
  save(repor, file=paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",dname,"/",model,"/report.RData"))
  save(summ_sdr, file=paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",dname,"/",model,"/summ_sdr.RData")) 
  save(SCAA_fit, file=paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",dname,"/",model,"/SCAA_fit.RData"))  
  }
  write(counter, file=paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",dname,"/",model,"/counter.txt"))
  write(convcounter, file=paste0("/blue/edvcamp/nfisch/Chapter_4/OMs/",dname,"/",model,"/convcounter.txt"))
 } 
}

#Call to fit function

#OM11
#run_TMB(N=1:100,model="OM11SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=45:50,model="OM11SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM12
#run_TMB(N=1:100,model="OM12SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM12SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM21
#run_TMB(N=1:100,model="OM21SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=44:44,model="OM21SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM22
#run_TMB(N=1:100,model="OM22SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM22SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM31
#run_TMB(N=1:100,model="OM31SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=51:55,model="OM31SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM32
#run_TMB(N=1:100,model="OM32SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM32SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM41
#run_TMB(N=1:100,model="OM41SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM41SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM41SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM41SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=86:90,model="OM41SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM41SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM41SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=46:50,model="OM41SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM41SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM42
#run_TMB(N=1:100,model="OM42SM11EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM12EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM13EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM21EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM22EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM23EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM31EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM32EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM42SM33EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");


#The information content runs 
#OM11
#run_TMB(N=1:100,model="OM11SM51EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM52EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM11SM53EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM21
#run_TMB(N=1:100,model="OM21SM51EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM52EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM21SM53EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM31
#run_TMB(N=1:100,model="OM31SM51EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM52EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM31SM53EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

#OM41
#run_TMB(N=1:100,model="OM41SM51EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
#run_TMB(N=1:100,model="OM41SM52EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");
run_TMB(N=1:100,model="OM41SM53EM3111_wsurv", surv=TRUE, Est_Rsd=TRUE, REML=FALSE, PE="PE");

