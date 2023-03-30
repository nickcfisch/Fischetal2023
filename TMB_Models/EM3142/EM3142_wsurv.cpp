//Example Catch at age assessment created by Nicholas Fisch to be used in simulation study

//EM3112 (Xu sel devs, Xu M devs, with Logistic-normal)

#include <TMB.hpp>

//this namespace is used for printing during iterations using CppAD::PrintFor
 namespace CppAD {
  void PrintFor(const char* before, const double& var) {  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  //TABLE OF CONTENTS
  //I. DATA INPUTS
  //II. PARAMETER DECLARATION
  //III. SETTING INTITAL VALUES FOR DATA AND BACK TRASFORMING PARAMETERS ("Pre-function" section)
  //IV. FUNCTIONS AND CALCULATIONS
  //    1. Selectivity
  //    2. Mortality
  //    3. Population
  //    4. Expected Catch and Catch at Age
  //    5. Expected Survey Catch at Age and index
  //    6. Projections
  //    7. Objective Functions
  //V. REPORT SECTION

  //I. DATA INPUTS
  //---------------------------------------------------------------------------------------------
  DATA_INTEGER(fyear); //First year
  DATA_INTEGER(lyear); //Last year
  DATA_INTEGER(fage); //First age
  DATA_INTEGER(lage); //Last age

  DATA_VECTOR(years); //Vector of years
  DATA_VECTOR(ages); //Vector of ages

  DATA_VECTOR(obs_harv);              //Observed harvest
  DATA_VECTOR(obs_fishery_cpue);      //Observed fishery CPUE
  DATA_MATRIX(obs_fishery_comp);      //Observed fishery Composition
  DATA_VECTOR(SS_fishery);            //Sample size for fishery compositions

  DATA_VECTOR(obs_FIM_CPUE);          //Observed survey index 
  DATA_MATRIX(obs_FIM_comp);          //Observed survey compositions
  DATA_VECTOR(SS_FIM);                //Sample size for survey compositions

  DATA_VECTOR(Fecund_aa);             //Fecundity at age vector 
  DATA_VECTOR(Wtage);                 //Weight at age vector

  //II. PARAMETER DECLARATION
  //---------------------------------------------------------------------------------------------
  PARAMETER(log_M_scalar);  //M scalar for Lorenzen M
  PARAMETER(log_q);         //Fishery Catchability
  PARAMETER(log_q_FIM);     //Survey Catchability
  
  PARAMETER_VECTOR(log_recruit_devs); //Log scale recruitment deviations
  PARAMETER(steepness);               //Steepness for recruitment
  PARAMETER(log_R0_FLA);              //Unfished recruitment
  PARAMETER(log_sigma_rec);           //Log sd for recruitment

  PARAMETER(log_cv_fishery);         //Log cv for fishery catch
  PARAMETER(log_cv_fishery_CPUE);    //Log cv for fishery cpue
  PARAMETER(log_cv_FIM_CPUE);        //Log cv for survey cpue

  PARAMETER(FIM_sellogis_k);         //logistic selectivity slope parameter for Survey 
  PARAMETER(FIM_sellogis_midpt);     //logistic selectivity midpoint parameter for Survey 
  
  PARAMETER(phi);                //Phi for fishery LN
  PARAMETER(lsig_LN);            //log sigma for fishery LN

  PARAMETER(phi_surv);                //Phi for survey LN
  PARAMETER(lsig_LN_surv);            //log sigma for survey LN
  
  PARAMETER(B1);                     //Double norml sel pars
  PARAMETER(B2);     
  PARAMETER(B3);         
  PARAMETER(B4);     
  PARAMETER(B5);         
  PARAMETER(B6);     

  PARAMETER_VECTOR(log_fint);       //Log scale fishing intensities (fully selected fishing mortalities)
  
  PARAMETER(log_lsel_sd);
  PARAMETER(rho_sela);                   //rho for age
  PARAMETER(rho_selt);                   //rho for time
  PARAMETER_VECTOR(log_sel_devs);

  PARAMETER(log_lM_sd);
  PARAMETER(rho_Ma);                   //rho for age
  PARAMETER(rho_Mt);                   //rho for time
  PARAMETER_VECTOR(log_M_devs);
  //---------------------------------------------------------------------------------------------

//DERIVED QUANTITIES
 
//back-transform log-scale params
  Type M_scalar=exp(log_M_scalar);
  Type q_scalar=exp(log_q);
  Type q_FIM=exp(log_q_FIM);
  Type R0=exp(log_R0_FLA);
  Type sd_rec=exp(log_sigma_rec);                //sigma for recruit deviations
  Type cv_fishery=exp(log_cv_fishery);           //sigma for fishery catch
  Type cv_fishery_CPUE=exp(log_cv_fishery_CPUE); //cv for fishery CPUE
  Type cv_FIM_CPUE=exp(log_cv_FIM_CPUE);         //sigma for survey CPUE

  vector<Type> q_time(years.size());                       //catchability, (fyear,lyear)
  vector<Type> fishery_sel(ages.size());                   //Fishery selectivity (fage,lage)
  vector<Type> fishery_sel_norm(ages.size());              //Fishery selectivity nonrmalized to max at 1 (fage,lage)
  vector<Type> FIM_sel(ages.size());                       //FIM selectivity (fage,lage)

  matrix<Type> M(years.size(),ages.size());                //Instantaneous natural mortality matrix
  matrix<Type> Sel(years.size(),ages.size());              //Selectivity
  matrix<Type> F(years.size(),ages.size());                //Instantaneous fishing mortality matrix
  matrix<Type> Z(years.size(),ages.size());                //Total mortality matrix
  matrix<Type> S(years.size(),ages.size());                //Annual survival matrix
  matrix<Type> A(years.size(),ages.size());                //Annual mortality matrix

  matrix<Type> N(years.size()+1,ages.size());	           //Predicted abundance at age, +1 is to incorporate final catch 
  matrix<Type> Nbar(years.size(),ages.size());	           //Nbar

  matrix<Type> pred_caa(years.size(),ages.size());           //Predicted fishery catch at age
  matrix<Type> pred_fishery_comp(years.size(),ages.size());  //Predicted age composition from the fishery catch
  vector<Type> pred_harv(years.size());                      //Predicted Harvest weight (kg), (fyear,lyear) 
  matrix<Type> pred_harv_aa(years.size(),ages.size());
  vector<Type> pred_totcatch_fishery(years.size());          //Predicted Fishery total catch, used as denominator to get age composition, (fyear,lyear)
  vector<Type> Biomass(years.size());                       //Biomass 
  matrix<Type> Biomass_aa(years.size(),ages.size());
  vector<Type> U(years.size());                             //Exploitation rate 

  matrix<Type> pred_FIM_caa(years.size(),ages.size());   //Predicted FIM catch at age
  matrix<Type> pred_FIM_comp(years.size(),ages.size());  //Predicted age composition from the survey
  vector<Type> pred_totcatch_FIM(years.size());          //Predicted FIM total catch, used as denominator to get age composition, (fyear,lyear)
  vector<Type> pred_fishery_cpue(years.size());          //Predicted fishery cpue, (fyear,lyear)
  matrix<Type> pred_fishery_cpue_aa(years.size(),ages.size());
  
  vector<Type> log_rec_devs(years.size()+ages.size());    //rec devs vector that adds in the current year, (fage,lyear+lage)

  vector<Type> spbiomass(years.size()+1);                  //Spawning biomass, (fyear,lyear+1)
  vector<Type> Depletion(years.size()+1);                  //Depletion, (fyear,lyear+1)
  matrix<Type> spbiomass_aa(years.size()+1,ages.size());   //Spawning biomass, (fyear,lyear+1)
  vector<Type> N0_FLA_age(ages.size());                    //Unfished numbers at, (fage,lage)
  vector<Type> lxo(ages.size());                           //Unfished numbers at age, (fage,lage)
  Type SSB0_FLA;                                           //Unfished spawning biomass
  
  vector<Type> LN_seg1(ages.size());
  matrix<Type> covmat(ages.size(),ages.size());
  matrix<Type> covmat_rho(ages.size(),ages.size());
  matrix<Type> indicator_mat(ages.size(),ages.size());
  matrix<Type> Kmat(ages.size()-1,ages.size());
  matrix<Type> Vmat(ages.size()-1,ages.size()-1);
  matrix<Type> Vinv(ages.size()-1,ages.size()-1);
  vector<Type> wts(years.size());
  vector<Type> ww(ages.size()-1);
  vector<Type> sig_LN_yr(years.size());

  Type SS_FIM_sum=0;
  vector<Type> LN_surv_seg1(ages.size());
  matrix<Type> covmat_surv(ages.size(),ages.size());
  matrix<Type> covmat_surv_rho(ages.size(),ages.size());
  matrix<Type> Vmat_surv(ages.size()-1,ages.size()-1);
  matrix<Type> Vinv_surv(ages.size()-1,ages.size()-1);
  vector<Type> wts_surv(years.size());
  vector<Type> ww_surv(ages.size()-1);
  vector<Type> sig_LN_surv_yr(years.size());
  
  Type rho_sela_logis; 
  Type rho_selt_logis;
//  matrix<Type> indicator_mat_a(ages.size(),ages.size());
  matrix<Type> indicator_mat_yr(years.size(),years.size());
//  matrix<Type> Ra_sel(ages.size(),ages.size());
  matrix<Type> Rt_sel(years.size(),years.size());
//  matrix<Type> Rtotal_sel(ages.size()*years.size(),ages.size()*years.size());
 
  matrix<Type> indicator_mat_a(9,9);                          //only doing ages 2-10
  matrix<Type> Ra_sel(9,9);
  matrix<Type> Rtotal_sel(9*years.size(),9*years.size());
  matrix<Type> varcov_sel(9*years.size(),9*years.size());

  Type rho_Ma_logis; 
  Type rho_Mt_logis;
//  matrix<Type> Ra_M(ages.size(),ages.size());
  matrix<Type> Rt_M(years.size(),years.size());
//  matrix<Type> Rtotal_M(ages.size()*years.size(),ages.size()*years.size());
 
  matrix<Type> Ra_M(9,9);
  matrix<Type> Rtotal_M(9*years.size(),9*years.size());
  matrix<Type> varcov_M(9*years.size(),9*years.size());

  //Growth params for Lorenzen M
  Type Linf=Type(85.64);       //L-Infinity (cm)
  Type K=Type(0.19);           //Brody Growth Coefficient
  Type tnot=Type(-0.39);      //T-not
  vector<Type> Lt(ages.size());
  
  Type peak2; 
  Type t1; 
  Type t2; 
  vector<Type> j1(ages.size());
  vector<Type> j2(ages.size());
  vector<Type> asc(ages.size());
  vector<Type> dsc(ages.size());

  //Likelihood components
  Type L1=Type(0);
  vector<Type> L2_comp(years.size());
  Type L2=Type(0);
  Type L3=Type(0);
  vector<Type> L4_comp(years.size());
  Type L4=Type(0);
  Type L5=Type(0);

  int i; 
  int j; 
  int k; 
   
  Type NLL = Type(0);
  Type NPRAND = Type(0);
  Type JNLL = Type(0);
  
//////////////////////////////////////
///////////////////////////
//Actual Model
///////////////////////////
/////////////////////////////////////

///////////////////////////
//FISHERY SELECTIVITY
///////////////////////////
  Type Amin=0;
  Type Amax=20;

  for(j=0;j<=ages.size()-1;j++){
   //Double normal sel
   peak2=B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)));
   t1=exp(-pow(Amin-B1,2)/exp(B3));
   t2=exp(-pow(Amax-peak2,2)/exp(B4));
  
   j1(j)=pow((1+exp(-20*((int(j)-B1)/(1+fabs(int(j)-B1))))),-1);
   j2(j)=pow((1+exp(-20*((int(j)-peak2)/(1+fabs(int(j)-peak2))))),-1);

   //Long form
   asc(j)=pow(1+exp(-B5),-1)+(1-pow(1+exp(-B5),-1))*((exp(-pow(int(j)-B1,2)/exp(B3))-t1)/(1-t1));
   dsc(j)=1+(pow(1+exp(-B6),-1)-1)*((exp(-(int(j)-peak2)/exp(B4))-1)/(t2-1));

   //Actual Sel
  fishery_sel(j)=asc(j)*(1-j1(j))+j1(j)*((1-j2(j))+j2(j)*dsc(j));
  }

  for(i=0;i<=years.size()-1;i++){
   for(j=2;j<=10;j++){
    Sel(i,j)=fishery_sel(j)*exp(log_sel_devs(i*9+(j-2))); 	
   }
   for(j=0;j<=1;j++){
    Sel(i,j)=fishery_sel(j)*exp(log_sel_devs(i*9)); 	  //fixing sel dev below age 2
   }
   for(j=11;j<=lage;j++){
    Sel(i,j)=fishery_sel(j)*exp(log_sel_devs(i*9+8)); 	//fixing sel dev above age 10
   }
  }
    
//////////////////////
//MORTALITY
//////////////////////
 //First getting size at age for Lorenzen M
  for(j=fage;j<=lage;j++){
   Lt(j)=Linf*(Type(1)-exp(-K*(int(j)-tnot)));
  }

 //Filling in natural mortality matrix
  for(i=0;i<=years.size()-1;i++){
   M(i,0)=2.0*exp(log_M_devs(i*9));   //Filling in age 0 M 
   M(i,1)=1.2*exp(log_M_devs(i*9));   //Filling in age 1 M, rest is based on Lorenzen
   for(j=2;j<=10;j++){
    M(i,j)=(M_scalar*Lt(7)/Lt(j))*exp(log_M_devs(i*9+(j-2))); 	
   }
   for(j=11;j<=lage;j++){
    M(i,j)=(M_scalar*Lt(7)/Lt(j))*exp(log_M_devs(i*9+8)); 	
   }
  }
  
  for(i=0;i<=years.size()-1;i++){
   for(j=fage;j<=lage;j++){
    Sel(i,j)=Sel(i,j)/Sel(i,3);               //Setting time-varying sel relative to age 3
    F(i,j) = Sel(i,j)*exp(log_fint(i)); //setting year and age specific fishing mortality as the product of selectivity and year specific fishing intensity
    Z(i,j) = M(i,j) + F(i,j);                 //Total instantaneous mortality
    S(i,j) = exp(Type(-1)*Z(i,j));                //Annual survival
    A(i,j)=Type(1.0)-S(i,j);                        //Annual mortality
   }
  }
 
//////////////////////
//POPULATION
//////////////////////
  //Unfished Spawning biomass calcs
  lxo(fage)=Type(1.0);
  for(j=fage+1;j<=lage;j++){
//   lxo(j)=lxo(j-1)*exp(-M(0,j-1));   //cumulative survival
   lxo(j)=lxo(j-1)*exp(-(M_scalar*Lt(7)/Lt(j-1)));   //cumulative survival
  }
//  lxo(lage)=lxo(lage)/(1-exp(-M(0,lage)));
  lxo(lage)=lxo(lage)/(1-exp(-(M_scalar*Lt(7)/Lt(lage))));
  N0_FLA_age=R0*lxo;
  SSB0_FLA=(N0_FLA_age*Fecund_aa).sum();   //Numbers at age * Fecundity

  //Filling in log_rec_devs
  for(i=0;i<=years.size()+lage-1;i++){
   log_rec_devs(i)=log_recruit_devs(i);
  }
  log_rec_devs(lyear+lage)=Type(0);
  
  //Abundance at age in the first year
  for(j=0;j<=ages.size()-1;j++){
   N(fyear-1,j)=R0*exp(log_recruit_devs(j))*lxo(j);           //set abundance in the fist year f
//  N(fyear-1,j)=R0*exp(log_recruit_devs(j)-Type(0.5)*pow(sd_rec,Type(2)))*lxo(j);           //set abundance in the fist year f
  }
  
  //Spawning Biomass in the first year
  for(j=0;j<=ages.size()-1;j++){
   spbiomass_aa(fyear-1,j)=N(fyear-1,j)*Fecund_aa(j);  // Getting spawning biomass for the first year (acounting for rec dev, unlike SSB0_FLA)
  }
  spbiomass(fyear-1)=spbiomass_aa.row(fyear-1).sum();
  
//Population loop
  for(i=fyear;i<=years.size();i++){  //TMB starts at zero so this is 81 years (or 41), starting at i+1
    for(j=fage+1;j<=lage;j++){
      N(i,j)=N(i-1,j-1)*S(i-1,j-1);
     }
   N(i,lage)+=N(i-1,lage)*S(i-1,lage);   //Plus group
  for(j=0;j<=ages.size()-1;j++){
   spbiomass_aa(i,j)=N(i,j)*Fecund_aa(j);  // Getting spawning biomass for the first year (acounting for rec dev, unlike SSB0_FLA)
  }
   spbiomass(i)=spbiomass_aa.row(i).sum();
  N(i,fage)=((Type(4.)*steepness*R0*spbiomass(i))/(SSB0_FLA*(Type(1.)-steepness)+spbiomass(i)*(Type(5.)*steepness-Type(1.))))*exp(log_rec_devs(i+lage));  //Recruitment
//   N(i,fage)=((Type(4.)*steepness*R0*spbiomass(i))/(SSB0_FLA*(Type(1.)-steepness)+spbiomass(i)*(Type(5.)*steepness-Type(1.))))*exp(log_rec_devs(i+lage)-Type(0.5)*pow(sd_rec,Type(2)));  //Recruitment
 }
 
////////////////
//CATCH
////////////////

  for(i=0;i<=years.size()-1;i++){  //TMB starts at zero so this is 81 years (or 41), starting at i+1
   for(j=fage;j<=lage;j++){  
    pred_caa(i,j)=F(i,j)/Z(i,j)*A(i,j)*N(i,j); //Baranov catch equation for predicting catch at age
   }
  }

  for(i=0;i<=lyear-1;i++){
   q_time(i)=exp(log_q);
   pred_totcatch_fishery(i)=pred_caa.row(i).sum();                        //total predicted catch by year
   pred_fishery_comp.row(i)=pred_caa.row(i)/(pred_totcatch_fishery(i)+Type(0.0001));  //calculating predicted catch age composition
   for(j=fage;j<=lage;j++){  
    pred_harv_aa(i,j)=pred_caa(i,j)*Wtage(j);                              //Predicted total harvest at age weight each year
    Biomass_aa(i,j)=N(i,j)*Wtage(j);
    Nbar(i,j)=(N(i,j)*((Type(1.0)-exp(-Type(1.0)*Z(i,j)))))*(Type(1.0)/Z(i,j));
    pred_fishery_cpue_aa(i,j)=q_time(i)*(Nbar(i,j)*Wtage(j))*Sel(i,j);
   }
   pred_harv(i)=pred_harv_aa.row(i).sum();                              //Predicted total harvest weight each year
   Biomass(i)=Biomass_aa.row(i).sum();                                    //Biomass each year
   U(i)=pred_harv(i)/Biomass(i);                                          //Exploitation rate
   Depletion(i)=spbiomass(i)/SSB0_FLA;
   pred_fishery_cpue(i)=pred_fishery_cpue_aa.row(i).sum();
  }
    
  Depletion(lyear)=spbiomass(lyear)/SSB0_FLA;
////////////////////
//SURVEY
////////////////////
  for (j=fage;j<=lage;j++){
   FIM_sel(j)=Type(1)/(Type(1)+exp(-FIM_sellogis_k*(int(j)-FIM_sellogis_midpt)));  //Logistic Selectivity
  }
  FIM_sel=FIM_sel/max(FIM_sel);

  for (i=0;i<=lyear-1;i++){
   for(j=fage;j<=lage;j++){  
    pred_FIM_caa(i,j)=exp(log_q_FIM)*FIM_sel(j)*N(i,j);	//Predicted survey catch-at-age each year 
   }
   pred_totcatch_FIM(i)=pred_FIM_caa.row(i).sum();          // total predicted survey cpue
   pred_FIM_comp.row(i)=pred_FIM_caa.row(i)/(pred_totcatch_FIM(i)+Type(0.0001));  //calculating predicted survey catch age composition
  }

/////////////////////////
//Objective function
/////////////////////////

//Composition Likelihoods
  Type pi = Type(3.14159265358979323844);
  Type sig_LN=exp(lsig_LN);
  Type sig_LN_surv=exp(lsig_LN_surv);             //Exponentiating sigma
  for(i=0;i<=lyear-1;i++){
   if(SS_fishery(i)>0){
    wts(i)=sqrt((SS_fishery.sum()/lyear)/SS_fishery(i));             //weights for each year of composition
    sig_LN_yr(i) = wts(i)*sig_LN; 
   }
  }

  for(i=(lyear/2);i<=lyear-1;i++){
   SS_FIM_sum += SS_FIM(i);
  }

  for(i=(lyear/2);i<=lyear-1;i++){
   if(SS_FIM(i)>0){
    wts_surv(i)=sqrt((SS_FIM_sum/(lyear/Type(2)))/SS_FIM(i));             //weights for each year of composition
    sig_LN_surv_yr(i) = wts_surv(i)*sig_LN_surv;     //Getting sigmas for each year 
   }
  }

  Type phi_new=atan(phi)*(1+1)/pi + 0.5*(-1+1);
  Type phi_surv_new=atan(phi_surv)*(1+1)/pi + 0.5*(-1+1);

   for(i=fage;i<=lage;i++){
    for(j=fage;j<=lage;j++){
     indicator_mat(i,j) = abs(int(i+1)-int(j+1));       //dummy matrix only needed for next calculation
     covmat_rho(i,j) =  pow(phi_new,indicator_mat(i,j));       //covariance matrix of only rho
     covmat_surv_rho(i,j) =  pow(phi_surv_new,indicator_mat(i,j));       //covariance matrix of only rho
    }
   }
   
  for(i=fage;i<=lage-1;i++){
   Kmat(i,i)=Type(1);                               //Making identity matrix
   Kmat(i,lage)=-Type(1);                           //adding a column of -1's
  }

//Now code for NLL of logistic normal
// Fishery
  for(i=0;i<=lyear-1;i++){
   if(SS_fishery(i)>0){
    for(j=fage;j<=lage;j++){
      for(k=fage;k<=lage;k++){
       covmat(j,k) = covmat_rho(j,k)*pow(sig_LN_yr(i),Type(2));   //bringing in sigma to Cov matrix
	  }
	}
    Vmat = Kmat * (covmat * Kmat.transpose());        //Vmat calc
    Vinv = atomic::matinv(Vmat);                            //inverse of Vmat
    for(j=fage;j<=lage-1;j++){
     ww(j) = log(obs_fishery_comp(i,j)/obs_fishery_comp(i,lage)) - log(pred_fishery_comp(i,j)/pred_fishery_comp(i,lage));
	}
    LN_seg1 = ww.matrix().transpose() * Vinv * ww.matrix();

    L2_comp(i) = Type(0.5)*(lage-fage)*log(Type(2)*pi)+log(obs_fishery_comp.array()).row(i).sum()+Type(0.5)*atomic::logdet(Vmat)+(lage-fage)*log(wts(i))+Type(0.5)*(LN_seg1.sum()/pow(wts(i),Type(2)));
   }
  }

  L2=L2_comp.sum();

//Survey 
  for(i=(lyear/2);i<=lyear-1;i++){
   if(SS_FIM(i)>0){
    for(j=fage;j<=lage;j++){
     for(k=fage;k<=lage;k++){
      covmat_surv(j,k) = covmat_surv_rho(j,k)*pow(sig_LN_surv_yr(i),Type(2));     //adding sig^2 into covariance matrix
	 }
	}
   Vmat_surv = Kmat * (covmat_surv * Kmat.transpose());        //Vmat calc
   Vinv_surv = atomic::matinv(Vmat_surv);                            //inverse of Vmat
   
    for(j=fage;j<=lage-1;j++){
     ww_surv(j) = log(obs_FIM_comp(i,j)/obs_FIM_comp(i,lage)) - log(pred_FIM_comp(i,j)/pred_FIM_comp(i,lage));
    }
    LN_surv_seg1 = ww_surv.matrix().transpose() * Vinv_surv * ww_surv.matrix();

    L4_comp(i) = Type(0.5)*(lage-fage)*log(Type(2)*pi)+log(obs_FIM_comp.array()).row(i).sum()+Type(0.5)*atomic::logdet(Vmat_surv)+(lage-fage)*log(wts_surv(i))+Type(0.5)*(LN_surv_seg1.sum()/pow(wts_surv(i),Type(2)));
   }
  }
  
  L4=L4_comp.sum();
  
//LIKELIHOODS
  //Catch Likelihood
  // Normal by CV
  for(i=0;i<=lyear-1;i++){
   L1 += log(cv_fishery*pred_harv(i)+Type(0.0001))+Type(0.5)*pow((obs_harv(i)-pred_harv(i))/(cv_fishery*pred_harv(i)+Type(0.0001)),Type(2));   //Normal likelihood using CV
  } 

//Fishery CPUE
  for(i=0;i<=lyear-1;i++){
   L3 += log(cv_fishery_CPUE*pred_fishery_cpue(i)+Type(0.0001))+Type(0.5)*pow((obs_fishery_cpue(i)-pred_fishery_cpue(i))/(cv_fishery_CPUE*pred_fishery_cpue(i)+Type(0.0001)),Type(2)); //Likelihood for Fishery CPUE
  }

//Fisheries-independent likelihoods
  for(i=(lyear/2);i<=lyear-1;i++){
   L5 += log(cv_FIM_CPUE*pred_totcatch_FIM(i)+Type(0.0001))+Type(0.5)*pow((obs_FIM_CPUE(i)-pred_totcatch_FIM(i))/(cv_FIM_CPUE*pred_totcatch_FIM(i)+Type(0.0001)),Type(2)); //Likelihood for FIM CPUE
  }
  
  using namespace density;

  rho_sela_logis=1/(1+exp(-rho_sela/1));
  rho_selt_logis=1/(1+exp(-rho_selt/1));

//   for(i=fage;i<=lage;i++){
//    for(j=fage;j<=lage;j++){
   for(i=fage;i<=8;i++){
    for(j=fage;j<=8;j++){
     indicator_mat_a(i,j) = abs(int(i+1)-int(j+1));       //dummy matrix only needed for next calculation
     Ra_sel(i,j) =  pow(rho_sela_logis,indicator_mat_a(i,j));       //covariance matrix of only rho
    }
   }

   for(i=0;i<=lyear-1;i++){
    for(j=0;j<=lyear-1;j++){
     indicator_mat_yr(i,j) = abs(int(i+1)-int(j+1));       //dummy matrix only needed for next calculation
     Rt_sel(i,j) =  pow(rho_selt_logis,indicator_mat_yr(i,j));       //covariance matrix of only rho
    }
   }

//  Rtotal_sel=kronecker(Ra_sel,Rt_sel);
  Rtotal_sel=kronecker(Rt_sel,Ra_sel);

  for(i=0;i<=9*years.size()-1;i++){
   for(j=0;j<=9*years.size()-1;j++){
    varcov_sel(i,j)=exp(log_lsel_sd)*Rtotal_sel(i,j)*exp(log_lsel_sd);
   }
  }
  
  //M Devs
  rho_Ma_logis=1/(1+exp(-rho_Ma/1));
  rho_Mt_logis=1/(1+exp(-rho_Mt/1));
  
//   for(i=fage;i<=lage;i++){
//    for(j=fage;j<=lage;j++){
   for(i=fage;i<=8;i++){
    for(j=fage;j<=8;j++){
     indicator_mat_a(i,j) = abs(int(i+1)-int(j+1));       //dummy matrix only needed for next calculation
     Ra_M(i,j) =  pow(rho_Ma_logis,indicator_mat_a(i,j));       //covariance matrix of only rho
    }
   }

   for(i=0;i<=lyear-1;i++){
    for(j=0;j<=lyear-1;j++){
     indicator_mat_yr(i,j) = abs(int(i+1)-int(j+1));       //dummy matrix only needed for next calculation
     Rt_M(i,j) =  pow(rho_Mt_logis,indicator_mat_yr(i,j));       //covariance matrix of only rho
    }
   }

//  Rtotal_M=kronecker(Ra_M,Rt_M);
  Rtotal_M=kronecker(Rt_M,Ra_M);

  for(i=0;i<=9*years.size()-1;i++){
   for(j=0;j<=9*years.size()-1;j++){
    varcov_M(i,j)=exp(log_lM_sd)*Rtotal_M(i,j)*exp(log_lM_sd);
   }
  }
  
  //Recruitment deviations
   NPRAND -= sum(dnorm(log_recruit_devs,0,sd_rec,true)); //Recruitment deviations
   NPRAND += MVNORM(varcov_sel)(log_sel_devs); //Sel deviations
   NPRAND += MVNORM(varcov_M)(log_M_devs); //M deviations

  NLL=L1+L2+L3+L4+L5;
//  NLL=L1+L2+L3;
  JNLL=NPRAND+NLL;

//  std::cout<<"Test of cout - current value of L2 second time is: "<<L2<<std::endl;
//  std::cout<<"Test of cout - current value of L4 second time is: "<<L4<<std::endl;

////////////////////////////
//Report
////////////////////////////

  REPORT(M_scalar);
  REPORT(q_scalar);
  REPORT(q_FIM);
  REPORT(R0);
  REPORT(sd_rec);
  REPORT(cv_fishery);
  REPORT(cv_fishery_CPUE);
  REPORT(cv_FIM_CPUE);
  REPORT(Fecund_aa);
  REPORT(log_recruit_devs);

  REPORT(fishery_sel);
  REPORT(Depletion);
  ADREPORT(Depletion);
  REPORT(Sel);
  ADREPORT(Sel);
  REPORT(N0_FLA_age);
  ADREPORT(N0_FLA_age);
  REPORT(Lt);
  REPORT(lxo);
  REPORT(F);
  ADREPORT(F);
  REPORT(M);
  ADREPORT(M);
  REPORT(Z);
  ADREPORT(Z);
  REPORT(A);
  REPORT(N);
  ADREPORT(N);
  REPORT(U);
  ADREPORT(U);
  REPORT(spbiomass);
  ADREPORT(spbiomass);
  REPORT(obs_harv);
  REPORT(pred_harv);
  REPORT(pred_caa);
  REPORT(pred_fishery_cpue);
  REPORT(obs_fishery_cpue);
  REPORT(pred_totcatch_FIM);
  REPORT(obs_FIM_CPUE);
  REPORT(cv_FIM_CPUE);
  REPORT(pred_fishery_comp);
  REPORT(obs_fishery_comp);
  REPORT(pred_FIM_comp);
  REPORT(obs_FIM_comp);

  REPORT(L1);
  REPORT(L2);
  REPORT(L3);
  REPORT(L4);
  REPORT(L5);
  REPORT(NLL);
  REPORT(NPRAND);

  return JNLL;

}

