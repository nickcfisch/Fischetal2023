//Example Catch at age assessment created by Nicholas Fisch to be used in simulation study

//EM1115 (Constant everything, with Multivariate-Tweedie [MVTW])

//Nov 2022

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
  
  PARAMETER(lphi);                //phi term for MVTW for fishery composition
  PARAMETER(est_psi);             //power term for MVTW for fishery composition
  PARAMETER(lphi_surv);           //phi term for MVTW for survey composition
  PARAMETER(est_psi_surv);        //power term for MVTW for survey composition
  
  PARAMETER(B1);                     //Double norml sel pars
  PARAMETER(B2);     
  PARAMETER(B3);         
  PARAMETER(B4);     
  PARAMETER(B5);         
  PARAMETER(B6);     

  PARAMETER_VECTOR(log_fint);       //Log scale fishing intensities (fully selected fishing mortalities)
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
  Type phi=exp(lphi);                            //Transforming MVTW params
  Type phi_surv=exp(lphi_surv);
  Type pi = Type(3.14159265358979323844);
  Type psi=atan(est_psi)*(2-1)/pi + 0.5*(1+2);
  Type psi_surv=atan(est_psi_surv)*(2-1)/pi + 0.5*(1+2);

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
  
  vector<Type> Beta(years.size());                    //Beta for DM,(fyear,lyear)
  vector<Type> Beta_FIM(years.size());                //Beta for survey DM,(fyear,lyear)

  matrix<Type> DM_component_1(years.size(),ages.size());
  matrix<Type> DM_component_2(years.size(),ages.size());
  matrix<Type> DM_FIM_component_1(years.size(),ages.size());
  matrix<Type> DM_FIM_component_2(years.size(),ages.size());
  
  //Growth params for Lorenzen M
  Type Linf=85.64;       //L-Infinity (cm)
  Type K=0.19;           //Brody Growth Coefficient
  Type tnot= -0.39;      //T-not
  vector<Type> Lt(ages.size());
  
  Type peak2; 
  Type t1; 
  Type t2; 
  vector<Type> j1(ages.size());
  vector<Type> j2(ages.size());
  vector<Type> asc(ages.size());
  vector<Type> dsc(ages.size());

  //Likelihood components
  Type L1;
  Type L2;
  Type L3;
  Type L4;
  Type L5;

  int i; 
  int j; 
  
  Type NLL = 0;
  Type NPRAND = 0;
  Type JNLL = 0;
  
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
   for(j=0;j<=lage;j++){
    Sel(i,j)=fishery_sel(j); 	
   }
  }
  
//////////////////////
//MORTALITY
//////////////////////
 //First getting size at age for Lorenzen M
  for(j=fage;j<=lage;j++){
   Lt(j)=Linf*(1-exp(-K*(int(j)-tnot)));
  }

 //Filling in natural mortality matrix
  for(i=0;i<=years.size()-1;i++){
   M(i,0)=2.0;   //Filling in age 0 M 
   M(i,1)=1.2;   //Filling in age 1 M, rest is based on Lorenzen
   for(j=2;j<=lage;j++){
    M(i,j)=M_scalar*Lt(7)/Lt(j); 	
   }
  }
  
  for(i=0;i<=years.size()-1;i++){
   for(j=fage;j<=lage;j++){
    F(i,j) = Sel(i,j)*exp(log_fint(i)); //setting year and age specific fishing mortality as the product of selectivity and year specific fishing intensity
    Z(i,j) = M(i,j) + F(i,j);                 //Total instantaneous mortality
    S(i,j) = exp(-1.0*Z(i,j));                //Annual survival
    A(i,j)=1.0-S(i,j);                        //Annual mortality
   }
  }
 
//////////////////////
//POPULATION
//////////////////////
  //Unfished Spawning biomass calcs
  lxo(fage)=1.0;
  for(j=fage+1;j<=lage;j++){
   lxo(j)=lxo(j-1)*exp(-M(0,j-1));   //cumulative survival
  }
  lxo(lage)=lxo(lage)/(1-exp(-M(0,lage)));
  N0_FLA_age=R0*lxo;
  SSB0_FLA=(N0_FLA_age*Fecund_aa).sum();   //Numbers at age * Fecundity

  //Filling in log_rec_devs
  for(i=0;i<=years.size()+lage-1;i++){
   log_rec_devs(i)=log_recruit_devs(i);
  }
  log_rec_devs(lyear+lage)=0;
  
  //Abundance at age in the first year
  for(j=0;j<=ages.size()-1;j++){
   N(fyear-1,j)=R0*exp(log_recruit_devs(j))*lxo(j);           //set abundance in the fist year f
//   N(fyear-1,j)=R0*exp(log_recruit_devs(j)-Type(0.5)*pow(sd_rec,Type(2)))*lxo(j);           //set abundance in the fist year f
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
   N(i,fage)=((4.*steepness*R0*spbiomass(i))/(SSB0_FLA*(1.-steepness)+spbiomass(i)*(5.*steepness-1.)))*exp(log_rec_devs(i+lage));  //Recruitment
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
   pred_fishery_comp.row(i)=pred_caa.row(i)/(pred_totcatch_fishery(i)+0.0001);  //calculating predicted catch age composition
   for(j=fage;j<=lage;j++){  
    pred_harv_aa(i,j)=pred_caa(i,j)*Wtage(j);                              //Predicted total harvest at age weight each year
    Biomass_aa(i,j)=N(i,j)*Wtage(j);
    Nbar(i,j)=(N(i,j)*((1.0-exp(-1.0*Z(i,j)))))*(1.0/Z(i,j));
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
   FIM_sel(j)=1/(1+exp(-FIM_sellogis_k*(int(j)-FIM_sellogis_midpt)));  //Logistic Selectivity
  }
  FIM_sel=FIM_sel/max(FIM_sel);

  for (i=0;i<=lyear-1;i++){
   for(j=fage;j<=lage;j++){  
    pred_FIM_caa(i,j)=exp(log_q_FIM)*FIM_sel(j)*N(i,j);	//Predicted survey catch-at-age each year 
   }
   pred_totcatch_FIM(i)=pred_FIM_caa.row(i).sum();          // total predicted survey cpue
   pred_FIM_comp.row(i)=pred_FIM_caa.row(i)/(pred_totcatch_FIM(i)+0.0001);  //calculating predicted survey catch age composition
  }

/////////////////////////
//Objective function
/////////////////////////

//LIKELIHOODS
  //Catch Likelihood
  // Normal by CV
  for(i=0;i<=lyear-1;i++){
   L1 += log(cv_fishery*pred_harv(i)+0.0001)+0.5*pow((obs_harv(i)-pred_harv(i))/(cv_fishery*pred_harv(i)+0.0001),2);   //Normal likelihood using CV
  } 

//MVTW - Multivariate-Tweedie
  for(i=0;i<=lyear-1;i++){
    for(j=0;j<=lage;j++){
     L2 += -1*log(dtweedie(obs_fishery_comp(i,j)*SS_fishery(i),pred_fishery_comp(i,j)*SS_fishery(i),phi,psi));
    }
   }

//Fishery CPUE
  for(i=0;i<=lyear-1;i++){
   L3 += log(cv_fishery_CPUE*pred_fishery_cpue(i)+0.0001)+0.5*pow((obs_fishery_cpue(i)-pred_fishery_cpue(i))/(cv_fishery_CPUE*pred_fishery_cpue(i)+0.0001),2); //Likelihood for Fishery CPUE
  }

//Fisheries-independent likelihoods
  for(i=(lyear/2);i<=lyear-1;i++){
//MVTW - Multivariate-Tweedie
   for(j=0;j<=lage;j++){
    L4 += -1*log(dtweedie(obs_FIM_comp(i,j)*SS_FIM(i),pred_FIM_comp(i,j)*SS_FIM(i),phi_surv,psi_surv));
   }
   L5 += log(cv_FIM_CPUE*pred_totcatch_FIM(i)+0.0001)+0.5*pow((obs_FIM_CPUE(i)-pred_totcatch_FIM(i))/(cv_FIM_CPUE*pred_totcatch_FIM(i)+0.0001),2); //Likelihood for FIM CPUE
  }

  //Recruitment deviations
   NPRAND -= sum(dnorm(log_recruit_devs,0,sd_rec,true)); //Recruitment deviations
   //Alternative way to do it in a loop
//  for(i=0;i<=lyear+lage-1;i++){
//   NPRAND += log(sd_rec*pow(2*3.14159265359,0.5))+0.5*pow(log_recruit_devs(i)/sd_rec,2); //Recruitment deviations
//  }

  NLL=L1+L2+L3+L4+L5;
//  NLL=L1+L2+L3;
  JNLL=NPRAND+NLL;

/////////////////////////
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

