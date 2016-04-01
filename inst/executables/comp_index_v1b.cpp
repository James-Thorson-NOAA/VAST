#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// dlnorm
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  //return 1/(sqrt(2*M_PI)*sd)*exp(-.5*pow((x-mean)/sd,2));
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

// CMP distribution
template<class Type>
Type dCMP(Type x, Type mu, Type nu, int give_log=0, int iter_max=30, int break_point=10){
  // Explicit
  Type ln_S_1 = nu*mu - ((nu-1)/2)*log(mu) - ((nu-1)/2)*log(2*M_PI) - 0.5*log(nu);
  // Recursive
  vector<Type> S_i(iter_max);
  S_i(0) = 1;
  for(int i=1; i<iter_max; i++) S_i(i) = S_i(i-1) * pow( mu/Type(i), nu );
  Type ln_S_2 = log( sum(S_i) );
  // Blend (breakpoint:  mu=10)
  Type prop_1 = invlogit( (mu-break_point)*5 );
  //Type S_comb = prop_1*exp(ln_S_1) + (1-prop_1)*exp(ln_S_2);
  Type log_S_comb = prop_1*ln_S_1 + (1-prop_1)*ln_S_2;
  // Likelihood
  Type loglike = nu*x*log(mu) - nu*lgamma(x+1) - log_S_comb;
  // Return
  if(give_log) return loglike; else return exp(loglike);
}

// Space time
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla;
  using namespace Eigen;
  using namespace density;

  // Dimensions
  DATA_INTEGER(n_i);         // Number of observations (stacked across all years)
  DATA_INTEGER(n_s);         // Number of "strata" (i.e., vectices in SPDE mesh) 
  DATA_INTEGER(n_x);         // Number of real "strata" (i.e., k-means locations) 
  DATA_INTEGER(n_t);         // Number of years
  DATA_INTEGER(n_c);         // Number of categories (e.g., length bins)
  DATA_INTEGER(n_j);         // Number of static covariates
  DATA_INTEGER(n_p);         // Number of dynamic covariates
  DATA_INTEGER(n_k);          // Number of catchability variables
  DATA_INTEGER(n_l);         // Number of indices to post-process
  DATA_INTEGER(n_m);         // Number of range metrics to use (probably 2 for Eastings-Northings)

  // Config
  DATA_IVECTOR( Options_vec );
  // Slot 0 -- Aniso: 0=No, 1=Yes
  // Slot 1 -- DEPRECATED
  // Slot 2 -- AR1 on betas (year intercepts) to deal with missing years: 0=No, 1=Yes
  // Slot 3 -- DEPRECATED
  // Slot 4 -- DEPRECATED
  // Slot 5 -- Upper limit constant of integration calculation for Conway-Maxwell-Poisson density function
  DATA_IVECTOR(FieldConfig);  // Input settings
  DATA_IVECTOR(ObsModel);    // Observation model
  DATA_IVECTOR(Options);    // Reporting options
  // Slot 0: DEPRECATED
  // Slot 1: DEPRECATED
  // Slot 2: DEPRECATED
  // Slot 3: DEPRECATED
  // Slot 4: Calculate mean_D_tl and effective_area_tl

  // Data vectors
  DATA_VECTOR(b_i);       	// Response (biomass) for each observation
  DATA_VECTOR(a_i);       	// Area swept for each observation (km^2)
  DATA_IVECTOR(c_i);         // Category for each observation
  DATA_IVECTOR(s_i)          // Station for each observation
  DATA_IVECTOR(t_i)          // Year for each observation
  DATA_MATRIX(a_xl);		     // Area for each "real" stratum(km^2) in each stratum
  DATA_MATRIX(X_xj);		    // Covariate design matrix (strata x covariate)
  DATA_ARRAY(X_xtp);		    // Covariate design matrix (strata x covariate)
  DATA_MATRIX(Q_ik);        // Catchability matrix (observations x variable)
  DATA_MATRIX(Z_xm);        // Derived quantity matrix

  // SPDE objects
  DATA_STRUCT(spde,spde_t);
  
  // Aniso objects
  DATA_STRUCT(spde_aniso,spde_aniso_t);
  
  // Parameters 
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  //  -- presence/absence
  PARAMETER_MATRIX(beta1_ct);       // Year effect
  PARAMETER_VECTOR(gamma1_j);        // Static covariate effect
  PARAMETER_ARRAY(gamma1_ctp);       // Dynamic covariate effect
  PARAMETER_VECTOR(lambda1_k);       // Catchability coefficients
  PARAMETER(logetaE1);      
  PARAMETER(logetaO1);
  PARAMETER(logkappa1);
  PARAMETER(Beta_mean1);  // mean-reversion for beta1_t
  PARAMETER(logsigmaB1);  // SD of beta1_t (default: not included in objective function)
  PARAMETER(Beta_rho1);  // AR1 for positive catch Epsilon component, Default=0
  PARAMETER(Epsilon_rho1);  // AR1 for presence/absence Epsilon component, Default=0
  PARAMETER(rho_c1);         // AR1 among categories
  // -- Gaussian random fields
  PARAMETER_ARRAY(Omegainput1_sc);      // Expectation
  PARAMETER_ARRAY(Epsiloninput1_sct);   // Annual variation
  //  -- positive catch rates
  PARAMETER_MATRIX(beta2_ct);  // Year effect
  PARAMETER_VECTOR(gamma2_j);        // Covariate effect
  PARAMETER_ARRAY(gamma2_ctp);       // Dynamic covariate effect
  PARAMETER_VECTOR(lambda2_k);       // Catchability coefficients
  PARAMETER(logetaE2);      
  PARAMETER(logetaO2);
  PARAMETER(logkappa2);
  PARAMETER(Beta_mean2);  // mean-reversion for beta2_t
  PARAMETER(logsigmaB2);  // SD of beta2_t (default: not included in objective function)
  PARAMETER(Beta_rho2);  // AR1 for positive catch Epsilon component, Default=0
  PARAMETER(Epsilon_rho2);  // AR1 for positive catch Epsilon component, Default=0
  PARAMETER(rho_c2);         // AR1 among categories
  PARAMETER_VECTOR(logSigmaM);   // Slots: 0=mix1 CV, 1=prob-of-mix1, 2=
  // -- Gaussian random fields
  PARAMETER_ARRAY(Omegainput2_sc);      // Expectation
  PARAMETER_ARRAY(Epsiloninput2_sct);   // Annual variation

  // Indices -- i=Observation; j=Covariate; v=Vessel; t=Year; s=Stratum
  int i,j,v,t,s,c;
  
  // Objective function
  vector<Type> jnll_comp(13);    // Slot 12=optional hyperparameters
  jnll_comp.setZero();
  Type jnll = 0;                
  
  // Derived parameters
  Type pi = 3.141592;
  Type logtauE1 = logetaE1 - logkappa1;
  Type logtauO1 = logetaO1 - logkappa1;
  Type kappa1_pow2 = exp(2.0*logkappa1);
  Type kappa1_pow4 = kappa1_pow2*kappa1_pow2;
  Type SigmaE1 = 1 / sqrt(4*pi*exp(2*logtauE1)*exp(2*logkappa1));
  Type SigmaO1 = 1 / sqrt(4*pi*exp(2*logtauO1)*exp(2*logkappa1));
  Type Range_raw1 = sqrt(8) / exp( logkappa1 );   // Range = approx. distance @ 10% correlation

  Type logtauE2 = logetaE2 - logkappa2;
  Type logtauO2 = logetaO2 - logkappa2;
  Type kappa2_pow2 = exp(2.0*logkappa2);
  Type kappa2_pow4 = kappa2_pow2*kappa2_pow2;
  Type SigmaE2 = 1 / sqrt(4*pi*exp(2*logtauE2)*exp(2*logkappa2));
  Type SigmaO2 = 1 / sqrt(4*pi*exp(2*logtauO2)*exp(2*logkappa2));
  Type Range_raw2 = sqrt(8) / exp( logkappa2 );     // Range = approx. distance @ 10% correlation
  vector<Type> SigmaM( logSigmaM.size() );
  SigmaM = exp( logSigmaM );

  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));

  // Derived fields
  matrix<Type> Omega1_sc(n_s, n_c);
  array<Type> Epsilon1_sct(n_s, n_c, n_t);
  matrix<Type> Omega2_sc(n_s, n_c);
  array<Type> Epsilon2_sct(n_s, n_c, n_t);
  for(s=0; s<n_s; s++){
  for(c=0; c<n_c; c++){
    Omega1_sc(s,c) = Omegainput1_sc(s,c) / exp(logtauO1);
    Omega2_sc(s,c) = Omegainput2_sc(s,c) / exp(logtauO2);
    for(t=0;t<n_t;t++){
      Epsilon1_sct(s,c,t) = Epsiloninput1_sct(s,c,t) / exp(logtauE1);
      Epsilon2_sct(s,c,t) = Epsiloninput2_sct(s,c,t) / exp(logtauE2);
    }
  }}
  
  // Random field probability                                                                                                                              
  Eigen::SparseMatrix<Type> Q1;
  Eigen::SparseMatrix<Type> Q2;
  if( Options_vec(0)==0 ){
    Q1 = Q_spde(spde, exp(logkappa1));
    Q2 = Q_spde(spde, exp(logkappa2));
  }
  if( Options_vec(0)==1 ){
    Q1 = Q_spde(spde_aniso, exp(logkappa1), H);
    Q2 = Q_spde(spde_aniso, exp(logkappa2), H);
  }
  GMRF_t<Type> Tmp1 = GMRF(Q1); // SEPARABLE( AR1(rho_c1), GMRF(Q1) );
  GMRF_t<Type> Tmp2 = GMRF(Q2); // SEPARABLE( AR1(rho_c2), GMRF(Q2) );
  if(FieldConfig(0)==1) jnll_comp(0) = SEPARABLE( AR1(rho_c1), GMRF(Q1) )(Omegainput1_sc);
  if(FieldConfig(1)==1){
    for(t=0;t<n_t;t++){
      if(t==0) jnll_comp(1) += SEPARABLE( AR1(rho_c1), GMRF(Q1) )(Epsiloninput1_sct.col(t));
      if(t>=1) jnll_comp(1) += SEPARABLE( AR1(rho_c1), GMRF(Q1) )(Epsiloninput1_sct.col(t) - Epsilon_rho1*Epsiloninput1_sct.col(t-1));
    }
  }
  if(FieldConfig(2)==1) jnll_comp(2) = SEPARABLE( AR1(rho_c2), GMRF(Q2) )(Omegainput2_sc);
  if(FieldConfig(3)==1){
    for(t=0;t<n_t;t++){
      if(t==0) jnll_comp(3) += SEPARABLE( AR1(rho_c2), GMRF(Q2) )(Epsiloninput2_sct.col(t));
      if(t>=1) jnll_comp(3) += SEPARABLE( AR1(rho_c2), GMRF(Q2) )(Epsiloninput2_sct.col(t) - Epsilon_rho2*Epsiloninput2_sct.col(t-1));
    }
  }

  // Possible structure on betas
  if( Options_vec(2)!=0 ){
    for(c=0; c<n_c; c++){
    for(t=1; t<n_t; t++){
      jnll_comp(8) -= dnorm( beta1_ct(c,t), Beta_rho1*beta1_ct(c,t-1) + Beta_mean1, exp(logsigmaB1), true );
      jnll_comp(9) -= dnorm( beta2_ct(c,t), Beta_rho2*beta2_ct(c,t-1) + Beta_mean2, exp(logsigmaB2), true );
    }}
  }
  
  // Covariates
  vector<Type> eta1_x = X_xj * gamma1_j.matrix();
  vector<Type> zeta1_i = Q_ik * lambda1_k.matrix();
  vector<Type> eta2_x = X_xj * gamma2_j.matrix();
  vector<Type> zeta2_i = Q_ik * lambda2_k.matrix();
  array<Type> eta1_xct(n_x, n_c, n_t);
  array<Type> eta2_xct(n_x, n_c, n_t);
  eta1_xct.setZero();
  eta2_xct.setZero();
  for(int x=0; x<n_x; x++){
  for(int c=0; c<n_c; c++){
  for(int t=0; t<n_t; t++){
  for(int p=0; p<n_p; p++){
    eta1_xct(x,c,t) += gamma1_ctp(c,t,p) * X_xtp(x,t,p);
    eta2_xct(x,c,t) += gamma2_ctp(c,t,p) * X_xtp(x,t,p);
  }}}}

  // Derived quantities
  vector<Type> var_i(n_i);
  Type var_y;
  // Linear predictor (pre-link) for presence/absence component
  vector<Type> P1_i(n_i);   
  // Response predictor (post-link)
  // ObsModel = 0:4 or 11:12: probability ("phi") that data is greater than zero
  // ObsModel = 5 (ZINB):  phi = 1-ZeroInflation_prob -> Pr[D=0] = NB(0|mu,var)*phi + (1-phi) -> Pr[D>0] = phi - NB(0|mu,var)*phi 
  vector<Type> R1_i(n_i);   
  vector<Type> LogProb1_i(n_i);
  // Linear predictor (pre-link) for positive component
  vector<Type> P2_i(n_i);   
  // Response predictor (post-link)
  // ObsModel = 0:3, 11:12:  expected value of data, given that data is greater than zero -> E[D] = mu*phi
  // ObsModel = 4 (ZANB):  expected value ("mu") of neg-bin PRIOR to truncating Pr[D=0] -> E[D] = mu/(1-NB(0|mu,var))*phi  ALSO  Pr[D] = NB(D|mu,var)/(1-NB(0|mu,var))*phi
  // ObsModel = 5 (ZINB):  expected value of data for non-zero-inflation component -> E[D] = mu*phi
  vector<Type> R2_i(n_i);   
  vector<Type> LogProb2_i(n_i);

  // Likelihood contribution from observations
  for (int i=0;i<n_i;i++){
    // Presence-absence prediction
    P1_i(i) =  beta1_ct(c_i(i),t_i(i)) + Omega1_sc(s_i(i),c_i(i)) + Epsilon1_sct(s_i(i),c_i(i),t_i(i)) + eta1_x(s_i(i)) + eta1_xct(s_i(i),c_i(i),t_i(i)) + zeta1_i(i);
    R1_i(i) = invlogit( P1_i(i) ); 
    // Positive density prediction
    if( b_i(i)>0 | ObsModel(0)==5 | ObsModel(0)==6 ){    // 1e-500 causes overflow on laptop
      P2_i(i) =  beta2_ct(c_i(i),t_i(i)) + Omega2_sc(s_i(i),c_i(i)) + Epsilon2_sct(s_i(i),c_i(i),t_i(i)) + eta2_x(s_i(i)) + eta2_xct(s_i(i),c_i(i),t_i(i)) + zeta2_i(i);
      R2_i(i) = exp( P2_i(i) );
    }else{
      P2_i(i) = 0;
      R2_i(i) = 0;
    }                                               
    // Likelihood for models with continuous positive support
    if(ObsModel(0)==0 | ObsModel(0)==1 | ObsModel(0)==2 | ObsModel(0)==11 | ObsModel(0)==12){ 
      // Presence-absence likelihood
      if( b_i(i) > 0 ){
        LogProb1_i(i) = log( R1_i(i) );
      }else{
        LogProb1_i(i) = log( 1-R1_i(i) );
      }
      // Positive density likelihood -- models with continuous positive support
      if( b_i(i) > 0 ){    // 1e-500 causes overflow on laptop
        if(ObsModel(0)==0) LogProb2_i(i) = dnorm(b_i(i), R2_i(i)*a_i(i), SigmaM(0), true);
        if(ObsModel(0)==1) LogProb2_i(i) = dlnorm(b_i(i), P2_i(i)+log(a_i(i))-pow(SigmaM(0),2)/2, SigmaM(0), true); // log-space
        if(ObsModel(0)==2) LogProb2_i(i) = dgamma(b_i(i), 1/pow(SigmaM(0),2), R2_i(i)*a_i(i)*pow(SigmaM(0),2), true); // shape = 1/CV^2, scale = mean*CV^2
      }else{
        LogProb2_i(i) = 0;
      }
    }
    // Likelihood for models with discrete support 
    if(ObsModel(0)==4 | ObsModel(0)==5 | ObsModel(0)==6){
      var_i(i) = R2_i(i)*a_i(i)*(1.0+SigmaM(0)) + pow(R2_i(i)*a_i(i),2.0)*SigmaM(1);
      if(ObsModel(0)==5){
        if( b_i(i)==0 ){
          LogProb2_i(i) = log( (1-R1_i(i)) + dnbinom2(Type(0.0), R2_i(i)*a_i(i), var_i(i), false)*R1_i(i) ); //  Pr[X=0] = 1-phi + NB(X=0)*phi
        }else{
          LogProb2_i(i) = dnbinom2(b_i(i), R2_i(i)*a_i(i), var_i(i), true) + log(R1_i(i)); // Pr[X=x] = NB(X=x)*phi
        }
      }
      if(ObsModel(0)==6){
        LogProb2_i(i) = dCMP(b_i(i), R2_i(i)*a_i(i), exp(P1_i(i)), true, Options_vec(5));
      }
      LogProb1_i(i) = 0;
    }
  }

  // Predictive distribution -- ObsModel(4) isn't implemented (it had a bug previously)
  Type a_average = a_i.sum()/a_i.size();
  array<Type> P1_xct(n_x, n_c, n_t);
  array<Type> R1_xct(n_x, n_c, n_t);
  array<Type> P2_xct(n_x, n_c, n_t);
  array<Type> R2_xct(n_x, n_c, n_t);
  array<Type> D_xct(n_x, n_c, n_t);
  for(int c=0; c<n_c; c++){
  for(int t=0; t<n_t; t++){
  for(int x=0; x<n_x; x++){
    P1_xct(x,c,t) = beta1_ct(c,t) + Omega1_sc(x,c) + Epsilon1_sct(x,c,t) + eta1_x(x) + eta1_xct(x,c,t);
    R1_xct(x,c,t) = invlogit( P1_xct(x,c,t) );
    P2_xct(x,c,t) =  beta2_ct(c,t) + Omega2_sc(x,c) + Epsilon2_sct(x,c,t) + eta2_x(x) + eta2_xct(x,c,t);
    if(ObsModel(0)==0 | ObsModel(0)==1 | ObsModel(0)==2 | ObsModel(0)==4 | ObsModel(0)==5) R2_xct(x,c,t) = exp( P2_xct(x,c,t) );
    // Expected value for predictive distribution in a grid cell
    if(ObsModel(0)==0 | ObsModel(0)==1 | ObsModel(0)==2 | ObsModel(0)==5) D_xct(x,c,t) = R1_xct(x,c,t) * R2_xct(x,c,t);
  }}}
  
  // Calculate indices
  array<Type> Index_xctl(n_x, n_c, n_t, n_l);
  array<Type> Index_ctl(n_c, n_t, n_l);
  array<Type> ln_Index_ctl(n_c, n_t, n_l);
  Index_ctl.setZero();
  for(int c=0; c<n_c; c++){
  for(int t=0; t<n_t; t++){
  for(int l=0; l<n_l; l++){
    for(int x=0; x<n_x; x++){
      Index_xctl(x,c,t,l) = D_xct(x,c,t) * a_xl(x,l) / 1000;  // Convert from kg to metric tonnes
      Index_ctl(c,t,l) += Index_xctl(x,c,t,l);
    }
  }}}
  ln_Index_ctl = log( Index_ctl );

  // Calculate other derived summaries
  // Each is the weighted-average X_xl over polygons (x) with weights equal to abundance in each polygon and time (where abundance is from the first index)
  array<Type> mean_Z_ctm(n_c, n_t, n_m);
  mean_Z_ctm.setZero();
  int report_summary_TF = false;
  for(int c=0; c<n_c; c++){
  for(int t=0; t<n_t; t++){
  for(int m=0; m<n_m; m++){
    for(int x=0; x<n_x; x++){
      if( Z_xm(x,m)!=0 ){
        mean_Z_ctm(t,m) += Z_xm(x,m) * Index_xctl(x,c,t,0)/Index_ctl(c,t,0);
        report_summary_TF = true;
      }
    }
  }}}

  // Calculate average density, weighted.mean( x=Abundance/Area, w=Abundance )
  // Doesn't require Z_xm, because it only depends upon Index_tl
  if( Options(4)==1 ){
    array<Type> mean_D_ctl(n_c, n_t, n_l);
    array<Type> log_mean_D_ctl(n_c, n_t, n_l);
    mean_D_ctl.setZero();
    for(int c=0; c<n_c; c++){
    for(int t=0; t<n_t; t++){
    for(int l=0; l<n_l; l++){
      for(int x=0; x<n_x; x++){
        mean_D_ctl(c,t,l) += D_xct(x,c,t) * Index_xctl(x,c,t,l)/Index_ctl(c,t,l);
      }
    }}}
    REPORT( mean_D_ctl );
    ADREPORT( mean_D_ctl );
    log_mean_D_ctl = log( log_mean_D_ctl );
    ADREPORT( log_mean_D_ctl );

    // Calculate effective area = Index / average density
    array<Type> effective_area_ctl(n_c, n_t, n_l);
    array<Type> log_effective_area_ctl(n_c, n_t, n_l);
    effective_area_ctl = Index_ctl / (mean_D_ctl/1000);  // Correct for different units of Index and density
    log_effective_area_ctl = log( effective_area_ctl );
    REPORT( effective_area_ctl );
    ADREPORT( effective_area_ctl );
    ADREPORT( log_effective_area_ctl );
  }

  // Joint likelihood
  jnll_comp(10) = -1*LogProb1_i.sum();
  jnll_comp(11) = -1*LogProb2_i.sum();
  jnll = jnll_comp.sum();

  // Diagnostic output
  REPORT( P1_i );
  REPORT( P2_i );
  REPORT( R1_i );
  REPORT( R2_i );
  REPORT( P1_xct );
  REPORT( P2_xct );
  REPORT( var_i );
  REPORT( LogProb1_i );
  REPORT( LogProb2_i );
  REPORT( a_average );
  REPORT( eta1_x );
  REPORT( eta2_x );
  REPORT( eta1_xct );
  REPORT( eta2_xct );
  REPORT( zeta1_i );
  REPORT( zeta2_i );
  
  REPORT( SigmaE1 );
  REPORT( SigmaO1 );
  REPORT( SigmaE2 );
  REPORT( SigmaO2 );
  REPORT( SigmaM );
  REPORT( Index_ctl );
  REPORT( D_xct );
  REPORT( R1_xct );
  REPORT( R2_xct );
  REPORT( Index_xctl );
  REPORT( Omega1_sc );
  REPORT( Omega2_sc );
  REPORT( Epsilon1_sct );
  REPORT( Epsilon2_sct );
  REPORT( H );
  REPORT( Range_raw1 );
  REPORT( Range_raw2 );
  REPORT( beta1_ct );
  REPORT( beta2_ct );
  REPORT( jnll_comp );
  REPORT( jnll );

  ADREPORT( Range_raw1 );
  ADREPORT( Range_raw2 );
  ADREPORT( Index_ctl );
  ADREPORT( ln_Index_ctl);
  ADREPORT( SigmaE1 );
  ADREPORT( SigmaO1 );
  ADREPORT( SigmaE2 );
  ADREPORT( SigmaO2 );
  ADREPORT( SigmaM );

  return jnll;
  
}
