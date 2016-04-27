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

// Generate loadings matrix
template<class Type>
matrix<Type> loadings_matrix( vector<Type> L_val, int n_rows, int n_cols ){
  matrix<Type> L_rc(n_rows, n_cols);
  int Count = 0;
  for(int r=0; r<n_rows; r++){
  for(int c=0; c<n_cols; c++){
    if(r>=c){
      L_rc(r,c) = L_val(Count);
      Count++;
    }else{
      L_rc(r,c) = 0.0;
    }
  }}
  return L_rc;
}

// IN: eta1_vf; n_f_input; L1_z
// OUT: jnll_comp; eta1_vc
template<class Type>
matrix<Type> overdispersion_by_category_nll( int n_f, int n_f_input, int n_c, int n_v, matrix<Type> eta_vf, vector<Type> L_z, Type &jnll_pointer){
  using namespace density;
  matrix<Type> eta_vc(n_v, n_c);
  vector<Type> Tmp_c;
  // Turn off
  if(n_f_input<0){
    eta_vc.setZero();
  }
  // AR1 structure
  if( n_f_input==0 ){
    for(int v=0; v<n_v; v++){
      Tmp_c = eta_vc.row(v);
      jnll_pointer += SCALE( AR1(L_z(1)), exp(L_z(0)) )( Tmp_c );
    }
    eta_vc = eta_vf;
  }
  // Factor analysis structure
  if( n_f_input>0 ){
    // Assemble the loadings matrix
    matrix<Type> L_cf = loadings_matrix( L_z, n_c, n_f );
    // Multiply out overdispersion
    eta_vc = eta_vf * L_cf.transpose();
    // Probability of overdispersion
    for(int v=0; v<n_v; v++){
    for(int f=0; f<n_f; f++){
      jnll_pointer -= dnorm( eta_vf(v,f), Type(0.0), Type(1.0), true );
    }}
  }
  return eta_vc;
}

// Input: L_omega1_z, Q1, Omegainput1_sf, n_f, n_s, n_c, FieldConfig(0)
// Output: jnll_comp(0), Omega1_sc
template<class Type>                                                                                        //
matrix<Type> gmrf_by_category_nll( int n_f, int n_s, int n_c, Type logkappa, array<Type> gmrf_input_sf, vector<Type> L_z, Eigen::SparseMatrix<Type> Q, Type &jnll_pointer){
  using namespace density;
  matrix<Type> gmrf_sc(n_s, n_c);
  Type logtau;
  // Turn off
  if(n_f<0){
    gmrf_sc.setZero();
  }
  // AR1 structure
  if(n_f==0){
    logtau = L_z(0) - logkappa;  //
    jnll_pointer += SEPARABLE( AR1(L_z(1)), GMRF(Q) )(gmrf_input_sf);
    gmrf_sc = gmrf_input_sf / exp(logtau);                                // Rescaling from comp_index_v1d.cpp
  }
  // Factor analysis structure
  if(n_f>0){
    logtau = log( 1 / (exp(logkappa) * sqrt(4*M_PI)) );
    for( int f=0; f<n_f; f++ ) jnll_pointer += SCALE( GMRF(Q), exp(-logtau))(gmrf_input_sf.col(f));  // Rescaling from spatial_vam_v13.cpp
    matrix<Type> L_cf = loadings_matrix( L_z, n_c, n_f );
    gmrf_sc = gmrf_input_sf.matrix() * L_cf.transpose();
  }
  return gmrf_sc;
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
  DATA_INTEGER(n_v);          // Number of tows/vessels (i.e., levels for the factor explaining overdispersion)
  DATA_INTEGER(n_f_input);          // Number of factors for tows/vessels effect (allowable range: from 0 to n_c)
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
  // Slot 0: Distribution for positive catches
  // Slot 1: Link function for encounter probabilities
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
  DATA_IVECTOR(s_i);          // Station for each observation
  DATA_IVECTOR(t_i);          // Year for each observation
  DATA_IVECTOR(v_i);          // tows/vessels for each observation (level of factor representing overdispersion)
  DATA_VECTOR(PredTF_i);          // tows/vessels for each observation (level of factor representing overdispersion)
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
  PARAMETER_VECTOR(L1_z);          // Overdispersion parameters
  PARAMETER_VECTOR(L_omega1_z);
  PARAMETER_VECTOR(L_epsilon1_z);
  PARAMETER(logkappa1);
  PARAMETER(Beta_mean1);  // mean-reversion for beta1_t
  PARAMETER(logsigmaB1);  // SD of beta1_t (default: not included in objective function)
  PARAMETER(Beta_rho1);  // AR1 for positive catch Epsilon component, Default=0
  PARAMETER(Epsilon_rho1);  // AR1 for presence/absence Epsilon component, Default=0
  // -- Gaussian random fields
  PARAMETER_MATRIX(eta1_vf);
  PARAMETER_ARRAY(Omegainput1_sf);      // Expectation
  PARAMETER_ARRAY(Epsiloninput1_sft);   // Annual variation
  //  -- positive catch rates
  PARAMETER_MATRIX(beta2_ct);  // Year effect
  PARAMETER_VECTOR(gamma2_j);        // Covariate effect
  PARAMETER_ARRAY(gamma2_ctp);       // Dynamic covariate effect
  PARAMETER_VECTOR(lambda2_k);       // Catchability coefficients
  PARAMETER_VECTOR(L2_z);          // Overdispersion parameters
  PARAMETER_VECTOR(L_omega2_z);
  PARAMETER_VECTOR(L_epsilon2_z);
  PARAMETER(logkappa2);
  PARAMETER(Beta_mean2);  // mean-reversion for beta2_t
  PARAMETER(logsigmaB2);  // SD of beta2_t (default: not included in objective function)
  PARAMETER(Beta_rho2);  // AR1 for positive catch Epsilon component, Default=0
  PARAMETER(Epsilon_rho2);  // AR1 for positive catch Epsilon component, Default=0
  PARAMETER_ARRAY(logSigmaM);   // Slots: 0=mix1 CV, 1=prob-of-mix1, 2=
  // -- Gaussian random fields
  PARAMETER_MATRIX(eta2_vf);
  PARAMETER_ARRAY(Omegainput2_sf);      // Expectation
  PARAMETER_ARRAY(Epsiloninput2_sft);   // Annual variation

  // Indices -- i=Observation; j=Covariate; v=Vessel; t=Year; s=Stratum
  int i,j,v,t,s,c;
  
  // Objective function
  vector<Type> jnll_comp(13);
  // Slot 0 -- spatial, encounter
  // Slot 1 -- spatio-temporal, encounter
  // Slot 2 -- spatial, positive catch
  // Slot 3 -- spatio-temporal, positive catch
  // Slot 4 -- tow/vessel overdispersion, encounter
  // Slot 5 -- tow/vessel overdispersion, positive catch
  // Slot 8 -- penalty on beta, encounter
  // Slot 9 -- penalty on beta, positive catch
  // Slot 10 -- likelihood of data, encounter
  // Slot 11 -- likelihood of data, positive catch
  jnll_comp.setZero();
  Type jnll = 0;                
  
  // Derived parameters
  Type Range_raw1 = sqrt(8) / exp( logkappa1 );   // Range = approx. distance @ 10% correlation
  Type Range_raw2 = sqrt(8) / exp( logkappa2 );     // Range = approx. distance @ 10% correlation
  array<Type> SigmaM( n_c, 3 );
  SigmaM = exp( logSigmaM );

  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));

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
  // Probability of encounter
  array<Type> Omega1_sc(n_s, n_c);
  Omega1_sc = gmrf_by_category_nll(FieldConfig(0), n_s, n_c, logkappa1, Omegainput1_sf, L_omega1_z, Q1, jnll_comp(0));
  array<Type> Epsilon1_sct(n_s, n_c, n_t);
  for(t=0; t<n_t; t++){
    if(t==0) Epsilon1_sct.col(t) = gmrf_by_category_nll(FieldConfig(1), n_s, n_c, logkappa1, Epsiloninput1_sft.col(t), L_epsilon1_z, Q1, jnll_comp(1));
    if(t>=1){
      Epsilon1_sct.col(t) = gmrf_by_category_nll(FieldConfig(1), n_s, n_c, logkappa1, Epsiloninput1_sft.col(t)-Epsilon_rho1*Epsiloninput1_sft.col(t-1), L_epsilon1_z, Q1, jnll_comp(1));
      Epsilon1_sct.col(t) += Epsilon_rho1 * Epsilon1_sct.col(t-1);
    }
  }
  // Positive catch rate
  array<Type> Omega2_sc(n_s, n_c);
  Omega2_sc = gmrf_by_category_nll(FieldConfig(2), n_s, n_c, logkappa2, Omegainput2_sf, L_omega2_z, Q2, jnll_comp(2));
  array<Type> Epsilon2_sct(n_s, n_c, n_t);
  for(t=0;t<n_t;t++){
    if(t==0) Epsilon2_sct.col(t) = gmrf_by_category_nll(FieldConfig(3), n_s, n_c, logkappa2, Epsiloninput2_sft.col(t), L_epsilon2_z, Q2, jnll_comp(3));
    if(t>=1){
      Epsilon2_sct.col(t) = gmrf_by_category_nll(FieldConfig(3), n_s, n_c, logkappa2, Epsiloninput2_sft.col(t)-Epsilon_rho2*Epsiloninput2_sft.col(t-1), L_epsilon2_z, Q2, jnll_comp(3));
      Epsilon2_sct.col(t) += Epsilon_rho2 * Epsilon2_sct.col(t-1);
    }
  }

  ////// Probability of correlated overdispersion among bins
  int n_f;
  if( n_f_input==0 ){
    n_f = n_c;     // AR1 correlation structure on overdispersion
  }else{
    n_f = n_f_input; // Factor correlation structure on overdispersion, with rank n_f
  }
  // IN: eta1_vf; n_f_input; L1_z
  // OUT: jnll_comp; eta1_vc
  matrix<Type> eta1_vc(n_v, n_c);
  eta1_vc = overdispersion_by_category_nll( n_f, n_f_input, n_c, n_v, eta1_vf, L1_z, jnll_comp(4) );
  matrix<Type> eta2_vc(n_v, n_c);
  eta2_vc = overdispersion_by_category_nll( n_f, n_f_input, n_c, n_v, eta2_vf, L2_z, jnll_comp(5) );

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
    // Linear predictors
    P1_i(i) =  beta1_ct(c_i(i),t_i(i)) + Omega1_sc(s_i(i),c_i(i)) + Epsilon1_sct(s_i(i),c_i(i),t_i(i)) + eta1_x(s_i(i)) + eta1_xct(s_i(i),c_i(i),t_i(i)) + zeta1_i(i) + eta1_vc(v_i(i),c_i(i));
    P2_i(i) =  beta2_ct(c_i(i),t_i(i)) + Omega2_sc(s_i(i),c_i(i)) + Epsilon2_sct(s_i(i),c_i(i),t_i(i)) + eta2_x(s_i(i)) + eta2_xct(s_i(i),c_i(i),t_i(i)) + zeta2_i(i) + eta2_vc(v_i(i),c_i(i));
    // Responses
    if( ObsModel(1)==0 ){
      // P1_i: Logit-Probability of occurrence
      // P2_i: Log-Positive density prediction
      R1_i(i) = invlogit( P1_i(i) );
      R2_i(i) = exp( P2_i(i) );
    }
    if( ObsModel(1)==1 ){
      // P1_i: Log-Density
      // P2_i: Log-Variation in positive catch rates in excess of density
      R1_i(i) = Type(1.0) - exp( -1*SigmaM(c_i(i),2)*exp(P1_i(i)) );
      R2_i(i) = exp(P1_i(i)) / R1_i(i) * exp( P2_i(i) );
    }
    // Likelihood for models with continuous positive support
    if(ObsModel(0)==0 | ObsModel(0)==1 | ObsModel(0)==2){
      // Presence-absence likelihood
      if( b_i(i) > 0 ){
        LogProb1_i(i) = log( R1_i(i) );
      }else{
        LogProb1_i(i) = log( 1-R1_i(i) );
      }
      // Positive density likelihood -- models with continuous positive support
      if( b_i(i) > 0 ){    // 1e-500 causes overflow on laptop
        if(ObsModel(0)==0) LogProb2_i(i) = dnorm(b_i(i), R2_i(i)*a_i(i), SigmaM(c_i(i),0), true);
        if(ObsModel(0)==1) LogProb2_i(i) = dlnorm(b_i(i), log(R2_i(i)*a_i(i))-pow(SigmaM(c_i(i),0),2)/2, SigmaM(c_i(i),0), true); // log-space
        if(ObsModel(0)==2) LogProb2_i(i) = dgamma(b_i(i), 1/pow(SigmaM(c_i(i),0),2), R2_i(i)*a_i(i)*pow(SigmaM(c_i(i),0),2), true); // shape = 1/CV^2, scale = mean*CV^2
      }else{
        LogProb2_i(i) = 0;
      }
    }
    // Likelihood for models with discrete support 
    if(ObsModel(0)==4 | ObsModel(0)==5 | ObsModel(0)==6){
      var_i(i) = R2_i(i)*a_i(i)*(1.0+SigmaM(c_i(i),0)) + pow(R2_i(i)*a_i(i),2.0)*SigmaM(c_i(i),1);
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
    P2_xct(x,c,t) =  beta2_ct(c,t) + Omega2_sc(x,c) + Epsilon2_sct(x,c,t) + eta2_x(x) + eta2_xct(x,c,t);
    if( ObsModel(1)==0 ){
      R1_xct(x,c,t) = invlogit( P1_xct(x,c,t) );
      R2_xct(x,c,t) = exp( P2_xct(x,c,t) );
    }
    if( ObsModel(1)==1 ){
      R1_xct(x,c,t) = Type(1.0) - exp( -SigmaM(c,2)*exp(P1_xct(x,c,t)) );
      R2_xct(x,c,t) = exp(P1_xct(x,c,t)) / R1_xct(x,c,t) * exp( P2_xct(x,c,t) );
    }
    // Expected value for predictive distribution in a grid cell
    D_xct(x,c,t) = R1_xct(x,c,t) * R2_xct(x,c,t);
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
  jnll_comp(10) = -1 * (LogProb1_i*(Type(1.0)-PredTF_i)).sum();
  jnll_comp(11) = -1 * (LogProb2_i*(Type(1.0)-PredTF_i)).sum();
  jnll = jnll_comp.sum();
  Type pred_jnll = -1 * ( LogProb1_i*PredTF_i + LogProb2_i*PredTF_i ).sum();
  REPORT( pred_jnll );

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
  REPORT( eta1_vc );
  REPORT( eta2_vc );
  REPORT( eta1_vf );
  REPORT( eta2_vf );
  REPORT( zeta1_i );
  REPORT( zeta2_i );
  
  REPORT( SigmaM );
  REPORT( Index_ctl );
  REPORT( D_xct );
  REPORT( R1_xct );
  REPORT( R2_xct );
  REPORT( Index_xctl );
  REPORT( Omega1_sc );
  REPORT( Omega2_sc );
  REPORT( Omegainput1_sf );
  REPORT( Omegainput2_sf );
  REPORT( Epsilon1_sct );
  REPORT( Epsilon2_sct );
  REPORT( Epsiloninput1_sft );
  REPORT( Epsiloninput2_sft );
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
  ADREPORT( SigmaM );

  return jnll;
  
}
