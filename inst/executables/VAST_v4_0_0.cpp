#include <TMB.hpp>
#include <Eigen/Eigenvalues>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Posfun
template<class Type>
Type posfun(Type x, Type lowerlimit, Type &pen){
  pen += CppAD::CondExpLt(x,lowerlimit,Type(0.01)*pow(x-lowerlimit,2),Type(0));
  return CppAD::CondExpGe(x,lowerlimit,x,lowerlimit/(Type(2)-x/lowerlimit));
}

// Variance
template<class Type>
Type var( array<Type> vec ){
  Type vec_mod = vec - (vec.sum()/vec.size());
  Type res = pow(vec_mod, 2).sum() / vec.size();
  return res;
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

// IN: eta1_vf; L1_z
// OUT: jnll_comp; eta1_vc
template<class Type>
matrix<Type> overdispersion_by_category_nll( int n_f, int n_v, int n_c, matrix<Type> eta_vf, vector<Type> L_z, Type &jnll_pointer){
  using namespace density;
  matrix<Type> eta_vc(n_v, n_c);
  vector<Type> Tmp_c;
  // Turn off
  if(n_f<0){
    eta_vc.setZero();
  }
  // AR1 structure
  if( n_f==0 ){
    for(int v=0; v<n_v; v++){
      Tmp_c = eta_vc.row(v);
      jnll_pointer += SCALE( AR1(L_z(1)), exp(L_z(0)) )( Tmp_c );
    }
    eta_vc = eta_vf;
  }
  // Factor analysis structure
  if( n_f>0 ){
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

template<class Type>                                                                                        //
matrix<Type> convert_upper_cov_to_cor( matrix<Type> cov ){
  int nrow = cov.row(0).size();
  for( int i=0; i<nrow; i++){
  for( int j=i+1; j<nrow; j++){
    cov(i,j) = cov(i,j) / pow(cov(i,i),0.5) / pow(cov(j,j),0.5);
  }}
  return cov;
}

// Input: L_omega1_z, Q1, Omegainput1_sf, n_f, n_s, n_c, FieldConfig(0)
// Output: jnll_comp(0), Omega1_sc
template<class Type>                                                                                        //
matrix<Type> gmrf_by_category_nll( int n_f, int method, int n_s, int n_c, Type logkappa, array<Type> gmrf_input_sf, array<Type> gmrf_mean_sf, vector<Type> L_z, density::GMRF_t<Type> gmrf_Q, Type &jnll_pointer){
  using namespace density;
  matrix<Type> gmrf_sc(n_s, n_c);
  Type logtau;
  // IID
  if(n_f == -2){
    for( int c=0; c<n_c; c++ ){
      if(method==0) logtau = log( 1 / (exp(logkappa) * sqrt(4*M_PI)) );
      if(method==1) logtau = log( 1 / sqrt(1-exp(logkappa*2)) );
      jnll_pointer += gmrf_Q(gmrf_input_sf.col(c) - gmrf_mean_sf.col(c));
      gmrf_sc.col(c) = gmrf_input_sf.col(c) / exp(logtau) * L_z(c);                                // Rescaling from comp_index_v1d.cpp
    }
  }
  // Turn off
  if(n_f == -1){
    gmrf_sc.setZero();
  }
  // AR1 structure
  if(n_f==0){
    logtau = L_z(0) - logkappa;  //
    jnll_pointer += SEPARABLE( AR1(L_z(1)), gmrf_Q )(gmrf_input_sf - gmrf_mean_sf);
    gmrf_sc = gmrf_input_sf / exp(logtau);                                // Rescaling from comp_index_v1d.cpp
  }
  // Factor analysis structure
  if(n_f>0){
    if(method==0) logtau = log( 1 / (exp(logkappa) * sqrt(4*M_PI)) );
    if(method==1) logtau = log( 1 / sqrt(1-exp(logkappa*2)) );
    for( int f=0; f<n_f; f++ ) jnll_pointer += gmrf_Q(gmrf_input_sf.col(f) - gmrf_mean_sf.col(f));  // Rescaling from spatial_vam_v13.cpp
    matrix<Type> L_cf = loadings_matrix( L_z, n_c, n_f );
    gmrf_sc = (gmrf_input_sf.matrix() * L_cf.transpose()) / exp(logtau);
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

// compound Poisson-Gamma ("Tweedie") distribution
template<class Type>
Type dPoisGam( Type x, Type shape, Type scale, Type intensity, vector<Type> &diag_z, int maxsum=50, int minsum=1, int give_log=0 ){
  // Maximum integration constant to prevent numerical overflow, but capped at value for maxsum to prevent numerical underflow when subtracting by a higher limit than is seen in the sequence
  Type max_log_wJ, z1, maxJ_bounded;
  if( x==0 ){
    diag_z(0) = 1;
    max_log_wJ = 0;
    diag_z(1) = 0;
  }else{
    z1 = log(intensity) + shape*log(x/scale) - shape*log(shape) + 1;
    diag_z(0) = exp( (z1 - 1) / (1 + shape) );
    maxJ_bounded = CppAD::CondExpGe(diag_z(0), Type(maxsum), Type(maxsum), diag_z(0));
    max_log_wJ = maxJ_bounded*log(intensity) + (maxJ_bounded*shape)*log(x/scale) - lgamma(maxJ_bounded+1) - lgamma(maxJ_bounded*shape);
    diag_z(1) = diag_z(0)*log(intensity) + (diag_z(0)*shape)*log(x/scale) - lgamma(diag_z(0)+1) - lgamma(diag_z(0)*shape);
  }
  // Integration constant
  Type W = 0;
  Type log_w_j, pos_penalty;
  for( int j=minsum; j<=maxsum; j++ ){
    Type j2 = j;
    //W += pow(intensity,j) * pow(x/scale, j2*shape) / exp(lgamma(j2+1)) / exp(lgamma(j2*shape)) / exp(max_log_w_j);
    log_w_j = j2*log(intensity) + (j2*shape)*log(x/scale) - lgamma(j2+1) - lgamma(j2*shape);
    //W += exp( posfun(log_w_j, Type(-30), pos_penalty) );
    W += exp( log_w_j - max_log_wJ );
    if(j==minsum) diag_z(2) = log_w_j;
    if(j==maxsum) diag_z(3) = log_w_j;
  }
  // Loglikelihood calculation
  Type loglike = 0;
  if( x==0 ){
    loglike = -intensity;
  }else{
    loglike = -x/scale - intensity - log(x) + log(W) + max_log_wJ;
  }
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
  DATA_INTEGER(n_t);         // Number of time-indices
  DATA_INTEGER(n_c);         // Number of categories (e.g., length bins)
  DATA_INTEGER(n_e);         // Number of error distributions
  DATA_INTEGER(n_j);         // Number of static covariates
  DATA_INTEGER(n_p);         // Number of dynamic covariates
  DATA_INTEGER(n_k);          // Number of catchability variables
  DATA_INTEGER(n_v);          // Number of tows/vessels (i.e., levels for the factor explaining overdispersion)
  DATA_INTEGER(n_l);         // Number of indices to post-process
  DATA_INTEGER(n_m);         // Number of range metrics to use (probably 2 for Eastings-Northings)

  // Config
  DATA_IVECTOR( Options_vec );
  // Slot 0 -- Aniso: 0=No, 1=Yes
  // Slot 1 -- DEPRECATED
  // Slot 2 -- AR1 on betas (year intercepts) to deal with missing years: 0=No, 1=Yes
  // Slot 3 -- DEPRECATED
  // Slot 4 -- DEPRECATED
  // Slot 5 -- Upper limit constant of integration calculation for infinite-series density functions (Conway-Maxwell-Poisson and Tweedie)
  // Slot 6 -- Breakpoint in CMP density function
  // Slot 7 -- Whether to use SPDE or 2D-AR1 hyper-distribution for spatial process: 0=SPDE; 1=2D-AR1
  DATA_IVECTOR(FieldConfig);  // Input settings (vector, length 4)
  DATA_IVECTOR(OverdispersionConfig);          // Input settings (vector, length 2)
  DATA_IMATRIX(ObsModel_ez);    // Observation model
  // Column 0: Probability distribution for data for each level of e_i
  // Column 1: Link function for linear predictors for each level of c_i
  // NOTE:  nlevels(c_i) must be <= nlevels(e_i)
  DATA_IVECTOR(Options);    // Reporting options
  // Slot 0: Calculate SD for Index_xctl
  // Slot 1: Calculate SD for log(Index_xctl)
  // Slot 2: Calculate mean_Z_ctm (i.e., center-of-gravity)
  // Slot 3: DEPRECATED
  // Slot 4: Calculate mean_D_tl and effective_area_tl
  // Slot 5: Calculate standard errors for Covariance and Correlation among categories using factor-analysis parameterization
  // Slot 6: Calculate synchrony for different periods specified via yearbounds_zz
  // Slot 7: Calculate coherence and variance for Epsilon1_sct and Epsilon2_sct
  DATA_IMATRIX(yearbounds_zz);
  // Two columns, and 1+ rows, specifying first and last t for each period used in calculating synchrony

  // Data vectors
  DATA_VECTOR(b_i);       	// Response (biomass) for each observation
  DATA_VECTOR(a_i);       	// Area swept for each observation (km^2)
  DATA_IMATRIX(c_iz);         // Category for each observation
  DATA_IVECTOR(e_i);         // Error distribution for each observation
  DATA_IVECTOR(s_i);          // Station for each observation
  DATA_IMATRIX(t_iz);          // Time-indices (year, season, etc.) for each observation
  DATA_IVECTOR(v_i);          // tows/vessels for each observation (level of factor representing overdispersion)
  DATA_VECTOR(PredTF_i);          // vector indicating whether an observatino is predictive (1=used for model evaluation) or fitted (0=used for parameter estimation)
  DATA_MATRIX(a_xl);		     // Area for each "real" stratum(km^2) in each stratum
  DATA_MATRIX(X_xj);		    // Covariate design matrix (strata x covariate)
  DATA_ARRAY(X_xtp);		    // Covariate design matrix (strata x covariate)
  DATA_MATRIX(Q_ik);        // Catchability matrix (observations x variable)
  DATA_IMATRIX(t_yz);        // Matrix for time-indices of calculating outputs (abundance index and "derived-quantity")
  DATA_MATRIX(Z_xm);        // Derived quantity matrix

  // SPDE objects
  DATA_STRUCT(spde,spde_t);
  
  // Aniso objects
  DATA_STRUCT(spde_aniso,spde_aniso_t);

  // Sparse matrices for precision matrix of 2D AR1 process
  // Q = M0*(1+rho^2)^2 + M1*(1+rho^2)*(-rho) + M2*rho^2
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Parameters 
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters

  //  -- presence/absence fixed effects
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
  PARAMETER_VECTOR(log_sigmaratio1_z);  // Ratio of variance for columns of t_iz

  // -- presence/absence random effects
  PARAMETER_MATRIX(eta1_vf);
  PARAMETER_ARRAY(Omegainput1_sf);      // Expectation
  PARAMETER_ARRAY(Epsiloninput1_sft);   // Annual variation

  //  -- positive catch rates fixed effects
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
  PARAMETER_VECTOR(log_sigmaratio2_z);  // Ratio of variance for columns of t_iz

  // Error distribution parameters
  PARAMETER_ARRAY(logSigmaM);
  // Columns: 0=CV, 1=[usually not used], 2=[usually not used]
  // Rows:  Each level of e_i and/or c_i
  // SigmaM[,0] indexed by e_i, e.g., SigmaM(e_i(i),0)
  // SigmaM[,1] and SigmaM[,2] indexed by c_i, e.g., SigmaM(c_i(i),2)

  // -- positive catch rates random effects
  PARAMETER_VECTOR(delta_i);
  PARAMETER_MATRIX(eta2_vf);
  PARAMETER_ARRAY(Omegainput2_sf);      // Expectation
  PARAMETER_ARRAY(Epsiloninput2_sft);   // Annual variation

  // Indices -- i=Observation; j=Covariate; v=Vessel; t=Year; s=Stratum
  int i,j,v,t,s,c,y,z;
  
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
  // Slot 12 -- Likelihood of Lognormal-Poisson overdispersion delta_i
  jnll_comp.setZero();
  Type jnll = 0;                

  // Derived parameters
  Type Range_raw1, Range_raw2;
  if( Options_vec(7)==0 ){
    Range_raw1 = sqrt(8) / exp( logkappa1 );   // Range = approx. distance @ 10% correlation
    Range_raw2 = sqrt(8) / exp( logkappa2 );     // Range = approx. distance @ 10% correlation
  }
  if( Options_vec(7)==1 ){
    Range_raw1 = log(0.1) / logkappa1;   // Range = approx. distance @ 10% correlation
    Range_raw2 = log(0.1) / logkappa2;     // Range = approx. distance @ 10% correlation
  }
  array<Type> SigmaM( n_e, 3 );
  SigmaM = exp( logSigmaM );

  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));

  ////////////////////////
  // Calculate joint likelihood
  ////////////////////////

  // Random field probability
  Eigen::SparseMatrix<Type> Q1;
  Eigen::SparseMatrix<Type> Q2;
  GMRF_t<Type> gmrf_Q;
  if( Options_vec(7)==0 & Options_vec(0)==0 ){
    Q1 = Q_spde(spde, exp(logkappa1));
    Q2 = Q_spde(spde, exp(logkappa2));
  }
  if( Options_vec(7)==0 & Options_vec(0)==1 ){
    Q1 = Q_spde(spde_aniso, exp(logkappa1), H);
    Q2 = Q_spde(spde_aniso, exp(logkappa2), H);
  }
  if( Options_vec(7)==1 ){
    Q1 = M0*pow(1+exp(logkappa1*2),2) + M1*(1+exp(logkappa1*2))*(-exp(logkappa1)) + M2*exp(logkappa1*2);
    Q2 = M0*pow(1+exp(logkappa2*2),2) + M1*(1+exp(logkappa2*2))*(-exp(logkappa2)) + M2*exp(logkappa2*2);
  }
  // Probability of encounter
  gmrf_Q = GMRF(Q1);
  // Omega1
  int n_f;
  n_f = Omegainput1_sf.cols();
  array<Type> Omegamean1_sf(n_s, n_f);
  Omegamean1_sf.setZero();
  // Epsilon1
  n_f = Epsiloninput1_sft.col(0).cols();
  array<Type> Epsilonmean1_sf(n_s, n_f);
  array<Type> Omega1_sc(n_s, n_c);
  Omega1_sc = gmrf_by_category_nll(FieldConfig(0), Options_vec(7), n_s, n_c, logkappa1, Omegainput1_sf, Omegamean1_sf, L_omega1_z, gmrf_Q, jnll_comp(0));
  array<Type> Epsilon1_sct(n_s, n_c, n_t);
  for(t=0; t<n_t; t++){
    if(t==0){
      Epsilonmean1_sf.setZero();
      Epsilon1_sct.col(t) = gmrf_by_category_nll(FieldConfig(1), Options_vec(7), n_s, n_c, logkappa1, Epsiloninput1_sft.col(t), Epsilonmean1_sf, L_epsilon1_z, gmrf_Q, jnll_comp(1));
    }
    if(t>=1){
      Epsilonmean1_sf = Epsilon_rho1 * Epsiloninput1_sft.col(t-1);
      Epsilon1_sct.col(t) = gmrf_by_category_nll(FieldConfig(1), Options_vec(7), n_s, n_c, logkappa1, Epsiloninput1_sft.col(t), Epsilonmean1_sf, L_epsilon1_z, gmrf_Q, jnll_comp(1));
    }
  }
  // Positive catch rate
  gmrf_Q = GMRF(Q2);
  // Omega2
  n_f = Omegainput2_sf.cols();
  array<Type> Omegamean2_sf(n_s, n_f);
  Omegamean2_sf.setZero();
  // Epsilon2
  n_f = Epsiloninput2_sft.col(0).cols();
  array<Type> Epsilonmean2_sf(n_s, n_f);
  array<Type> Omega2_sc(n_s, n_c);
  Omega2_sc = gmrf_by_category_nll(FieldConfig(2), Options_vec(7), n_s, n_c, logkappa2, Omegainput2_sf, Omegamean2_sf, L_omega2_z, gmrf_Q, jnll_comp(2));
  array<Type> Epsilon2_sct(n_s, n_c, n_t);
  for(t=0; t<n_t; t++){
    if(t==0){
      Epsilonmean2_sf.setZero();
      Epsilon2_sct.col(t) = gmrf_by_category_nll(FieldConfig(3), Options_vec(7), n_s, n_c, logkappa2, Epsiloninput2_sft.col(t), Epsilonmean2_sf, L_epsilon2_z, gmrf_Q, jnll_comp(3));
    }
    if(t>=1){
      Epsilonmean2_sf = Epsilon_rho2 * Epsiloninput2_sft.col(t-1);
      Epsilon2_sct.col(t) = gmrf_by_category_nll(FieldConfig(3), Options_vec(7), n_s, n_c, logkappa2, Epsiloninput2_sft.col(t), Epsilonmean2_sf, L_epsilon2_z, gmrf_Q, jnll_comp(3));
    }
  }

  ////// Probability of correlated overdispersion among bins
  // IN: eta1_vf; n_f; L1_z
  // OUT: jnll_comp; eta1_vc
  matrix<Type> eta1_vc(n_v, n_c);
  eta1_vc = overdispersion_by_category_nll( OverdispersionConfig(0), n_v, n_c, eta1_vf, L1_z, jnll_comp(4) );
  matrix<Type> eta2_vc(n_v, n_c);
  eta2_vc = overdispersion_by_category_nll( OverdispersionConfig(1), n_v, n_c, eta2_vf, L2_z, jnll_comp(5) );

  // Possible structure on betas
  if( Options_vec(2)!=0 ){
    for(c=0; c<n_c; c++){
    for(t=1; t<n_t; t++){
      jnll_comp(8) -= dnorm( beta1_ct(c,t), Beta_rho1*beta1_ct(c,t-1) + Beta_mean1, exp(logsigmaB1), true );
      jnll_comp(9) -= dnorm( beta2_ct(c,t), Beta_rho2*beta2_ct(c,t-1) + Beta_mean2, exp(logsigmaB2), true );
    }}
  }
  
  // Penalty on lognormal-Poisson overdispesrion delta_i
  for(i=0; i<delta_i.size(); i++){
    if( ObsModel_ez(e_i(i),0)==11 ){
      jnll_comp(12) -= dnorm( delta_i(i), Type(0.0), Type(1.0), true );
    }
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
  Type tmp_calc1;
  Type tmp_calc2;
  // Linear predictor (pre-link) for presence/absence component
  matrix<Type> P1_iz(n_i,c_iz.row(0).size());
  // Response predictor (post-link)
  // ObsModel_ez(e,0) = 0:4 or 11:12: probability ("phi") that data is greater than zero
  // ObsModel_ez(e,0) = 5 (ZINB):  phi = 1-ZeroInflation_prob -> Pr[D=0] = NB(0|mu,var)*phi + (1-phi) -> Pr[D>0] = phi - NB(0|mu,var)*phi
  vector<Type> R1_i(n_i);   
  vector<Type> log_one_minus_R1_i(n_i);
  vector<Type> LogProb1_i(n_i);
  // Linear predictor (pre-link) for positive component
  matrix<Type> P2_iz(n_i,c_iz.row(0).size());
  // Response predictor (post-link)
  // ObsModel_ez(e,0) = 0:3, 11:12:  expected value of data, given that data is greater than zero -> E[D] = mu*phi
  // ObsModel_ez(e,0) = 4 (ZANB):  expected value ("mu") of neg-bin PRIOR to truncating Pr[D=0] -> E[D] = mu/(1-NB(0|mu,var))*phi  ALSO  Pr[D] = NB(D|mu,var)/(1-NB(0|mu,var))*phi
  // ObsModel_ez(e,0) = 5 (ZINB):  expected value of data for non-zero-inflation component -> E[D] = mu*phi
  vector<Type> R2_i(n_i);   
  vector<Type> LogProb2_i(n_i);
  vector<Type> maxJ_i(n_i);
  vector<Type> diag_z(4);
  matrix<Type> diag_iz(n_i,4);
  diag_iz.setZero();  // Used to track diagnostics for Tweedie distribution (columns: 0=maxJ; 1=maxW; 2=lowerW; 3=upperW)
  P1_iz.setZero();
  P2_iz.setZero();

  // Likelihood contribution from observations
  for(int i=0; i<n_i; i++){
    if( !isNA(b_i(i)) ){
      // Linear predictors
      for( int zc=0; zc<c_iz.row(0).size(); zc++ ){
        if( c_iz(i,zc)>=0 & c_iz(i,zc)<n_c ){
        //if( !isNA(c_iz(i,zc)) ){
          P1_iz(i,zc) = Omega1_sc(s_i(i),c_iz(i,zc)) + eta1_x(s_i(i)) + zeta1_i(i) + eta1_vc(v_i(i),c_iz(i,zc));
          P2_iz(i,zc) = Omega2_sc(s_i(i),c_iz(i,zc)) + eta2_x(s_i(i)) + zeta2_i(i) + eta2_vc(v_i(i),c_iz(i,zc));
          for( int zt=0; zt<t_iz.row(0).size(); zt++ ){
            if( t_iz(i,zt)>=0 & t_iz(i,zt)<n_t ){  // isNA doesn't seem to work for IMATRIX type
              P1_iz(i,zc) += beta1_ct(c_iz(i,zc),t_iz(i,zt)) + Epsilon1_sct(s_i(i),c_iz(i,zc),t_iz(i,zt))*exp(log_sigmaratio1_z(zt)) + eta1_xct(s_i(i),c_iz(i,zc),t_iz(i,zt));
              P2_iz(i,zc) += beta2_ct(c_iz(i,zc),t_iz(i,zt)) + Epsilon2_sct(s_i(i),c_iz(i,zc),t_iz(i,zt))*exp(log_sigmaratio2_z(zt)) + eta2_xct(s_i(i),c_iz(i,zc),t_iz(i,zt));
            }
          }
        }
      }
      // Responses
      if( ObsModel_ez(c_iz(i,0),1)==0 | ObsModel_ez(c_iz(i,0),1)==3 ){
        // Log and logit-link, where area-swept only affects positive catch rate exp(P2_i(i))
        // P1_i: Logit-Probability of occurrence;  R1_i:  Probability of occurrence
        // P2_i: Log-Positive density prediction;  R2_i:  Positive density prediction
        R1_i(i) = invlogit( P1_iz(i,0) );
        R2_i(i) = a_i(i) * exp( P2_iz(i,0) );
      }
      if( ObsModel_ez(c_iz(i,0),1)==1 ){
        // Poisson-process link, where area-swept affects numbers density exp(P1_i(i))
        // P1_i: Log-numbers density;  R1_i:  Probability of occurrence
        // P2_i: Log-average weight;  R2_i:  Positive density prediction
        tmp_calc1 = 0;
        tmp_calc2 = 0;
        for( int zc=0; zc<c_iz.row(0).size(); zc++ ){
          if( c_iz(i,zc)>=0 & c_iz(i,zc)<n_c ){
          //if( !isNA(c_iz(i,zc)) ){
            tmp_calc1 += exp(P1_iz(i,zc));
            tmp_calc2 += exp(P1_iz(i,zc)) * exp(P2_iz(i,zc));
          }
        }
        R1_i(i) = Type(1.0) - exp( -1*a_i(i)*tmp_calc1 );
        R2_i(i) = a_i(i) * tmp_calc2 / R1_i(i);
        // log_one_minus_R1_i is useful to prevent numerical underflow e.g., for 1 - exp(-40)
        log_one_minus_R1_i(i) = -1*a_i(i)*tmp_calc1;
      }
      if( ObsModel_ez(c_iz(i,0),1)==2 ){
        // Tweedie link, where area-swept affects numbers density exp(P1_i(i))
        // P1_i: Log-numbers density;  R1_i:  Expected numbers
        // P2_i: Log-average weight;  R2_i:  Expected average weight
        R1_i(i) = a_i(i) * exp( P1_iz(i,0) );
        R2_i(i) = exp( P2_iz(i,0) );
      }
      // Likelihood for delta-models with continuous positive support
      if(ObsModel_ez(e_i(i),0)==0 | ObsModel_ez(e_i(i),0)==1 | ObsModel_ez(e_i(i),0)==2){
        // Presence-absence likelihood
        if( b_i(i) > 0 ){
          LogProb1_i(i) = log( R1_i(i) );
        }else{
          if( ObsModel_ez(e_i(i),1)==1 ){
            LogProb1_i(i) = log_one_minus_R1_i(i);
          }else{
            LogProb1_i(i) = log( 1-R1_i(i) );
          }
        }
        // Positive density likelihood -- models with continuous positive support
        if( b_i(i) > 0 ){    // 1e-500 causes overflow on laptop
          if(ObsModel_ez(e_i(i),0)==0) LogProb2_i(i) = dnorm(b_i(i), R2_i(i), SigmaM(e_i(i),0), true);
          if(ObsModel_ez(e_i(i),0)==1) LogProb2_i(i) = dlnorm(b_i(i), log(R2_i(i))-pow(SigmaM(e_i(i),0),2)/2, SigmaM(e_i(i),0), true); // log-space
          if(ObsModel_ez(e_i(i),0)==2) LogProb2_i(i) = dgamma(b_i(i), 1/pow(SigmaM(e_i(i),0),2), R2_i(i)*pow(SigmaM(e_i(i),0),2), true); // shape = 1/CV^2, scale = mean*CV^2
        }else{
          LogProb2_i(i) = 0;
        }
      }
      // Likelihood for Tweedie model with continuous positive support
      if(ObsModel_ez(e_i(i),0)==8){
        LogProb1_i(i) = 0;
        //dPoisGam( Type x, Type shape, Type scale, Type intensity, Type &max_log_w_j, int maxsum=50, int minsum=1, int give_log=0 )
        LogProb2_i(i) = dPoisGam( b_i(i), SigmaM(e_i(i),0), R2_i(i), R1_i(i), diag_z, Options_vec(5), Options_vec(6), true );
        diag_iz.row(i) = diag_z;
      }
      if(ObsModel_ez(e_i(i),0)==10){
        // Packaged code
        LogProb1_i(i) = 0;
        // dtweedie( Type y, Type mu, Type phi, Type p, int give_log=0 )
        // R1*R2 = mean
        LogProb2_i(i) = dtweedie( b_i(i), R1_i(i)*R2_i(i), R1_i(i), invlogit(SigmaM(e_i(i),0))+1.0, true );
      }
      // Likelihood for models with discrete support
      if(ObsModel_ez(e_i(i),0)==4 | ObsModel_ez(e_i(i),0)==5 | ObsModel_ez(e_i(i),0)==6 | ObsModel_ez(e_i(i),0)==7 | ObsModel_ez(e_i(i),0)==9 | ObsModel_ez(e_i(i),0)==11){
        if(ObsModel_ez(e_i(i),0)==5){
          // Zero-inflated negative binomial (not numerically stable!)
          var_i(i) = R2_i(i)*(1.0+SigmaM(e_i(i),0)) + pow(R2_i(i),2.0)*SigmaM(c_iz(i,0),1);
          if( b_i(i)==0 ){
            //LogProb2_i(i) = log( (1-R1_i(i)) + dnbinom2(Type(0.0), R2_i(i), var_i(i), false)*R1_i(i) ); //  Pr[X=0] = 1-phi + NB(X=0)*phi
            LogProb2_i(i) = logspace_add( log(1-R1_i(i)), dnbinom2(Type(0.0),R2_i(i),var_i(i),true)+log(R1_i(i)) ); //  Pr[X=0] = 1-phi + NB(X=0)*phi
          }else{
            LogProb2_i(i) = dnbinom2(b_i(i), R2_i(i), var_i(i), true) + log(R1_i(i)); // Pr[X=x] = NB(X=x)*phi
          }
        }
        if(ObsModel_ez(e_i(i),0)==6){
          // Conway-Maxwell-Poisson
          LogProb2_i(i) = dCMP(b_i(i), R2_i(i), exp(P1_iz(i,0)), true, Options_vec(5));
        }
        if(ObsModel_ez(e_i(i),0)==7){
          // Zero-inflated Poisson
          if( b_i(i)==0 ){
            //LogProb2_i(i) = log( (1-R1_i(i)) + dpois(Type(0.0), R2_i(i), false)*R1_i(i) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
            LogProb2_i(i) = logspace_add( log(1-R1_i(i)), dpois(Type(0.0),R2_i(i),true)+log(R1_i(i)) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
          }else{
            LogProb2_i(i) = dpois(b_i(i), R2_i(i), true) + log(R1_i(i)); // Pr[X=x] = Pois(X=x)*phi
          }
        }
        if(ObsModel_ez(e_i(i),0)==9){
          // Binned Poisson (for REEF data: 0=none; 1=1; 2=2-10; 3=>11)
          /// Doesn't appear stable given spatial or spatio-temporal variation
          vector<Type> logdBinPois(4);
          logdBinPois(0) = logspace_add( log(1-R1_i(i)), dpois(Type(0.0), R2_i(i), true) + log(R1_i(i)) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
          logdBinPois(1) = dpois(Type(1.0), R2_i(i), true) + log(R1_i(i));                                 //  Pr[X | X>0] = Pois(X)*phi
          logdBinPois(2) = dpois(Type(2.0), R2_i(i), true) + log(R1_i(i));                                 // SUM_J( Pr[X|X>0] ) = phi * SUM_J( Pois(J) )
          for(int j=3; j<=10; j++){
            logdBinPois(2) += logspace_add( logdBinPois(2), dpois(Type(j), R2_i(i), true) + log(R1_i(i)) );
          }
          logdBinPois(3) = logspace_sub( log(Type(1.0)), logdBinPois(0) );
          logdBinPois(3) = logspace_sub( logdBinPois(3), logdBinPois(1) );
          logdBinPois(3) = logspace_sub( logdBinPois(3), logdBinPois(2) );
          if( b_i(i)==0 ) LogProb2_i(i) = logdBinPois(0);
          if( b_i(i)==1 ) LogProb2_i(i) = logdBinPois(1);
          if( b_i(i)==2 ) LogProb2_i(i) = logdBinPois(2);
          if( b_i(i)==3 ) LogProb2_i(i) = logdBinPois(3);
        }
        if(ObsModel_ez(e_i(i),0)==11){
          // Zero-inflated Poisson
          if( b_i(i)==0 ){
            //LogProb2_i(i) = log( (1-R1_i(i)) + dpois(Type(0.0), R2_i(i), false)*R1_i(i) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
            LogProb2_i(i) = logspace_add( log(1-R1_i(i)), dpois(Type(0.0),R2_i(i)*exp(SigmaM(e_i(i),0)*delta_i(i)-pow(SigmaM(e_i(i),0),2)/2),true)+log(R1_i(i)) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
          }else{
            LogProb2_i(i) = dpois(b_i(i), R2_i(i)*exp(SigmaM(e_i(i),0)*delta_i(i)-pow(SigmaM(e_i(i),0),2)/2), true) + log(R1_i(i)); // Pr[X=x] = Pois(X=x)*phi
          }
        }
        LogProb1_i(i) = 0;
      }
    }
  }
  REPORT( diag_iz );

  // Joint likelihood
  jnll_comp(10) = -1 * (LogProb1_i * (Type(1.0)-PredTF_i)).sum();
  jnll_comp(11) = -1 * (LogProb2_i * (Type(1.0)-PredTF_i)).sum();
  jnll = jnll_comp.sum();
  Type pred_jnll = -1 * ( LogProb1_i*PredTF_i + LogProb2_i*PredTF_i ).sum();
  REPORT( pred_jnll );
  REPORT( tmp_calc1 );
  REPORT( tmp_calc2 );

  ////////////////////////
  // Calculate outputs
  ////////////////////////

  // Number of output-years
  int n_y = t_yz.col(0).size();

  // Predictive distribution -- ObsModel_ez(e,0)==4 isn't implemented (it had a bug previously)
  Type a_average = a_i.sum()/a_i.size();
  array<Type> P1_xcy(n_x, n_c, n_y);
  array<Type> R1_xcy(n_x, n_c, n_y);
  array<Type> P2_xcy(n_x, n_c, n_y);
  array<Type> R2_xcy(n_x, n_c, n_y);
  array<Type> D_xcy(n_x, n_c, n_y);
  for(int c=0; c<n_c; c++){
  for(int y=0; y<n_y; y++){
  for(int x=0; x<n_x; x++){
    // Calculate linear predictors
    P1_xcy(x,c,y) = Omega1_sc(x,c) + eta1_x(x);
    P2_xcy(x,c,y) =  Omega2_sc(x,c) + eta2_x(x);
    for( int z=0; z<t_yz.row(0).size(); z++ ){
      if( t_yz(y,z)>=0 & t_yz(y,z)<n_t ){    // isNA doesn't seem to work for IMATRIX type
        P1_xcy(x,c,y) += beta1_ct(c,t_yz(y,z)) + Epsilon1_sct(x,c,t_yz(y,z))*exp(log_sigmaratio1_z(z)) + eta1_xct(x,c,t_yz(y,z));
        P2_xcy(x,c,y) += beta2_ct(c,t_yz(y,z)) + Epsilon2_sct(x,c,t_yz(y,z))*exp(log_sigmaratio2_z(z)) + eta2_xct(x,c,t_yz(y,z));
      }
    }
    // Calculate predictors in link-space
    if( ObsModel_ez(c,1)==0 ){
      R1_xcy(x,c,y) = invlogit( P1_xcy(x,c,y) );
      R2_xcy(x,c,y) = exp( P2_xcy(x,c,y) );
      D_xcy(x,c,y) = R1_xcy(x,c,y) * R2_xcy(x,c,y);
    }
    if( ObsModel_ez(c,1)==1 ){
      R1_xcy(x,c,y) = Type(1.0) - exp( -exp(P1_xcy(x,c,y)) );
      R2_xcy(x,c,y) = exp(P1_xcy(x,c,y)) / R1_xcy(x,c,y) * exp( P2_xcy(x,c,y) );
      D_xcy(x,c,y) = exp( P1_xcy(x,c,y) ) * exp( P2_xcy(x,c,y) );        // Use this line to prevent numerical over/underflow
    }
    if( ObsModel_ez(c,1)==2 ){
      R1_xcy(x,c,y) = exp( P1_xcy(x,c,y) );
      R2_xcy(x,c,y) = exp( P2_xcy(x,c,y) );
      D_xcy(x,c,y) = R1_xcy(x,c,y) * R2_xcy(x,c,y);
    }
  }}}

  // Calculate indices
  array<Type> Index_xcyl(n_x, n_c, n_y, n_l);
  array<Type> Index_cyl(n_c, n_y, n_l);
  array<Type> ln_Index_cyl(n_c, n_y, n_l);
  Index_cyl.setZero();
  for(int c=0; c<n_c; c++){
  for(int y=0; y<n_y; y++){
  for(int l=0; l<n_l; l++){
    for(int x=0; x<n_x; x++){
      Index_xcyl(x,c,y,l) = D_xcy(x,c,y) * a_xl(x,l) / 1000;  // Convert from kg to metric tonnes
      Index_cyl(c,y,l) += Index_xcyl(x,c,y,l);
    }
  }}}
  ln_Index_cyl = log( Index_cyl );

  // Calculate other derived summaries
  // Each is the weighted-average X_xl over polygons (x) with weights equal to abundance in each polygon and time (where abundance is from the first index)
  array<Type> mean_Z_cym(n_c, n_y, n_m);
  if( Options(2)==1 ){
    mean_Z_cym.setZero();
    int report_summary_TF = false;
    for(int c=0; c<n_c; c++){
    for(int y=0; y<n_y; y++){
    for(int m=0; m<n_m; m++){
      for(int x=0; x<n_x; x++){
        if( Z_xm(x,m)!=0 ){
          mean_Z_cym(c,y,m) += Z_xm(x,m) * Index_xcyl(x,c,y,0)/Index_cyl(c,y,0);
          report_summary_TF = true;
        }
      }
    }}}
    if( report_summary_TF==true ){
      REPORT( mean_Z_cym );
      ADREPORT( mean_Z_cym );
    }
  }

  // Calculate average density, weighted.mean( x=Abundance/Area, w=Abundance )
  // Doesn't require Z_xm, because it only depends upon Index_tl
  if( Options(4)==1 ){
    array<Type> mean_D_cyl(n_c, n_y, n_l);
    array<Type> log_mean_D_cyl(n_c, n_y, n_l);
    mean_D_cyl.setZero();
    for(int c=0; c<n_c; c++){
    for(int y=0; y<n_y; y++){
    for(int l=0; l<n_l; l++){
      for(int x=0; x<n_x; x++){
        mean_D_cyl(c,y,l) += D_xcy(x,c,y) * Index_xcyl(x,c,y,l)/Index_cyl(c,y,l);
      }
    }}}
    REPORT( mean_D_cyl );
    ADREPORT( mean_D_cyl );
    log_mean_D_cyl = log( mean_D_cyl );
    ADREPORT( log_mean_D_cyl );

    // Calculate effective area = Index / average density
    array<Type> effective_area_cyl(n_c, n_y, n_l);
    array<Type> log_effective_area_cyl(n_c, n_y, n_l);
    effective_area_cyl = Index_cyl / (mean_D_cyl/1000);  // Correct for different units of Index and density
    log_effective_area_cyl = log( effective_area_cyl );
    REPORT( effective_area_cyl );
    ADREPORT( effective_area_cyl );
    ADREPORT( log_effective_area_cyl );
  }

  // Reporting and standard-errors for covariance and correlation matrices
  if( Options(5)==1 ){
    if( FieldConfig(0)>0 ){
      matrix<Type> L1_omega_cf = loadings_matrix( L_omega1_z, n_c, FieldConfig(0) );
      matrix<Type> lowercov_uppercor_omega1 = L1_omega_cf * L1_omega_cf.transpose();
      lowercov_uppercor_omega1 = convert_upper_cov_to_cor( lowercov_uppercor_omega1 );
      REPORT( lowercov_uppercor_omega1 );
      ADREPORT( lowercov_uppercor_omega1 );
    }
    if( FieldConfig(1)>0 ){
      matrix<Type> L1_epsilon_cf = loadings_matrix( L_epsilon1_z, n_c, FieldConfig(1) );
      matrix<Type> lowercov_uppercor_epsilon1 = L1_epsilon_cf * L1_epsilon_cf.transpose();
      lowercov_uppercor_epsilon1 = convert_upper_cov_to_cor( lowercov_uppercor_epsilon1 );
      REPORT( lowercov_uppercor_epsilon1 );
      ADREPORT( lowercov_uppercor_epsilon1 );
    }
    if( FieldConfig(2)>0 ){
      matrix<Type> L2_omega_cf = loadings_matrix( L_omega2_z, n_c, FieldConfig(2) );
      matrix<Type> lowercov_uppercor_omega2 = L2_omega_cf * L2_omega_cf.transpose();
      lowercov_uppercor_omega2 = convert_upper_cov_to_cor( lowercov_uppercor_omega2 );
      REPORT( lowercov_uppercor_omega2 );
      ADREPORT( lowercov_uppercor_omega2 );
    }
    if( FieldConfig(3)>0 ){
      matrix<Type> L2_epsilon_cf = loadings_matrix( L_epsilon2_z, n_c, FieldConfig(3) );
      matrix<Type> lowercov_uppercor_epsilon2 = L2_epsilon_cf * L2_epsilon_cf.transpose();
      lowercov_uppercor_epsilon2 = convert_upper_cov_to_cor( lowercov_uppercor_epsilon2 );
      REPORT( lowercov_uppercor_epsilon2 );
      ADREPORT( lowercov_uppercor_epsilon2 );
    }
  }

  // Synchrony
  if( Options(6)==1 ){
    int n_z = yearbounds_zz.col(0).size();
    // Density ("D") or area-expanded total biomass ("B") for each category (use B when summing across sites)
    matrix<Type> D_xy( n_x, n_y );
    matrix<Type> B_cy( n_c, n_y );
    vector<Type> B_y( n_y );
    D_xy.setZero();
    B_cy.setZero();
    B_y.setZero();
    // Sample variance in category-specific density ("D") and biomass ("B")
    array<Type> varD_xcz( n_x, n_c, n_z );
    array<Type> varD_xz( n_x, n_z );
    array<Type> varB_cz( n_c, n_z );
    vector<Type> varB_z( n_z );
    vector<Type> varB_xbar_z( n_z );
    vector<Type> varB_cbar_z( n_z );
    vector<Type> ln_varB_z( n_z );
    vector<Type> ln_varB_xbar_z( n_z );
    vector<Type> ln_varB_cbar_z( n_z );
    array<Type> maxsdD_xz( n_x, n_z );
    array<Type> maxsdB_cz( n_c, n_z );
    vector<Type> maxsdB_z( n_z );
    varD_xcz.setZero();
    varD_xz.setZero();
    varB_cz.setZero();
    varB_z.setZero();
    varB_xbar_z.setZero();
    varB_cbar_z.setZero();
    maxsdD_xz.setZero();
    maxsdB_cz.setZero();
    maxsdB_z.setZero();
    // Proportion of total biomass ("P") for each location or each category
    matrix<Type> propB_xz( n_x, n_z );
    matrix<Type> propB_cz( n_c, n_z );
    propB_xz.setZero();
    propB_cz.setZero();
    // Synchrony indices
    matrix<Type> phi_xz( n_x, n_z );
    matrix<Type> phi_cz( n_c, n_z );
    vector<Type> phi_xbar_z( n_z );
    vector<Type> phi_cbar_z( n_z );
    vector<Type> phi_z( n_z );
    phi_xbar_z.setZero();
    phi_cbar_z.setZero();
    phi_z.setZero();
    // Calculate total biomass for different categories
    for( int y=0; y<n_y; y++ ){
      for( int c=0; c<n_c; c++ ){
        for( int x=0; x<n_x; x++ ){
          D_xy(x,y) += D_xcy(x,c,y);
          B_cy(c,y) += D_xcy(x,c,y) * a_xl(x,0);
          B_y(y) += D_xcy(x,c,y) * a_xl(x,0);
        }
      }
    }
    // Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
    Type temp_mean;
    for( int z=0; z<n_z; z++ ){
      for( int x=0; x<n_x; x++ ){
        // Variance for biomass in each category, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
        for( int c=0; c<n_c; c++ ){
          temp_mean = 0;
          for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) temp_mean += D_xcy(x,c,y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
          for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ){
            varD_xcz(x,c,z) += pow(D_xcy(x,c,y)-temp_mean,2) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0));
          }
        }
        // Variance for combined biomass across categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
        temp_mean = 0;
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) temp_mean += D_xy(x,y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
          varD_xz(x,z) += pow(D_xy(x,y)-temp_mean,2) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0));
        }
      }
      for( int c=0; c<n_c; c++ ){
        // Variance for combined biomass across sites, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
        temp_mean = 0;
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) temp_mean += B_cy(c,y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
          varB_cz(c,z) += pow(B_cy(c,y)-temp_mean,2) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0));
        }
      }
      // Variance for combined biomass across sites and categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
      temp_mean = 0;
      for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) temp_mean += B_y(y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
      for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
        varB_z(z) += pow(B_y(y)-temp_mean,2) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0));
      }
      // Proportion in each site
      for( int x=0; x<n_x; x++ ){
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
          propB_xz(x,z) += a_xl(x,0) * D_xy(x,y) / B_y(y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
        }
      }
      // Proportion in each category
      for( int c=0; c<n_c; c++ ){
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
          propB_cz(c,z) += B_cy(c,y) / B_y(y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
        }
      }
      // Species-buffering index (calculate in Density so that areas with zero area are OK)
      for( int x=0; x<n_x; x++ ){
        for( int c=0; c<n_c; c++ ){
          maxsdD_xz(x,z) += pow(varD_xcz(x,c,z), 0.5);
        }
        phi_xz(x,z) = varD_xz(x,z) / pow( maxsdD_xz(x,z), 2);
        varB_xbar_z(z) += pow(a_xl(x,0),2) * varD_xz(x,z) * propB_xz(x,z);
        phi_xbar_z(z) += phi_xz(x,z) * propB_xz(x,z);
      }
      // Spatial-buffering index
      for( int c=0; c<n_c; c++ ){
        for( int x=0; x<n_x; x++ ){
          maxsdB_cz(c,z) += a_xl(x,0) * pow(varD_xcz(x,c,z), 0.5);
        }
        phi_cz(c,z) = varB_cz(c,z) / pow( maxsdB_cz(c,z), 2);
        varB_cbar_z(z) += varB_cz(c,z) * propB_cz(c,z);
        phi_cbar_z(z) += phi_cz(c,z) * propB_cz(c,z);
      }
      // Spatial and species-buffering index
      for( int c=0; c<n_c; c++ ){
        for( int x=0; x<n_x; x++ ){
          maxsdB_z(z) += a_xl(x,0) * pow(varD_xcz(x,c,z), 0.5);
        }
      }
      phi_z(z) = varB_z(z) / pow( maxsdB_z(z), 2);
    }
    ln_varB_xbar_z = log( varB_xbar_z );
    ln_varB_cbar_z = log( varB_cbar_z );
    ln_varB_z = log( varB_z );
    REPORT( B_y );
    REPORT( D_xy );
    REPORT( B_cy );
    REPORT( phi_xz );
    REPORT( phi_xbar_z );
    REPORT( phi_cz );
    REPORT( phi_cbar_z );
    REPORT( phi_z );
    REPORT( propB_xz );
    REPORT( propB_cz );
    REPORT( varD_xcz );
    REPORT( varD_xz );
    REPORT( varB_cz );
    REPORT( varB_z );
    REPORT( varB_xbar_z );
    REPORT( varB_cbar_z );
    REPORT( maxsdB_z );
    REPORT( maxsdD_xz );
    REPORT( maxsdB_cz );
    ADREPORT( varB_xbar_z );
    ADREPORT( varB_cbar_z );
    ADREPORT( varB_z );
    ADREPORT( B_y );
    ADREPORT( ln_varB_xbar_z );
    ADREPORT( ln_varB_cbar_z );
    ADREPORT( ln_varB_z );
    ADREPORT( phi_xbar_z );
    ADREPORT( phi_cbar_z );
    ADREPORT( phi_z );
  }

  // Calculate coherence and variance and covariance matrices
  // psi: "Coherence" = degree to which covariance is explained by one or many factors (1/n_c, 1), 1=Single factor; 1/n_c=even factors
  if( Options(7)==1 ){
    // Eigendecomposition see: https://github.com/kaskr/adcomp/issues/144#issuecomment-228426834
    using namespace Eigen;
    // Spatio-temporal covariance (summation only works when ObsModel_ez[c,1]==1)
    //Matrix<Type,Dynamic,Dynamic> CovHat( n_c, n_c );
    matrix<Type> CovHat( n_c, n_c );
    CovHat.setIdentity();
    CovHat *= pow(0.0001, 2);
    if( FieldConfig(1)>0 ) CovHat += loadings_matrix(L_epsilon1_z, n_c, FieldConfig(1)) * loadings_matrix(L_epsilon1_z, n_c, FieldConfig(1)).transpose();
    if( FieldConfig(3)>0 ) CovHat += loadings_matrix(L_epsilon2_z, n_c, FieldConfig(3)) * loadings_matrix(L_epsilon2_z, n_c, FieldConfig(3)).transpose();
    // Coherence ranges from 0 (all factors are equal) to 1 (first factor explains all variance)
    SelfAdjointEigenSolver<Matrix<Type,Dynamic,Dynamic> > es(CovHat);
    vector<Type> eigenvalues_c = es.eigenvalues();       // Ranked from lowest to highest for some reason
    Type psi = 0;
    for(int c=0; c<n_c; c++) psi += eigenvalues_c(n_c-c-1) * (n_c - c);
    psi = 2 * ((psi / eigenvalues_c.sum() / n_c) - 0.5);
    // Total variance
    vector<Type> diag_CovHat( n_c );
    vector<Type> log_diag_CovHat( n_c );
    for(int c=0; c<n_c; c++) diag_CovHat(c) = CovHat(c,c);
    Type totalvar_CovHat = diag_CovHat.sum();
    Type log_totalvar_CovHat = log(totalvar_CovHat);
    log_diag_CovHat = log(diag_CovHat);
    // Reporting
    REPORT( CovHat );
    REPORT( psi );
    REPORT( eigenvalues_c );
    ADREPORT( psi );
    ADREPORT( diag_CovHat );
    ADREPORT( totalvar_CovHat );
    ADREPORT( log_totalvar_CovHat );
    ADREPORT( log_diag_CovHat );
    ADREPORT( eigenvalues_c );
  }

  // Calculate proportion of biomass for each category
  if( Options(8)==1 ){
    array<Type> PropIndex_cyl(n_c, n_y, n_l);
    array<Type> ln_PropIndex_cyl(n_c, n_y, n_l);
    Type sumtemp;
    for(int y=0; y<n_y; y++){
    for(int l=0; l<n_l; l++){                 // .col(0).cols()
      sumtemp = 0;
      for(int c=0; c<n_c; c++){
        sumtemp += Index_cyl(c,y,l);
      }
      for(int c=0; c<n_c; c++){
        PropIndex_cyl(c,y,l) = Index_cyl(c,y,l) / sumtemp;
      }
    }}
    ln_PropIndex_cyl = log( PropIndex_cyl );
    REPORT( PropIndex_cyl );
    REPORT( ln_PropIndex_cyl );
    ADREPORT( PropIndex_cyl );
    ADREPORT( ln_PropIndex_cyl );
  }

  // Diagnostic output
  REPORT( Q1 );
  REPORT( Q2 );
  REPORT( P1_iz );
  REPORT( P2_iz );
  REPORT( R1_i );
  REPORT( R2_i );
  REPORT( P1_xcy );
  REPORT( P2_xcy );
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
  REPORT( Index_cyl );
  REPORT( D_xcy );
  REPORT( R1_xcy );
  REPORT( R2_xcy );
  REPORT( Index_xcyl );
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
  ADREPORT( Index_cyl );
  ADREPORT( ln_Index_cyl );
  ADREPORT( SigmaM );

  // Additional miscellaneous outputs
  if( Options(0)==1 ){
    ADREPORT( Index_xcyl );
  }
  if( Options(1)==1 ){
    ADREPORT( log(Index_xcyl) );
  }

  return jnll;
  
}
