#include <TMB.hpp>
#include <Eigen/Eigenvalues>

// Function to import R list for user-defined Options_vec and Options, packaged as list Options_list in TmbData
template<class Type>
struct options_list {
  vector<int> Options_vec;
  vector<int> Options;
  matrix<int> yearbounds_zz;
  matrix<int> Expansion_cz;
  options_list(SEXP x){ // Constructor
    Options_vec = asVector<int>(getListElement(x,"Options_vec"));
    Options = asVector<int>(getListElement(x,"Options"));
    yearbounds_zz = asMatrix<int>(getListElement(x,"yearbounds_zz"));
    Expansion_cz = asMatrix<int>(getListElement(x,"Expansion_cz"));
  }
};

// Needed for returning SparseMatrix
template<class Type>
Eigen::SparseMatrix<Type> Q_network( Type log_theta, int n_s, vector<int> parent_s, vector<int> child_s, vector<Type> dist_s ){
  Eigen::SparseMatrix<Type> Q( n_s, n_s );
  Type theta = exp( log_theta );
  for(int s=0; s<n_s; s++){
    Q.coeffRef( s, s ) = Type(1.0);
  }
  for(int s=1; s<parent_s.size(); s++){
    if( exp(-dist_s(s))!=0 ){
      Q.coeffRef( parent_s(s), child_s(s) ) = -exp(-theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
      Q.coeffRef( child_s(s), parent_s(s) ) = Q.coeffRef( parent_s(s), child_s(s) );
      Q.coeffRef( parent_s(s), parent_s(s) ) += exp(-2*theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
      Q.coeffRef( child_s(s), child_s(s) ) += exp(-2*theta*dist_s(s)) / (1-exp(-2*theta*dist_s(s)));
    }
  }
  return Q;
}

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
// eta_jf could be either eta_vf (for overdispersion) or eta_tf (for year effects)
template<class Type>
matrix<Type> covariation_by_category_nll( int n_f, int n_j, int n_c, matrix<Type> eta_jf, matrix<Type> eta_mean_jf, vector<Type> L_z, Type &jnll_pointer, objective_function<Type>* of){
  using namespace density;
  matrix<Type> eta_jc(n_j, n_c);
  vector<Type> Tmp_c;
  // IID
  if( n_f == -2 ){
    for( int j=0; j<n_j; j++ ){
    for( int c=0; c<n_c; c++ ){
      int f = c;
      jnll_pointer -= dnorm( eta_jf(j,f), eta_mean_jf(j,f), Type(1.0), true );
      // Simulate new values when using obj.simulate()
      if(isDouble<Type>::value && of->do_simulate) {
        eta_jf(j,f) = rnorm( eta_mean_jf(j,f), Type(1.0) );
      }
      // Rescale
      eta_jc(j,c) = eta_jf(j,f) * L_z(f);
    }}
  }
  // Turn off
  if( n_f == -1 ){
    eta_jc.setZero();
  }
  // AR1 structure
  if( n_f==0 ){
    for( int j=0; j<n_j; j++ ){
      Tmp_c = eta_jf.row(j);
      jnll_pointer += SCALE( AR1(L_z(1)), exp(L_z(0)) )( Tmp_c );
      // Simulate new values when using obj.simulate()
      if(isDouble<Type>::value && of->do_simulate){
        SCALE( AR1(L_z(1)), exp(L_z(0)) ).simulate(Tmp_c);
        eta_jf.row(j) = Tmp_c;
      }
    }
    eta_jc = eta_jf;
  }
  // Factor analysis structure
  if( n_f>0 ){
    // Assemble the loadings matrix
    matrix<Type> L_cf = loadings_matrix( L_z, n_c, n_f );
    // Probability of overdispersion
    for( int j=0; j<n_j; j++ ){
    for( int f=0; f<n_f; f++ ){
      jnll_pointer -= dnorm( eta_jf(j,f), eta_mean_jf(j,f), Type(1.0), true );
      // Simulate new values when using obj.simulate()
      if(isDouble<Type>::value && of->do_simulate){
        eta_jf(j,f) = rnorm( eta_mean_jf(j,f), Type(1.0) );
      }
    }}
    // Multiply out overdispersion
    eta_jc = eta_jf * L_cf.transpose();
  }
  return eta_jc;
}

template<class Type>                                                                                        //
matrix<Type> convert_upper_cov_to_cor( matrix<Type> cov ){
  int nrow = cov.rows();
  for( int i=0; i<nrow; i++){
  for( int j=i+1; j<nrow; j++){
    cov(i,j) = cov(i,j) / pow(cov(i,i),0.5) / pow(cov(j,j),0.5);
  }}
  return cov;
}

// Input:  n_g, n_f, n_t, logical_flag, Ags_ij, Ags_x, Mat_sc
template<class Type>                                                                                        //
array<Type> project_knots( int n_g, int n_f, int n_t, int is_epsilon, array<Type> Mat_sft, matrix<int> A_ij, vector<Type> A_x ){
  array<Type> Mat_gf(n_g, n_f);
  array<Type> Mat_gft(n_g, n_f, n_t);
  if( is_epsilon!=1 ) Mat_gf.setZero();
  if( is_epsilon==1 ) Mat_gft.setZero();
  for( int t=0; t<n_t; t++ ){
  for( int Arow=0; Arow<A_ij.rows(); Arow++ ){
  for( int f=0; f<n_f; f++ ){
    int g = A_ij(Arow,0);
    int s = A_ij(Arow,1);
    if( is_epsilon!=1 ) Mat_gf(g,f) += A_x(Arow) * Mat_sft(s,f);
    if( is_epsilon==1 ) Mat_gft(g,f,t) += A_x(Arow) * Mat_sft(s,f,t);
  }}}
  if( is_epsilon!=1 ){return Mat_gf;}else{return Mat_gft;}
}


// Input: L_omega1_z, Q1, Omegainput1_sf, n_f, n_s, n_c, FieldConfig(0,0)
// Output: jnll_comp(0), Omega1_sc
template<class Type>                                                                                        //
matrix<Type> gmrf_by_category_nll( int n_f, int method, int timing, int n_s, int n_c, Type logkappa, array<Type> gmrf_input_sf, array<Type> gmrf_mean_sf, vector<Type> L_z, density::GMRF_t<Type> gmrf_Q, Type &jnll_pointer, objective_function<Type>* of){
  using namespace density;
  matrix<Type> gmrf_sc(n_s, n_c);
  vector<Type> gmrf_s(n_s);
  matrix<Type> Cov_cc(n_c,n_c);
  array<Type> diff_gmrf_sc(n_s, n_c); // Requires an array
  Type logtau;
  if(method==0) logtau = log( 1 / (exp(logkappa) * sqrt(4*M_PI)) );
  if(method==1) logtau = log( 1 / sqrt(1-exp(logkappa*2)) );
  if( (method!=0) & (method!=1) ) logtau = Type(0.0);
  // IID
  if(n_f == -2){
    for( int c=0; c<n_c; c++ ){
      int f = c;
      jnll_pointer += gmrf_Q(gmrf_input_sf.col(f) - gmrf_mean_sf.col(f));
      // Simulate new values when using obj.simulate()
      if(isDouble<Type>::value && of->do_simulate) {
        gmrf_Q.simulate(gmrf_s);
        gmrf_input_sf.col(f) = gmrf_s + gmrf_mean_sf.col(f);
      }
      // Rescale
      gmrf_sc.col(c) = gmrf_input_sf.col(f) / exp(logtau) * L_z(f);                                // Rescaling from comp_index_v1d.cpp
    }
  }
  // Turn off
  if(n_f == -1){
    gmrf_sc.setZero();
  }
  // AR1 structure
  if(n_f==0){
    jnll_pointer += SEPARABLE( AR1(L_z(1)), gmrf_Q )(gmrf_input_sf - gmrf_mean_sf);
    // Simulate new values when using obj.simulate()
    if(isDouble<Type>::value && of->do_simulate) {
      SEPARABLE( AR1(L_z(1)), gmrf_Q ).simulate(gmrf_input_sf);
      gmrf_input_sf += gmrf_input_sf;
    }
    // Rescale
    logtau = L_z(0) - logkappa;  //
    gmrf_sc = gmrf_input_sf / exp(logtau);                                // Rescaling from comp_index_v1d.cpp
  }
  // Factor analysis structure
  if(n_f>0){
    // PDF if density-dependence/interactions occurs prior to correlated dynamics
    if( timing==0 ){
      for( int f=0; f<n_f; f++ ){
        // Calculate likelihood
        jnll_pointer += gmrf_Q(gmrf_input_sf.col(f) - gmrf_mean_sf.col(f));  // Rescaling from spatial_vam_v13.cpp
        // Simulate new values when using obj.simulate()
        if(isDouble<Type>::value && of->do_simulate) {
          gmrf_Q.simulate(gmrf_s);
          gmrf_input_sf.col(f) = gmrf_s + gmrf_mean_sf.col(f);
        }
      }
      // Rescale
      matrix<Type> L_cf = loadings_matrix( L_z, n_c, n_f );
      gmrf_sc = (gmrf_input_sf.matrix() * L_cf.transpose()) / exp(logtau);
    }
    // PDF if density-dependence/interactions occurs after correlated dynamics (Only makes sense if n_f == n_c)
    if( timing==1 ){
      // Calculate difference without rescaling
      gmrf_sc = gmrf_input_sf.matrix();
      for( int s=0; s<n_s; s++){
      for( int c=0; c<n_c; c++){
        diff_gmrf_sc(s,c) = gmrf_sc(s,c) - gmrf_mean_sf(s,c);
      }}
      // Calculate likelihood
      matrix<Type> L_cf = loadings_matrix( L_z, n_c, n_f );
      Cov_cc = L_cf * L_cf.transpose();
      jnll_pointer += SCALE(SEPARABLE(MVNORM(Cov_cc), gmrf_Q), exp(-logtau))( diff_gmrf_sc );
      //gmrf_sc = gmrf_sc / exp(logtau);
      // Simulate new values when using obj.simulate()
      if(isDouble<Type>::value && of->do_simulate) {
        SEPARABLE(MVNORM(Cov_cc), gmrf_Q).simulate( diff_gmrf_sc );
        gmrf_sc = gmrf_mean_sf + diff_gmrf_sc/exp(logtau);
      }
    }
  }
  return gmrf_sc;
}

// Used to calculate GMRF PDF for initial condition given covariance Cov_cc
// Only makes sense given:
// 1. full-rank factor model
// 2. Spatial Gompertz model conditions
// 3. Timing = 1
template<class Type>
matrix<Type> gmrf_stationary_nll( int method, int n_s, int n_c, Type logkappa, array<Type> gmrf_input_sc, matrix<Type> Cov_cc, density::GMRF_t<Type> gmrf_Q, Type &jnll_pointer, objective_function<Type>* of){
  using namespace density;
  array<Type> gmrf_sc(n_s, n_c);
  Type logtau;
  if(method==0) logtau = log( 1 / (exp(logkappa) * sqrt(4*M_PI)) );
  if(method==1) logtau = log( 1 / sqrt(1-exp(logkappa*2)) );
  if( (method!=0) & (method!=1) ) logtau = Type(0.0);
  // PDF if density-dependence/interactions occurs after correlated dynamics (Only makes sense if n_f == n_c)
  gmrf_sc = gmrf_input_sc.matrix();
  // Calculate likelihood
  jnll_pointer += SCALE(SEPARABLE(MVNORM(Cov_cc), gmrf_Q), exp(-logtau))( gmrf_sc );
  // Simulate new values when using obj.simulate()
  if(isDouble<Type>::value && of->do_simulate) {
    SEPARABLE(MVNORM(Cov_cc), gmrf_Q).simulate( gmrf_sc );
    gmrf_sc = gmrf_sc / exp(logtau);
  }
  return gmrf_sc.matrix();
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
  Type log_w_j;
  //Type pos_penalty;
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

// Calculate B_cc
template<class Type>
matrix<Type> calculate_B( int method, int n_f, int n_r, matrix<Type> Chi_fr, matrix<Type> Psi_fr, Type &jnll_pointer ){
  matrix<Type> B_ff( n_f, n_f );
  matrix<Type> BplusI_ff( n_f, n_f );
  matrix<Type> Chi_rf = Chi_fr.transpose();
  matrix<Type> Psi_rf = Psi_fr.transpose();
  matrix<Type> Identity_ff( n_f, n_f );
  Identity_ff.setIdentity();

  // No interactions (default)
  if( method==0 ){
    B_ff.setZero();
  }
  // Simple co-integration -- complex unbounded eigenvalues
  if( method==1 ){
    B_ff = Chi_fr * Psi_rf;
  }
  // Real eigenvalues
  if( method==2 ){
    matrix<Type> Chi_ff( n_f, n_f );
    Chi_ff = Identity_ff;
    // Make Chi_ff
    vector<Type> colnorm_r( n_r );
    colnorm_r.setZero();
    for(int f=0; f<n_f; f++){
    for(int r=0; r<n_r; r++){
      Chi_ff(f,r) = Chi_fr(f,r);
      colnorm_r(r) += pow( Chi_ff(f,r), 2 );
    }}
    for(int f=0; f<n_f; f++){
    for(int r=0; r<n_r; r++){
      Chi_ff(f,r) /= pow( colnorm_r(r), 0.5 );
    }}
    // Make Psi_ff
    matrix<Type> Psi_ff( n_f, n_f );
    Psi_ff = Identity_ff;
    for(int f=n_r; f<n_f; f++){
    for(int r=0; r<n_r; r++){
      Psi_ff(f,r) = Psi_fr(f,r);
    }}
    // Make L_ff
    matrix<Type> L_ff(n_f, n_f);
    L_ff.setZero();
    for(int r=0; r<n_r; r++){
      L_ff(r,r) = Psi_fr(r,r);
    }
    // Build B_ff
    matrix<Type> invChi_ff = atomic::matinv( Chi_ff );
    matrix<Type> trans_Psi_ff = Psi_ff.transpose();
    matrix<Type> trans_invPsi_ff = atomic::matinv( Psi_ff ).transpose();
    B_ff = Chi_ff * trans_Psi_ff;
    B_ff = B_ff * L_ff;
    B_ff = B_ff * trans_invPsi_ff;
    B_ff = B_ff * invChi_ff;
    // Penalize colnorm_r
    jnll_pointer += ( log(colnorm_r)*log(colnorm_r) ).sum();
  }
  // Complex bounded eigenvalues
  if( method==3 ){
    BplusI_ff = Chi_fr * Psi_rf + Identity_ff;
    // Extract eigenvalues
    vector< std::complex<Type> > eigenvalues_B_ff = B_ff.eigenvalues();
    vector<Type> real_eigenvalues_B_ff = eigenvalues_B_ff.real();
    vector<Type> imag_eigenvalues_B_ff = eigenvalues_B_ff.imag();
    vector<Type> mod_eigenvalues_B_ff( n_f );
    // Calculate maximum eigenvalues
    Type MaxEigen = 1;
    for(int f=0; f<n_f; f++){
      mod_eigenvalues_B_ff(f) = pow( pow(real_eigenvalues_B_ff(f),2) + pow(imag_eigenvalues_B_ff(f),2), 0.5 );
      MaxEigen = CppAD::CondExpGt(mod_eigenvalues_B_ff(f), MaxEigen, mod_eigenvalues_B_ff(f), MaxEigen);
    }
    // Rescale interaction matrix
    BplusI_ff = BplusI_ff / MaxEigen;
    B_ff = BplusI_ff - Identity_ff;
    jnll_pointer += CppAD::CondExpGe( MaxEigen, Type(1.0), pow(MaxEigen-Type(1.0),2), Type(0.0) );
  }
  return B_ff;
}

// Calculate variance of stationary distributino
template<class Type>
matrix<Type> stationary_variance( int n_c, matrix<Type> B_cc, matrix<Type> Cov_cc ){
  int n2_c = n_c*n_c;
  matrix<Type> Kronecker_c2c2(n2_c,n2_c);
  matrix<Type> InvDiff_c2c2(n2_c, n2_c);
  matrix<Type> Vinf_cc(n_c, n_c);
  Kronecker_c2c2 = kronecker( B_cc, B_cc );
  InvDiff_c2c2.setIdentity();
  InvDiff_c2c2 = InvDiff_c2c2 - Kronecker_c2c2;
  InvDiff_c2c2 = atomic::matinv( InvDiff_c2c2 );
  Vinf_cc.setZero();
  for(int i=0; i<n_c; i++){
  for(int j=0; j<n_c; j++){
    int k = i + j*n_c;
    for(int i2=0; i2<n_c; i2++){
    for(int j2=0; j2<n_c; j2++){
      int k2 = i2 + j2*n_c;
      Vinf_cc(i,j) += InvDiff_c2c2(k,k2) * Cov_cc(i2,j2);
    }}
  }}
  return Vinf_cc;
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
  DATA_INTEGER(n_g);         // Number of extrapolation-grid cells
  DATA_INTEGER(n_t);         // Number of time-indices
  DATA_INTEGER(n_c);         // Number of categories (e.g., length bins)
  DATA_INTEGER(n_e);         // Number of error distributions
  DATA_INTEGER(n_p);         // Number of dynamic covariates
  DATA_INTEGER(n_v);          // Number of tows/vessels (i.e., levels for the factor explaining overdispersion)
  DATA_INTEGER(n_l);         // Number of indices to post-process
  DATA_INTEGER(n_m);         // Number of range metrics to use (probably 2 for Eastings-Northings)

  // Config
  DATA_STRUCT( Options_list, options_list );
  // Options_list.Options_vec
    // Slot 0 -- Aniso: 0=No, 1=Yes
    // Slot 1 -- DEPRECATED
    // Slot 2 -- AR1 on beta1 (year intercepts for 1st linear predictor) to deal with missing years: 0=No, 1=Yes
    // Slot 3 -- AR1 on beta2 (year intercepts for 2nd linear predictor) to deal with missing years: 0=No, 1=Yes
    // Slot 4 -- DEPRECATED
    // Slot 5 -- Upper limit constant of integration calculation for infinite-series density functions (Conway-Maxwell-Poisson and Tweedie)
    // Slot 6 -- Breakpoint in CMP density function
    // Slot 7 -- Whether to use SPDE or 2D-AR1 hyper-distribution for spatial process: 0=SPDE; 1=2D-AR1; 2=Stream-network
    // Slot 8 -- Whether to use F_ct or ignore it for speedup
  // Options_list.Options
    // Slot 0: Calculate SE for Index_xctl
    // Slot 1: Calculate SE for log(Index_xctl)
    // Slot 2: Calculate mean_Z_ctm (i.e., center-of-gravity)
    // Slot 3: Calculate SE for D_i (expected density for every observation)
    // Slot 4: Calculate mean_D_tl and effective_area_tl
    // Slot 5: Calculate standard errors for Covariance and Correlation among categories using factor-analysis parameterization
    // Slot 6: Calculate synchrony for different periods specified via yearbounds_zz
    // Slot 7: Calculate coherence and variance for Epsilon1_sct and Epsilon2_sct
    // Slot 8: Calculate proportions and SE
    // Slot 9: Include normalization in GMRF PDF
    // Slot 10: Calculate Fratio as F_ct divided by F achieving 40% of B0
    // Slot 11: Calculate B0 and Bratio
    // Slot 12: Calculate Omegainput1_gf, Omegainput2_gf, Epsiloninput1_gft, Epsiloninput1_gft
  // Options_list.yearbounds_zz
    // Two columns, and 1+ rows, specifying first and last t for each period used in calculating synchrony
  // Options_list.Expansion_cz
    // Two columns and n_c rows.  1st column:  Type of expansion (0=area-expansion; 1=biomass-expansion);  2nd column:  Category used for biomass-expansion
  DATA_IMATRIX(FieldConfig);  // Input settings (vector, length 4)
  DATA_IVECTOR(RhoConfig);
  DATA_IVECTOR(OverdispersionConfig);          // Input settings (vector, length 2)
  DATA_IMATRIX(ObsModel_ez);    // Observation model
  // Column 0: Probability distribution for data for each level of e_i
  // Column 1: Link function for linear predictors for each level of c_i
  // NOTE:  nlevels(c_i) must be <= nlevels(e_i)
  DATA_IVECTOR(VamConfig);
  // Slot 0 -- method for calculating n_c-by-n_c interaction matrix, B_ff
  // Slot 1 -- rank of interaction matrix B_ff
  // Slot 2 -- Timing of interactions;  0=Before correlated dynamics;  1=After correlated dynamics
  // Current implementation only makes sense when (1) intercepts are constant among years; (2) using a Poisson-link delta model; (3) n_f=n_c for spatio-temporal variation; (4) starts near equilibrium manifold
  DATA_IARRAY(Xconfig_zcp);
  // Row 0 -- Methods for 1st component for each covariate in X_xtp (0=Off;  1=Estimate;  2=Estimate with spatially varying coefficient)
  // Row 1 -- Methods for 2nd component for each covariate in X_xtp (0=Off;  1=Estimate;  2=Estimate with spatially varying coefficient)
  DATA_INTEGER(include_data);   // Always use TRUE except for internal usage to extract GRMF normalization when turn off GMRF normalization in CPP

  // Data vectors
  DATA_VECTOR(b_i);       	// Response (biomass) for each observation
  DATA_VECTOR(a_i);       	// Area swept for each observation (km^2)
  DATA_IMATRIX(c_iz);         // Category for each observation
  DATA_IVECTOR(e_i);         // Error distribution for each observation
  DATA_IMATRIX(t_iz);          // Time-indices (year, season, etc.) for each observation
  DATA_IVECTOR(v_i);          // tows/vessels for each observation (level of factor representing overdispersion)
  DATA_VECTOR(PredTF_i);          // vector indicating whether an observatino is predictive (1=used for model evaluation) or fitted (0=used for parameter estimation)
  DATA_MATRIX(a_gl);		     // Area for each "real" stratum(km^2) in each stratum
  DATA_ARRAY(X_itp);		    // Covariate design matrix (strata x covariate)
  DATA_ARRAY(X_gtp);		    // Covariate design matrix (strata x covariate)
  DATA_MATRIX(Q_ik);        // Catchability matrix (observations x variable)
  DATA_IMATRIX(t_yz);        // Matrix for time-indices of calculating outputs (abundance index and "derived-quantity")
  DATA_MATRIX(Z_gm);        // Derived quantity matrix
  DATA_MATRIX(F_ct);         // Matrix of annual fishing mortality for each category

  // Spatial network inputs
  DATA_IVECTOR(parent_s);  // Columns:  0=Parent index, 1=Child index, 2=Distance from parent to child
  DATA_IVECTOR(child_s);  // Columns:  0=Parent index, 1=Child index, 2=Distance from parent to child
  DATA_VECTOR(dist_s);  // Columns:  0=Parent index, 1=Child index, 2=Distance from parent to child

  // SPDE objects
  DATA_STRUCT(spde,spde_t);
  
  // Aniso objects
  DATA_STRUCT(spde_aniso,spde_aniso_t);

  // Sparse matrices for precision matrix of 2D AR1 process
  // Q = M0*(1+rho^2)^2 + M1*(1+rho^2)*(-rho) + M2*rho^2
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);

  // Projection matrices from knots s to data i or extrapolation-grid cells x
  //DATA_SPARSE_MATRIX(A_is);
  //DATA_SPARSE_MATRIX(A_gs);
  DATA_IMATRIX( Ais_ij );
  DATA_VECTOR( Ais_x );
  DATA_IMATRIX( Ags_ij );
  DATA_VECTOR( Ags_x );

  // Parameters 
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  PARAMETER_MATRIX(Chi_fr);   // error correction responses
  PARAMETER_MATRIX(Psi_fr);   // error correction loadings, B_ff = Chi_fr %*% t(Psi_fr)

  //  -- presence/absence fixed effects
  PARAMETER_MATRIX(beta1_ft);       // Year effect
  PARAMETER_ARRAY(gamma1_ctp);       // Dynamic covariate effect
  PARAMETER_VECTOR(lambda1_k);       // Catchability coefficients
  PARAMETER_VECTOR(L1_z);          // Overdispersion parameters
  PARAMETER_VECTOR(L_omega1_z);
  PARAMETER_VECTOR(L_epsilon1_z);
  PARAMETER_VECTOR(L_beta1_z);
  PARAMETER(logkappa1);
  PARAMETER_VECTOR(Beta_mean1_c);  // mean-reversion for beta1_ft
  PARAMETER_VECTOR(Beta_rho1_f);  // AR1 for presence/absence Beta component, Default=0
  PARAMETER_VECTOR(Epsilon_rho1_f);  // AR1 for presence/absence Epsilon component, Default=0
  PARAMETER_ARRAY(log_sigmaXi1_cp);  // log-SD of Xi1_scp
  PARAMETER_VECTOR(log_sigmaratio1_z);  // Ratio of variance for columns of t_iz

  // -- presence/absence random effects
  PARAMETER_MATRIX(eta1_vf);
  PARAMETER_ARRAY(Xiinput1_scp);         // spatially varying coefficient
  PARAMETER_ARRAY(Omegainput1_sf);      // Expectation
  PARAMETER_ARRAY(Epsiloninput1_sft);   // Annual variation

  //  -- positive catch rates fixed effects
  PARAMETER_MATRIX(beta2_ft);  // Year effect
  PARAMETER_ARRAY(gamma2_ctp);       // Dynamic covariate effect
  PARAMETER_VECTOR(lambda2_k);       // Catchability coefficients
  PARAMETER_VECTOR(L2_z);          // Overdispersion parameters
  PARAMETER_VECTOR(L_omega2_z);
  PARAMETER_VECTOR(L_epsilon2_z);
  PARAMETER_VECTOR(L_beta2_z);
  PARAMETER(logkappa2);
  PARAMETER_VECTOR(Beta_mean2_c);  // mean-reversion for beta2_t
  PARAMETER_VECTOR(Beta_rho2_f);  // AR1 for positive catch Beta component, Default=0
  PARAMETER_VECTOR(Epsilon_rho2_f);  // AR1 for positive catch Epsilon component, Default=0
  PARAMETER_ARRAY(log_sigmaXi2_cp);  // log-SD of Xi2_scp
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
  PARAMETER_ARRAY(Xiinput2_scp);         // spatially varying coefficient
  PARAMETER_ARRAY(Omegainput2_sf);      // Expectation
  PARAMETER_ARRAY(Epsiloninput2_sft);   // Annual variation

  ////////////////////////
  // Preparatory bookkeeping
  ////////////////////////

  // Indices -- i=Observation; t=Year; c=Category; p=Dynamic-covariate
  int i,t,c,p,s,g;
  
  // Objective function
  vector<Type> jnll_comp(16);
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
  // Slot 13 -- penalty on estimate_B structure
  // Slot 14 -- Spatially varying coefficient, encounter
  // Slot 15 -- Spatially varying coefficient, positive catch
  jnll_comp.setZero();
  Type jnll = 0;                

  // Unpack Options_list
  vector<int> Options_vec( Options_list.Options_vec.size() );
  Options_vec = Options_list.Options_vec;
  vector<int> Options( Options_list.Options.size() );
  Options = Options_list.Options;
  matrix<int> yearbounds_zz( Options_list.yearbounds_zz.rows(), 2 );
  yearbounds_zz = Options_list.yearbounds_zz;
  matrix<int> Expansion_cz( n_c, 2 );
  Expansion_cz = Options_list.Expansion_cz;

  // Derived parameters
  Type Range_raw1, Range_raw2;
  if( Options_vec(7)==0 ){
    Range_raw1 = sqrt(8) / exp( logkappa1 );   // Range = approx. distance @ 10% correlation
    Range_raw2 = sqrt(8) / exp( logkappa2 );     // Range = approx. distance @ 10% correlation
  }
  if( (Options_vec(7)==1) | (Options_vec(7)==2) ){
    Range_raw1 = log(0.1) / logkappa1;   // Range = approx. distance @ 10% correlation
    Range_raw2 = log(0.1) / logkappa2;     // Range = approx. distance @ 10% correlation
  }
  array<Type> SigmaM( n_e, 3 );
  array<Type> sigmaXi1_cp( n_c, n_p );
  array<Type> sigmaXi2_cp( n_c, n_p );
  SigmaM = exp( logSigmaM );
  sigmaXi1_cp = exp( log_sigmaXi1_cp );
  sigmaXi2_cp = exp( log_sigmaXi2_cp );

  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));

  // Overwrite parameters when mirroring them
  if( RhoConfig(1)==6 ){
    Beta_rho2_f = Beta_rho1_f;
  }
  if( RhoConfig(3)==6 ){
    Epsilon_rho2_f = Epsilon_rho1_f;
  }

  ////////////////////////
  // Interactions and fishing mortality
  ////////////////////////

  // Define interaction matrix for Epsilon1, and also the impact of F_ct on intercepts
  int n_f1;
  n_f1 = Epsiloninput1_sft.col(0).cols();
  int n_f2;
  n_f2 = Epsiloninput2_sft.col(0).cols();
  matrix<Type> B_ff( n_f1, n_f1 );          // Interactions among factors
  B_ff = calculate_B( VamConfig(0), n_f1, VamConfig(1), Chi_fr, Psi_fr, jnll_comp(13) );
  matrix<Type> iota_ct( n_c, n_t );       // Cumulative impact of fishing mortality F_ct in years <= current year t
  matrix<Type> B1_cc( n_c, n_c );        // Interactions among categories
  matrix<Type> covE1_cc( n_c, n_c );
  matrix<Type> B2_cc( n_c, n_c );        // Interactions among categories
  matrix<Type> covE2_cc( n_c, n_c );
  matrix<Type> I_cc( n_c, n_c );
  matrix<Type> IminusB_cc( n_c, n_c );
  I_cc.setIdentity();
  B1_cc.setZero();
  B2_cc.setZero();
  covE1_cc.setZero();
  covE2_cc.setZero();
  // Calculate interaction matrix B_cc for categories if feasible
  if( (n_c==n_f1) & (n_c==n_f2) & (FieldConfig(1,0)>0) & (FieldConfig(1,1)>0) ){
    matrix<Type> L_epsilon1_cf = loadings_matrix( L_epsilon1_z, n_c, n_f1 );
    matrix<Type> Cov_epsilon1_cc = L_epsilon1_cf * L_epsilon1_cf.transpose();
    matrix<Type> L_epsilon2_cf = loadings_matrix( L_epsilon2_z, n_c, n_f2 );
    matrix<Type> Cov_epsilon2_cc = L_epsilon2_cf * L_epsilon2_cf.transpose();
    matrix<Type> Btemp_cc( n_c, n_c );
    // Assemble interaction matrix
    B1_cc = B_ff;
    for( c=0; c<n_c; c++ ){
      B1_cc(c,c) += Epsilon_rho1_f(c);
    }
    // If Timing=0, transform from interaction among factors to interaction among categories
    if( VamConfig(2)==0 ){
      Btemp_cc = L_epsilon1_cf * B1_cc;
      B1_cc = Btemp_cc * L_epsilon1_cf.inverse();
    }
    // Assemble interaction matrix
    B2_cc = B_ff;
    for( c=0; c<n_c; c++ ){
      B2_cc(c,c) += Epsilon_rho2_f(c);
    }
    // If Timing=0, transform from interaction among factors to interaction among categories
    if( VamConfig(2)==0 ){
      Btemp_cc = L_epsilon2_cf * B2_cc;
      B2_cc = Btemp_cc * L_epsilon2_cf.inverse();
    }
    REPORT( B1_cc );
    REPORT( L_epsilon1_cf );
    REPORT( B2_cc );
    REPORT( L_epsilon2_cf );
    ADREPORT( B1_cc );
    // Calculate F resulting in 40% of B0 if requested (only makes sense when B1_cc = B2_cc or Epsilon2 is turned off)
    if( Options(10)==1 ){
      vector<Type> Btarg_c( n_c );
      vector<Type> Ftarg_c( n_c );
      matrix<Type> Fratio_ct( n_c, n_t );
      IminusB_cc = I_cc - B1_cc;
      Btarg_c = log( 0.4 );  // 40% target, transformed for log-link
      Ftarg_c = -1 * ( IminusB_cc * Btarg_c );
      for( t=0; t<n_t; t++ ){
      for( c=0; c<n_c; c++ ){
        Fratio_ct(c,t) = F_ct(c,t) / Ftarg_c(c);
      }}
      REPORT( Ftarg_c );
      REPORT( Fratio_ct );
      ADREPORT( Ftarg_c );
      ADREPORT( Fratio_ct );
    }
    // Calculate variance of stationary distribution only if necessary to calculate B0
    if( Options(11)==1 ){
      covE1_cc = stationary_variance( n_c, B1_cc, Cov_epsilon1_cc );
      REPORT( covE1_cc );
      covE2_cc = stationary_variance( n_c, B2_cc, Cov_epsilon2_cc );
      REPORT( covE2_cc );
    }
    // Define impact of F_ct on intercepts
    if( Options_vec(8)==0 ){
      iota_ct.setZero();
    }
    // Use F_ct in first year as initial condition...
    if( (Options_vec(8)==1) ){
      iota_ct.col(0) = -1 * F_ct.col(0);
    }
    // ... or use median of stationary distribution given F_ct in first year as initial condition
    if( Options_vec(8)==2 ){
      matrix<Type> sumB1_cc( n_c, n_c );
      IminusB_cc = I_cc - B1_cc;
      sumB1_cc = IminusB_cc.inverse();
      iota_ct.col(0) -= sumB1_cc * F_ct.col(0);
    }
    if( (Options_vec(8)==1) | (Options_vec(8)==2) ){
      // Project forward effect of F_ct from initial year through current year
      for( t=1; t<n_t; t++ ){
        iota_ct.col(t) = B1_cc * iota_ct.col(t-1) - F_ct.col(t);
      }
    }
  }else{
    iota_ct.setZero();
  }

  ////////////////////////
  // Random effects
  //  Spatial and spatio-temporal variation
  // (Should be first to optimize speedup involving normalization of GMRFs, but must be after interactions)
  ////////////////////////

  // Random field probability
  Eigen::SparseMatrix<Type> Q1( n_s, n_s );
  Eigen::SparseMatrix<Type> Q2( n_s, n_s );
  GMRF_t<Type> gmrf_Q;
  if( (Options_vec(7)==0) & (Options_vec(0)==0) ){
    Q1 = Q_spde(spde, exp(logkappa1));
    Q2 = Q_spde(spde, exp(logkappa2));
  }
  if( (Options_vec(7)==0) & (Options_vec(0)==1) ){
    Q1 = Q_spde(spde_aniso, exp(logkappa1), H);
    Q2 = Q_spde(spde_aniso, exp(logkappa2), H);
  }
  if( Options_vec(7)==1 ){
    Q1 = M0*pow(1+exp(logkappa1*2),2) + M1*(1+exp(logkappa1*2))*(-exp(logkappa1)) + M2*exp(logkappa1*2);
    Q2 = M0*pow(1+exp(logkappa2*2),2) + M1*(1+exp(logkappa2*2))*(-exp(logkappa2)) + M2*exp(logkappa2*2);
  }
  if( Options_vec(7)==2 ){
    Q1 = Q_network( logkappa1, n_s, parent_s, child_s, dist_s );
    Q2 = Q_network( logkappa2, n_s, parent_s, child_s, dist_s );
  }

  /////
  // 1st component
  /////
  gmrf_Q = GMRF( Q1, bool(Options(9)) );

  // Omega1
  array<Type> Omegamean1_sf(n_s, Omegainput1_sf.cols() );
  Omegamean1_sf.setZero();
  array<Type> Omega1_sc(n_s, n_c);
  Omega1_sc = gmrf_by_category_nll(FieldConfig(0,0), Options_vec(7), VamConfig(2), n_s, n_c, logkappa1, Omegainput1_sf, Omegamean1_sf, L_omega1_z, gmrf_Q, jnll_comp(0), this);

  // Projection for Omega1
  array<Type> Omega1_iz(n_i, c_iz.cols());
  array<Type> Omega1_gc(n_g, n_c);
  Omega1_iz.setZero();
  Omega1_gc.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
  for( int zc=0; zc<c_iz.cols(); zc++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
      Omega1_iz(i,zc) += Ais_x(Arow) * Omega1_sc(s,c_iz(i,zc));
    }
  }}
  Omega1_gc = project_knots( n_g, n_c, int(1), int(0), Omega1_sc, Ags_ij, Ags_x );

  // Epsilon1
  array<Type> Epsilonmean1_sf(n_s, n_f1 );
  // PDF for Epsilon1
  array<Type> Epsilon1_sct(n_s, n_c, n_t);
  for(t=0; t<n_t; t++){
    // PDF for B0 (not tied to autoregressive variation)
    if( (Options(11)==1) & (t==(Options(11)-1)) ){
      Epsilon1_sct.col(t) = gmrf_stationary_nll( Options_vec(7), n_s, n_c, logkappa1, Epsiloninput1_sft.col(t), covE1_cc, gmrf_Q, jnll_comp(1), this);
    }
    // PDF for first year of autoregression
    if( t==(Options(11)+0) ){
      Epsilonmean1_sf.setZero();
      Epsilon1_sct.col(t) = gmrf_by_category_nll(FieldConfig(1,0), Options_vec(7), VamConfig(2), n_s, n_c, logkappa1, Epsiloninput1_sft.col(t), Epsilonmean1_sf, L_epsilon1_z, gmrf_Q, jnll_comp(1), this);
    }
    // PDF for subsequent years of autoregression
    if( t>=(Options(11)+1) ){
      // Prediction for spatio-temporal component
      // Default, and also necessary whenever VamConfig(2)==1 & n_f1!=n_c
      if( (VamConfig(0)==0) | ((n_f1!=n_c) & (VamConfig(2)==1)) ){
        // If no interactions, then just autoregressive for factors
        for(s=0; s<n_s; s++){
        for(int f=0; f<n_f1; f++){
          Epsilonmean1_sf(s,f) = Epsilon_rho1_f(f) * Epsiloninput1_sft(s,f,t-1);
        }}
      }else{
        // Impact of interactions, B_ff
        Epsilonmean1_sf.setZero();
        for(s=0; s<n_s; s++){
        for(int f1=0; f1<n_f1; f1++){
        for(int f2=0; f2<n_f1; f2++){
          if( VamConfig(2)==0 ){
            Epsilonmean1_sf(s,f1) += B_ff(f1,f2) * Epsiloninput1_sft(s,f2,t-1);
            if( f1==f2 ) Epsilonmean1_sf(s,f1) += Epsilon_rho1_f(f1) * Epsiloninput1_sft(s,f2,t-1);
          }
          if( VamConfig(2)==1 ){
            Epsilonmean1_sf(s,f1) += B_ff(f1,f2) * Epsilon1_sct(s,f2,t-1);
            if( f1==f2 ) Epsilonmean1_sf(s,f1) += Epsilon_rho1_f(f1) * Epsilon1_sct(s,f2,t-1);
          }
        }}}
      }
      // Hyperdistribution for spatio-temporal component
      Epsilon1_sct.col(t) = gmrf_by_category_nll(FieldConfig(1,0), Options_vec(7), VamConfig(2), n_s, n_c, logkappa1, Epsiloninput1_sft.col(t), Epsilonmean1_sf, L_epsilon1_z, gmrf_Q, jnll_comp(1), this);
    }
  }

  // Projection for Epsilon1
  array<Type> Epsilon1_izz(n_i, c_iz.cols(), t_iz.cols());
  array<Type> Epsilon1_gct(n_g, n_c, n_t);
  Epsilon1_izz.setZero();
  Epsilon1_gct.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
  for( int zc=0; zc<c_iz.cols(); zc++ ){
  for( int zt=0; zt<t_iz.cols(); zt++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
    if( (t_iz(i,zt)>=0) & (t_iz(i,zt)<n_t) ){
      Epsilon1_izz(i,zc,zt) += Ais_x(Arow) * Epsilon1_sct(s,c_iz(i,zc),t_iz(i,zt));
    }}
  }}}
  Epsilon1_gct = project_knots( n_g, n_c, n_t, int(1), Epsilon1_sct, Ags_ij, Ags_x );

  // Xi1_scp
  array<Type> Ximean1_sc(n_s, 1);
  array<Type> Xi1_scp(n_s, n_c, n_p);
  vector<Type> Sigma1(1);
  array<Type> Tmp1_sc(n_s, 1);
  Ximean1_sc.setZero();
  Xi1_scp.setZero();
  for(p=0; p<n_p; p++){
  for(c=0; c<n_c; c++){
    // Hyperdistribution for spatially varying coefficients (uses IID option)
    if( (Xconfig_zcp(0,c,p)==2) | (Xconfig_zcp(0,c,p)==3) ){
      Sigma1(0) = sigmaXi1_cp(c,p);
      Tmp1_sc.col(0) = Xiinput1_scp.col(p).col(c);
      Xi1_scp.col(p).col(c) = gmrf_by_category_nll( int(-2), Options_vec(7), VamConfig(2), n_s, int(1), logkappa1, Tmp1_sc, Ximean1_sc, Sigma1, gmrf_Q, jnll_comp(14), this);
    }
  }}

  // Projection for Xi1
  array<Type> Xi1_izp(n_i, c_iz.cols(), n_p);
  array<Type> Xi1_gcp(n_g, n_c, n_p);
  Xi1_izp.setZero();
  Xi1_gcp.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
  for( int zc=0; zc<c_iz.cols(); zc++ ){
  for( p=0; p<n_p; p++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
      Xi1_izp(i,zc,p) += Ais_x(Arow) * Xi1_scp(s,c_iz(i,zc),p);
    }
  }}}
  Xi1_gcp = project_knots( n_g, n_c, n_p, int(1), Xi1_scp, Ags_ij, Ags_x );

  /////
  // 2nd component
  /////
  gmrf_Q = GMRF( Q2, bool(Options(9)) );

  // Omega2
  array<Type> Omegamean2_sf(n_s, Omegainput2_sf.cols() );
  Omegamean2_sf.setZero();
  array<Type> Omega2_sc(n_s, n_c);
  Omega2_sc = gmrf_by_category_nll(FieldConfig(0,1), Options_vec(7), VamConfig(2), n_s, n_c, logkappa2, Omegainput2_sf, Omegamean2_sf, L_omega2_z, gmrf_Q, jnll_comp(2), this);

  // Projection for Omega2
  array<Type> Omega2_iz(n_i, c_iz.cols());
  array<Type> Omega2_gc(n_g, n_c);
  Omega2_iz.setZero();
  Omega2_gc.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
  for( int zc=0; zc<c_iz.cols(); zc++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
      Omega2_iz(i,zc) += Ais_x(Arow) * Omega2_sc(s,c_iz(i,zc));
    }
  }}
  Omega2_gc = project_knots( n_g, n_c, int(1), int(0), Omega2_sc, Ags_ij, Ags_x );

  // Epsilon2
  array<Type> Epsilonmean2_sf(n_s, n_f2);
  // PDF for Epsilon2
  array<Type> Epsilon2_sct(n_s, n_c, n_t);
  for(t=0; t<n_t; t++){
    // PDF for B0 (not tied to autoregressive variation)
    if( (Options(11)==1) & (t==(Options(11)-1)) ){
      Epsilon2_sct.col(t) = gmrf_stationary_nll( Options_vec(7), n_s, n_c, logkappa2, Epsiloninput2_sft.col(t), covE2_cc, gmrf_Q, jnll_comp(3), this);
    }
    // PDF for first year of autoregression
    if( t==(Options(11)+0) ){
      Epsilonmean2_sf.setZero();
      Epsilon2_sct.col(t) = gmrf_by_category_nll(FieldConfig(1,1), Options_vec(7), VamConfig(2), n_s, n_c, logkappa2, Epsiloninput2_sft.col(t), Epsilonmean2_sf, L_epsilon2_z, gmrf_Q, jnll_comp(3), this);
    }
    // PDF for subsequent years of autoregression
    if( t>=(Options(11)+1) ){
      // Prediction for spatio-temporal component
      // Default, and also necessary whenever VamConfig(2)==1 & n_f2!=n_c
      if( (VamConfig(0)==0) | ((n_f2!=n_c) & (VamConfig(2)==1)) ){
        // If no interactions, then just autoregressive for factors
        for(s=0; s<n_s; s++){
        for(int f=0; f<n_f2; f++){
          Epsilonmean2_sf(s,f) = Epsilon_rho2_f(f) * Epsiloninput2_sft(s,f,t-1);
        }}
      }else{
        // Impact of interactions, B_ff
        Epsilonmean2_sf.setZero();
        for(s=0; s<n_s; s++){
        for(int f1=0; f1<n_f2; f1++){
        for(int f2=0; f2<n_f2; f2++){
          if( VamConfig(2)==0 ){
            Epsilonmean2_sf(s,f1) += B_ff(f1,f2) * Epsiloninput2_sft(s,f2,t-1);
            if( f1==f2 ) Epsilonmean2_sf(s,f1) += Epsilon_rho2_f(f1) * Epsiloninput2_sft(s,f2,t-1);
          }
          if( VamConfig(2)==1 ){
            Epsilonmean2_sf(s,f1) += B_ff(f1,f2) * Epsilon2_sct(s,f2,t-1);
            if( f1==f2 ) Epsilonmean2_sf(s,f1) += Epsilon_rho2_f(f1) * Epsilon2_sct(s,f2,t-1);
          }
        }}}
      }
      // Hyperdistribution for spatio-temporal component
      Epsilon2_sct.col(t) = gmrf_by_category_nll(FieldConfig(1,1), Options_vec(7), VamConfig(2), n_s, n_c, logkappa2, Epsiloninput2_sft.col(t), Epsilonmean2_sf, L_epsilon2_z, gmrf_Q, jnll_comp(3), this);
    }
  }

  // Projection for Epsilon2
  array<Type> Epsilon2_izz(n_i, c_iz.cols(), t_iz.cols());
  array<Type> Epsilon2_gct(n_g, n_c, n_t);
  Epsilon2_izz.setZero();
  Epsilon2_gct.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
  for( int zc=0; zc<c_iz.cols(); zc++ ){
  for( int zt=0; zt<t_iz.cols(); zt++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
    if( (t_iz(i,zt)>=0) & (t_iz(i,zt)<n_t) ){
      Epsilon2_izz(i,zc,zt) += Ais_x(Arow) * Epsilon2_sct(s,c_iz(i,zc),t_iz(i,zt));
    }}
  }}}
  Epsilon2_gct = project_knots( n_g, n_c, n_t, int(1), Epsilon2_sct, Ags_ij, Ags_x );

  // Xi2_scp
  array<Type> Ximean2_sc(n_s, 1);
  array<Type> Xi2_scp(n_s, n_c, n_p);
  vector<Type> Sigma2(1);
  array<Type> Tmp2_sc(n_s, 1);
  Ximean2_sc.setZero();
  Xi2_scp.setZero();
  for(p=0; p<n_p; p++){
  for(c=0; c<n_c; c++){
    // Hyperdistribution for spatially varying coefficients (uses IID option)
    if( (Xconfig_zcp(1,c,p)==2) | (Xconfig_zcp(1,c,p)==3) ){
      Tmp2_sc.col(0) = Xiinput2_scp.col(p).col(c);
      Sigma2(0) = sigmaXi2_cp(c,p);
      Xi2_scp.col(p).col(c) = gmrf_by_category_nll( int(-2), Options_vec(7), VamConfig(2), n_s, int(1), logkappa2, Tmp2_sc, Ximean2_sc, Sigma2, gmrf_Q, jnll_comp(15), this);
    }
  }}

  // Projection for Xi2
  array<Type> Xi2_izp(n_i, c_iz.cols(), n_p);
  array<Type> Xi2_gcp(n_g, n_c, n_p);
  Xi2_izp.setZero();
  Xi2_gcp.setZero();
  for( int Arow=0; Arow<Ais_ij.rows(); Arow++ ){
  for( int zc=0; zc<c_iz.cols(); zc++ ){
  for( p=0; p<n_p; p++ ){
    i = Ais_ij(Arow,0);
    s = Ais_ij(Arow,1);
    if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
      Xi2_izp(i,zc,p) += Ais_x(Arow) * Xi2_scp(s,c_iz(i,zc),p);
    }
  }}}
  Xi2_gcp = project_knots( n_g, n_c, n_p, int(1), Xi2_scp, Ags_ij, Ags_x );

  // Normalization of GMRFs to normalize during outer-optimization step in R
  Type jnll_GMRF = jnll_comp(0) + jnll_comp(1) + jnll_comp(2) + jnll_comp(3);
  if( include_data == 0 ){
    return( jnll_GMRF );
  }

  ////////////////////////
  // Random effects
  //  1. Correlated overdispersion
  //  2. Intercepts
  //  3. Lognormal-Poisson overdispersion
  ////////////////////////

  ////// Probability of correlated overdispersion among bins
  // 1st component
  int n_eta_f1;
  n_eta_f1 = eta1_vf.cols();
  matrix<Type> eta1_mean_vf(n_v, n_eta_f1);
  eta1_mean_vf.setZero();
  matrix<Type> eta1_vc(n_v, n_c);
  eta1_vc = covariation_by_category_nll( OverdispersionConfig(0), n_v, n_c, eta1_vf, eta1_mean_vf, L1_z, jnll_comp(4), this );
  // 1st component
  int n_eta_f2;
  n_eta_f2 = eta2_vf.cols();
  matrix<Type> eta2_mean_vf(n_v, n_eta_f2);
  eta2_mean_vf.setZero();
  matrix<Type> eta2_vc(n_v, n_c);
  eta2_vc = covariation_by_category_nll( OverdispersionConfig(1), n_v, n_c, eta2_vf, eta2_mean_vf, L2_z, jnll_comp(5), this );

  ////// Probability of correlated innovations on intercepts
  // 1st component
  Type jnll_beta1 = 0;
  int n_beta_f1;
  n_beta_f1 = beta1_ft.rows();
  matrix<Type> beta1_mean_tf(n_t, n_beta_f1);
  matrix<Type> beta1_tf( n_t, n_beta_f1 );
  beta1_tf = beta1_ft.transpose();
  for( int f=0; f<n_beta_f1; f++ ){
    beta1_mean_tf(0,f) = Type(0.0);
    for( t=1; t<n_t; t++ ){
      beta1_mean_tf(t,f) = beta1_tf(t-1,f) * Beta_rho1_f(f);
    }
  }
  matrix<Type> beta1_tc(n_t, n_c);
  beta1_tc = covariation_by_category_nll( FieldConfig(2,0), n_t, n_c, beta1_tf, beta1_mean_tf, L_beta1_z, jnll_beta1, this );
  for( c=0; c<n_c; c++ ){
  for( t=0; t<n_t; t++ ){
    beta1_tc(t,c) += Beta_mean1_c(c);
  }}
  if( (RhoConfig(0)==1) | (RhoConfig(0)==2) | (RhoConfig(0)==4) ){
    jnll_comp(8) = jnll_beta1;
  }

  // 2nd component
  Type jnll_beta2 = 0;
  int n_beta_f2;
  n_beta_f2 = beta2_ft.rows();
  matrix<Type> beta2_mean_tf(n_t, n_beta_f2);
  matrix<Type> beta2_tf( n_t, n_beta_f2 );
  beta2_tf = beta2_ft.transpose();
  for( int f=0; f<n_beta_f2; f++ ){
    beta2_mean_tf(0,f) = Type(0.0);
    for( t=1; t<n_t; t++ ){
      beta2_mean_tf(t,f) = beta2_tf(t-1,f) * Beta_rho2_f(f);
    }
  }
  matrix<Type> beta2_tc(n_t, n_c);
  beta2_tc = covariation_by_category_nll( FieldConfig(2,1), n_t, n_c, beta2_tf, beta2_mean_tf, L_beta2_z, jnll_beta2, this );
  for( c=0; c<n_c; c++ ){
  for( t=0; t<n_t; t++ ){
    beta2_tc(t,c) += Beta_mean2_c(c);
  }}
  if( (RhoConfig(1)==1) | (RhoConfig(1)==2) | (RhoConfig(1)==4) | (RhoConfig(1)==6) ){
    jnll_comp(9) = jnll_beta2;
  }

  ////// Penalty on lognormal-Poisson overdispesrion delta_i
  for(i=0; i<delta_i.size(); i++){
    if( (ObsModel_ez(e_i(i),0)==11) | (ObsModel_ez(e_i(i),0)==14) ){
      jnll_comp(12) -= dnorm( delta_i(i), Type(0.0), Type(1.0), true );
      // Simulate new values when using obj.simulate()
      SIMULATE{
        delta_i(i) = rnorm( Type(0.0), Type(1.0) );
      }
    }
  }

  ////////////////////////
  // Covariate effects
  ////////////////////////

  vector<Type> zeta1_i = Q_ik * lambda1_k.matrix();
  vector<Type> zeta2_i = Q_ik * lambda2_k.matrix();
  array<Type> eta1_izz(n_i, c_iz.cols(), t_iz.cols());
  array<Type> eta2_izz(n_i, c_iz.cols(), t_iz.cols());
  array<Type> eta1_gct(n_g, n_c, n_t);
  array<Type> eta2_gct(n_g, n_c, n_t);
  eta1_izz.setZero();
  eta2_izz.setZero();
  eta1_gct.setZero();
  eta2_gct.setZero();
  for(p=0; p<n_p; p++){
    for(c=0; c<n_c; c++){
    for(t=0; t<n_t; t++){
    for(g=0; g<n_g; g++){
      eta1_gct(g,c,t) += (gamma1_ctp(c,t,p) + Xi1_gcp(g,c,p)) * X_gtp(g,t,p);
      eta2_gct(g,c,t) += (gamma2_ctp(c,t,p) + Xi2_gcp(g,c,p)) * X_gtp(g,t,p);
    }}}
    for( i=0; i<n_i; i++ ){
    for( int zc=0; zc<c_iz.cols(); zc++ ){
    for( int zt=0; zt<t_iz.cols(); zt++ ){
      if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
      if( (t_iz(i,zt)>=0) & (t_iz(i,zt)<n_t) ){
        eta1_izz(i,zc,zt) += (gamma1_ctp(c_iz(i,zc),t_iz(i,zt),p) + Xi1_izp(i,zc,p)) * X_itp(i,t_iz(i,zt),p);
        eta2_izz(i,zc,zt) += (gamma2_ctp(c_iz(i,zc),t_iz(i,zt),p) + Xi2_izp(i,zc,p)) * X_itp(i,t_iz(i,zt),p);
      }}
    }}}
  }

  ////////////////////////
  // Likelihood for data
  ////////////////////////

  // Derived quantities
  vector<Type> var_i(n_i);
  Type tmp_calc1;
  Type tmp_calc2;
  Type log_tmp_calc2;
  // Linear predictor (pre-link) for presence/absence component
  matrix<Type> P1_iz(n_i,c_iz.cols());
  // Response predictor (post-link)
  // ObsModel_ez(e,0) = 0:4 or 11:12: probability ("phi") that data is greater than zero
  // ObsModel_ez(e,0) = 5 (ZINB):  phi = 1-ZeroInflation_prob -> Pr[D=0] = NB(0|mu,var)*phi + (1-phi) -> Pr[D>0] = phi - NB(0|mu,var)*phi
  vector<Type> R1_i(n_i);
  vector<Type> log_one_minus_R1_i(n_i);
  vector<Type> log_R1_i(n_i);
  vector<Type> LogProb1_i(n_i);
  // Linear predictor (pre-link) for positive component
  matrix<Type> P2_iz(n_i,c_iz.cols());
  // Response predictor (post-link)
  // ObsModel_ez(e,0) = 0:3, 11:12:  expected value of data, given that data is greater than zero -> E[D] = mu*phi
  // ObsModel_ez(e,0) = 4 (ZANB):  expected value ("mu") of neg-bin PRIOR to truncating Pr[D=0] -> E[D] = mu/(1-NB(0|mu,var))*phi  ALSO  Pr[D] = NB(D|mu,var)/(1-NB(0|mu,var))*phi
  // ObsModel_ez(e,0) = 5 (ZINB):  expected value of data for non-zero-inflation component -> E[D] = mu*phi
  vector<Type> R2_i(n_i);
  vector<Type> log_R2_i(n_i);
  vector<Type> LogProb2_i(n_i);
  vector<Type> maxJ_i(n_i);
  vector<Type> diag_z(4);
  matrix<Type> diag_iz(n_i,4);
  diag_iz.setZero();  // Used to track diagnostics for Tweedie distribution (columns: 0=maxJ; 1=maxW; 2=lowerW; 3=upperW)
  P1_iz.setZero();
  P2_iz.setZero();

  // Likelihood contribution from observations
  LogProb1_i.setZero();
  LogProb2_i.setZero();
  for(i=0; i<n_i; i++){
    if( !isNA(b_i(i)) ){
      // Linear predictors
      for( int zc=0; zc<c_iz.cols(); zc++ ){
        if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
          P1_iz(i,zc) = Omega1_iz(i,zc) + zeta1_i(i) + eta1_vc(v_i(i),c_iz(i,zc));
          P2_iz(i,zc) = Omega2_iz(i,zc) + zeta2_i(i) + eta2_vc(v_i(i),c_iz(i,zc));
          for( int zt=0; zt<t_iz.cols(); zt++ ){
            if( (t_iz(i,zt)>=0) & (t_iz(i,zt)<n_t) ){  // isNA doesn't seem to work for IMATRIX type
              P1_iz(i,zc) += beta1_tc(t_iz(i,zt),c_iz(i,zc)) + Epsilon1_izz(i,zc,zt)*exp(log_sigmaratio1_z(zt)) + eta1_izz(i,zc,zt) + iota_ct(c_iz(i,zc),t_iz(i,zt));
              P2_iz(i,zc) += beta2_tc(t_iz(i,zt),c_iz(i,zc)) + Epsilon2_izz(i,zc,zt)*exp(log_sigmaratio2_z(zt)) + eta2_izz(i,zc,zt);
            }
          }
        }
      }
      // Apply link function to calculate responses
      if( (ObsModel_ez(c_iz(i,0),1)==0) | (ObsModel_ez(c_iz(i,0),1)==3) ){
        // Log and logit-link, where area-swept only affects positive catch rate exp(P2_i(i))
        // P1_i: Logit-Probability of occurrence;  R1_i:  Probability of occurrence
        // P2_i: Log-Positive density prediction;  R2_i:  Positive density prediction
        R1_i(i) = invlogit( P1_iz(i,0) );
        R2_i(i) = a_i(i) * exp( P2_iz(i,0) );
        // Calulate in logspace to prevent numerical over/under-flow
        log_R1_i(i) = log(Type(1.0)) - logspace_add( log(Type(1.0)), -1.0*P1_iz(i,0) );
        log_one_minus_R1_i(i) = log(Type(1.0)) - logspace_add( log(Type(1.0)), P1_iz(i,0) );
        log_R2_i(i) = log(a_i(i)) + P2_iz(i,0);
      }
      if( (ObsModel_ez(c_iz(i,0),1)==1) | (ObsModel_ez(c_iz(i,0),1)==4) ){
        // Poisson-process link, where area-swept affects numbers density exp(P1_i(i))
        // P1_i: Log-numbers density;  R1_i:  Probability of occurrence
        // P2_i: Log-average weight;  R2_i:  Positive density prediction
        tmp_calc1 = 0;
        tmp_calc2 = 0;
        log_tmp_calc2 = 0;
        for( int zc=0; zc<c_iz.cols(); zc++ ){
          if( (c_iz(i,zc)>=0) & (c_iz(i,zc)<n_c) ){
            tmp_calc1 += exp(P1_iz(i,zc));
            tmp_calc2 += exp(P1_iz(i,zc)) * exp(P2_iz(i,zc));
            if( zc==0 ) log_tmp_calc2 = P1_iz(i,zc) + P2_iz(i,zc);
            if( zc>=1 ) log_tmp_calc2 = logspace_add( log_tmp_calc2, P1_iz(i,zc) + P2_iz(i,zc) );
          }
        }
        R1_i(i) = Type(1.0) - exp( -1*a_i(i)*tmp_calc1 );
        R2_i(i) = a_i(i) * tmp_calc2 / R1_i(i);
        // Calulate in logspace to prevent numerical over/under-flow
        log_R1_i(i) = logspace_sub( log(Type(1.0)), -1*a_i(i)*tmp_calc1 );
        log_one_minus_R1_i(i) = -1*a_i(i)*tmp_calc1;
        log_R2_i(i) = log(a_i(i)) + log_tmp_calc2 - log_R1_i(i);
      }
      if( ObsModel_ez(c_iz(i,0),1)==2 ){
        // Tweedie link, where area-swept affects numbers density exp(P1_i(i))
        // P1_i: Log-numbers density;  R1_i:  Expected numbers
        // P2_i: Log-average weight;  R2_i:  Expected average weight
        R1_i(i) = a_i(i) * exp( P1_iz(i,0) );
        R2_i(i) = exp( P2_iz(i,0) );
        // Calulate in logspace to prevent numerical over/under-flow
        log_R1_i(i) = log(a_i(i)) + P1_iz(i,0);
        log_R2_i(i) = P2_iz(i,0);
      }
      // Likelihood for delta-models with continuous positive support
      if( (ObsModel_ez(e_i(i),0)==0) | (ObsModel_ez(e_i(i),0)==1) | (ObsModel_ez(e_i(i),0)==2) ){
        // Presence-absence likelihood
        if( (ObsModel_ez(e_i(i),1)==0) | (ObsModel_ez(e_i(i),1)==1) | (ObsModel_ez(e_i(i),1)==3) | (ObsModel_ez(e_i(i),1)==4) ){
          if( b_i(i) > 0 ){
            LogProb1_i(i) = log_R1_i(i);
          }else{
            LogProb1_i(i) = log_one_minus_R1_i(i);
          }
        }else{
          if( b_i(i) > 0 ){
            LogProb1_i(i) = log( R1_i(i) );
          }else{
            LogProb1_i(i) = log( 1-R1_i(i) );
          }
        }
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = rbinom( Type(1), R1_i(i) );
        }
        // Positive density likelihood -- models with continuous positive support
        if( b_i(i) > 0 ){    // 1e-500 causes overflow on laptop
          if(ObsModel_ez(e_i(i),0)==0){
            LogProb2_i(i) = dnorm(b_i(i), R2_i(i), SigmaM(e_i(i),0), true);
            // Simulate new values when using obj.simulate()
            SIMULATE{
              b_i(i) = rnorm( R2_i(i), SigmaM(e_i(i),0) );
            }
          }
          if(ObsModel_ez(e_i(i),0)==1){
            LogProb2_i(i) = dlnorm(b_i(i), log_R2_i(i)-pow(SigmaM(e_i(i),0),2)/2, SigmaM(e_i(i),0), true); // log-space
            // Simulate new values when using obj.simulate()
            SIMULATE{
              b_i(i) = exp(rnorm( log_R2_i(i)-pow(SigmaM(e_i(i),0),2)/2, SigmaM(e_i(i),0) ));
            }
          }
          if(ObsModel_ez(e_i(i),0)==2){
            LogProb2_i(i) = dgamma(b_i(i), 1/pow(SigmaM(e_i(i),0),2), R2_i(i)*pow(SigmaM(e_i(i),0),2), true); // shape = 1/CV^2, scale = mean*CV^2
            // Simulate new values when using obj.simulate()
            SIMULATE{
              b_i(i) = rgamma( 1/pow(SigmaM(e_i(i),0),2), R2_i(i)*pow(SigmaM(e_i(i),0),2) );
            }
          }
        }else{
          LogProb2_i(i) = 0;
        }
      }
      // Likelihood for Tweedie model with continuous positive support
      if(ObsModel_ez(e_i(i),0)==8){
        LogProb1_i(i) = 0;
        //dPoisGam( Type x, Type shape, Type scale, Type intensity, Type &max_log_w_j, int maxsum=50, int minsum=1, give_log=0 )
        LogProb2_i(i) = dPoisGam( b_i(i), SigmaM(e_i(i),0), R2_i(i), R1_i(i), diag_z, Options_vec(5), Options_vec(6), true );
        diag_iz.row(i) = diag_z;
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = 0;   // Option not available
        }
      }
      // Likelihood #2 for Tweedie model with continuous positive support
      if(ObsModel_ez(e_i(i),0)==10){
        // Packaged code
        LogProb1_i(i) = 0;
        // dtweedie( Type y, Type mu, Type phi, Type p, int give_log=0 )
        // R1*R2 = mean
        LogProb2_i(i) = dtweedie( b_i(i), R1_i(i)*R2_i(i), R1_i(i), invlogit(SigmaM(e_i(i),0))+1.0, true );
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = 0;   // Option not available
        }
      }
      ///// Likelihood for models with discrete support
      // Zero-inflated negative binomial (not numerically stable!)
      if(ObsModel_ez(e_i(i),0)==5){
        var_i(i) = R2_i(i)*(1.0+SigmaM(e_i(i),0)) + pow(R2_i(i),2.0)*SigmaM(c_iz(i,0),1);
        if( b_i(i)==0 ){
          //LogProb2_i(i) = log( (1-R1_i(i)) + dnbinom2(Type(0.0), R2_i(i), var_i(i), false)*R1_i(i) ); //  Pr[X=0] = 1-phi + NB(X=0)*phi
          LogProb2_i(i) = logspace_add( log(1-R1_i(i)), dnbinom2(Type(0.0),R2_i(i),var_i(i),true)+log(R1_i(i)) ); //  Pr[X=0] = 1-phi + NB(X=0)*phi
        }else{
          LogProb2_i(i) = dnbinom2(b_i(i), R2_i(i), var_i(i), true) + log(R1_i(i)); // Pr[X=x] = NB(X=x)*phi
        }
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = rbinom( Type(1), R1_i(i) );
          if( b_i(i)>0 ){
            b_i(i) = rnbinom2( R2_i(i), var_i(i) );
          }
        }
      }
      // Conway-Maxwell-Poisson
      if(ObsModel_ez(e_i(i),0)==6){
        LogProb2_i(i) = dCMP(b_i(i), R2_i(i), exp(P1_iz(i,0)), true, Options_vec(5));
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = 0;   // Option not available
        }
      }
      // Zero-inflated Poisson
      if(ObsModel_ez(e_i(i),0)==7){
        if( b_i(i)==0 ){
          //LogProb2_i(i) = log( (1-R1_i(i)) + dpois(Type(0.0), R2_i(i), false)*R1_i(i) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
          LogProb2_i(i) = logspace_add( log(1-R1_i(i)), dpois(Type(0.0),R2_i(i),true)+log(R1_i(i)) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
        }else{
          LogProb2_i(i) = dpois(b_i(i), R2_i(i), true) + log(R1_i(i)); // Pr[X=x] = Pois(X=x)*phi
        }
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = rbinom( Type(1), R1_i(i) );
          if( b_i(i)>0 ){
            b_i(i) = rpois( R2_i(i) );
          }
        }
      }
      // Binned Poisson (for REEF data: 0=none; 1=1; 2=2-10; 3=>11)
        /// Doesn't appear stable given spatial or spatio-temporal variation
      if(ObsModel_ez(e_i(i),0)==9){
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
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = 0;   // Option not available
        }
      }
      // Zero-inflated Lognormal Poisson
      if(ObsModel_ez(e_i(i),0)==11){
        if( b_i(i)==0 ){
          //LogProb2_i(i) = log( (1-R1_i(i)) + dpois(Type(0.0), R2_i(i), false)*R1_i(i) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
          LogProb2_i(i) = logspace_add( log(1-R1_i(i)), dpois(Type(0.0),R2_i(i)*exp(SigmaM(e_i(i),0)*delta_i(i)-0.5*pow(SigmaM(e_i(i),0),2)),true)+log(R1_i(i)) ); //  Pr[X=0] = 1-phi + Pois(X=0)*phi
        }else{
          LogProb2_i(i) = dpois(b_i(i), R2_i(i)*exp(SigmaM(e_i(i),0)*delta_i(i)-0.5*pow(SigmaM(e_i(i),0),2)), true) + log(R1_i(i)); // Pr[X=x] = Pois(X=x)*phi
        }
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = rbinom( Type(1), R1_i(i) );
          if( b_i(i)>0 ){
            b_i(i) = rpois( R2_i(i)*exp(SigmaM(e_i(i),0)*delta_i(i)-0.5*pow(SigmaM(e_i(i),0),2)) );
          }
        }
      }
      // Non-zero-inflated Poisson using log link from 1st linear predictor
      if(ObsModel_ez(e_i(i),0)==12){
        LogProb2_i(i) = dpois(b_i(i), R1_i(i), true);
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = rpois( R1_i(i) );
        }
      }
      // Non-zero-inflated Bernoulli using cloglog link from 1st lilnear predict
      if(ObsModel_ez(e_i(i),0)==13){
        if( b_i(i)==0 ){
          LogProb2_i(i) = dpois(Type(0), R1_i(i), true);
        }else{
          LogProb2_i(i) = logspace_sub( log(Type(1.0)), dpois(Type(0), R1_i(i), true) );
        }
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = rpois( R1_i(i) );
          if( b_i(i)>0 ){
            b_i(i) = 1;
          }
        }
      }
      // Non-zero-inflated Lognormal-Poisson using log link from 1st linear predictor
      if(ObsModel_ez(e_i(i),0)==14){
        LogProb2_i(i) = dpois(b_i(i), R1_i(i)*exp(SigmaM(e_i(i),0)*delta_i(i)-0.5*pow(SigmaM(e_i(i),0),2)), true);
        // Simulate new values when using obj.simulate()
        SIMULATE{
          b_i(i) = rpois( R1_i(i)*exp(SigmaM(e_i(i),0)*delta_i(i)-0.5*pow(SigmaM(e_i(i),0),2)) );
        }
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
  // Calculate index of abundance and density
  ////////////////////////

  // Number of output-years
  int n_y = t_yz.rows();

  // Predictive distribution -- ObsModel_ez(e,0)==4 isn't implemented (it had a bug previously)
  Type a_average = a_i.sum()/a_i.size();
  array<Type> P1_gcy(n_g, n_c, n_y);
  array<Type> R1_gcy(n_g, n_c, n_y);
  array<Type> P2_gcy(n_g, n_c, n_y);
  array<Type> R2_gcy(n_g, n_c, n_y);
  array<Type> D_gcy(n_g, n_c, n_y);
  for(c=0; c<n_c; c++){
  for(int y=0; y<n_y; y++){
  for(g=0; g<n_g; g++){
    // Calculate linear predictors
    P1_gcy(g,c,y) = Omega1_gc(g,c);
    P2_gcy(g,c,y) =  Omega2_gc(g,c);
    for( int z=0; z<t_yz.cols(); z++ ){
      if( (t_yz(y,z)>=0) & (t_yz(y,z)<n_t) ){    // isNA doesn't seem to work for IMATRIX type
        P1_gcy(g,c,y) += beta1_tc(t_yz(y,z),c) + Epsilon1_gct(g,c,t_yz(y,z))*exp(log_sigmaratio1_z(z)) + eta1_gct(g,c,t_yz(y,z)) + iota_ct(c,t_yz(y,z));
        P2_gcy(g,c,y) += beta2_tc(t_yz(y,z),c) + Epsilon2_gct(g,c,t_yz(y,z))*exp(log_sigmaratio2_z(z)) + eta2_gct(g,c,t_yz(y,z));
      }
    }
    // Calculate predictors in link-space
    if( (ObsModel_ez(c,1)==0) | (ObsModel_ez(c,1)==3) ){
      R1_gcy(g,c,y) = invlogit( P1_gcy(g,c,y) );
      R2_gcy(g,c,y) = exp( P2_gcy(g,c,y) );
      D_gcy(g,c,y) = R1_gcy(g,c,y) * R2_gcy(g,c,y);
    }
    if( (ObsModel_ez(c,1)==1) | (ObsModel_ez(c,1)==4) ){
      R1_gcy(g,c,y) = Type(1.0) - exp( -exp(P1_gcy(g,c,y)) );
      R2_gcy(g,c,y) = exp(P1_gcy(g,c,y)) / R1_gcy(g,c,y) * exp( P2_gcy(g,c,y) );
      D_gcy(g,c,y) = exp( P1_gcy(g,c,y) ) * exp( P2_gcy(g,c,y) );        // Use this line to prevent numerical over/underflow
    }
    if( ObsModel_ez(c,1)==2 ){
      R1_gcy(g,c,y) = exp( P1_gcy(g,c,y) );
      R2_gcy(g,c,y) = exp( P2_gcy(g,c,y) );
      D_gcy(g,c,y) = R1_gcy(g,c,y) * R2_gcy(g,c,y);
    }
  }}}

  // Calculate indices
  array<Type> Index_gcyl(n_g, n_c, n_y, n_l);
  array<Type> Index_cyl(n_c, n_y, n_l);
  array<Type> ln_Index_cyl(n_c, n_y, n_l);
  Index_cyl.setZero();
  for(int y=0; y<n_y; y++){
  for(int l=0; l<n_l; l++){
    // Expand by area and convert from kg to metric tonnes
    for(c=0; c<n_c; c++){
      if( Expansion_cz(c,0)==0 ){
        for(g=0; g<n_g; g++){
          Index_gcyl(g,c,y,l) = D_gcy(g,c,y) * a_gl(g,l) / 1000;
          Index_cyl(c,y,l) += Index_gcyl(g,c,y,l);
        }
      }
    }
    // Expand by biomass for another category and convert from kg to metric tonnes
    for(c=0; c<n_c; c++){
      if( Expansion_cz(c,0)==1 ){
        for(g=0; g<n_g; g++){
          Index_gcyl(g,c,y,l) = D_gcy(g,c,y) * Index_gcyl(g,Expansion_cz(c,1),y,l);    // Had Index_gcyl(g,Expansion_cz(c,1)+1,y,l) in original draft but I can't remember why
          Index_cyl(c,y,l) += Index_gcyl(g,c,y,l);
        }
      }
    }
  }}
  ln_Index_cyl = log( Index_cyl );

  ////////////////////////
  // Calculate optional derived quantities
  ////////////////////////

  // Calculate B / B0
  if( Options(11)==1 ){
    array<Type> Bratio_cyl(n_c, n_y, n_l);
    array<Type> ln_Bratio_cyl(n_c, n_y, n_l);
    for(c=0; c<n_c; c++){
    for(int y=0; y<n_y; y++){
    for(int l=0; l<n_l; l++){
      Bratio_cyl(c,y,l) = Index_cyl(c,y,l) / Index_cyl(c,0,l);
    }}}
    ln_Bratio_cyl = log( Bratio_cyl );
    REPORT( Bratio_cyl );
    ADREPORT( Bratio_cyl );
    ADREPORT( ln_Bratio_cyl );
  }

  // Calculate other derived summaries
  // Each is the weighted-average X_xl over polygons (x) with weights equal to abundance in each polygon and time (where abundance is from the first index)
  array<Type> mean_Z_cym(n_c, n_y, n_m);
  if( Options(2)==1 ){
    mean_Z_cym.setZero();
    int report_summary_TF = false;
    for(c=0; c<n_c; c++){
    for(int y=0; y<n_y; y++){
    for(int m=0; m<n_m; m++){
      for(g=0; g<n_g; g++){
        if( Z_gm(g,m)!=0 ){
          mean_Z_cym(c,y,m) += Z_gm(g,m) * Index_gcyl(g,c,y,0)/Index_cyl(c,y,0);
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
    for(c=0; c<n_c; c++){
    for(int y=0; y<n_y; y++){
    for(int l=0; l<n_l; l++){
      for(g=0; g<n_g; g++){
        mean_D_cyl(c,y,l) += D_gcy(g,c,y) * Index_gcyl(g,c,y,l)/Index_cyl(c,y,l);
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
    if( FieldConfig(0,0)>0 ){
      matrix<Type> L1_omega_cf = loadings_matrix( L_omega1_z, n_c, FieldConfig(0,0) );
      matrix<Type> lowercov_uppercor_omega1 = L1_omega_cf * L1_omega_cf.transpose();
      lowercov_uppercor_omega1 = convert_upper_cov_to_cor( lowercov_uppercor_omega1 );
      REPORT( lowercov_uppercor_omega1 );
      ADREPORT( lowercov_uppercor_omega1 );
    }
    if( FieldConfig(1,0)>0 ){
      matrix<Type> L1_epsilon_cf = loadings_matrix( L_epsilon1_z, n_c, FieldConfig(1,0) );
      matrix<Type> lowercov_uppercor_epsilon1 = L1_epsilon_cf * L1_epsilon_cf.transpose();
      lowercov_uppercor_epsilon1 = convert_upper_cov_to_cor( lowercov_uppercor_epsilon1 );
      REPORT( lowercov_uppercor_epsilon1 );
      ADREPORT( lowercov_uppercor_epsilon1 );
    }
    if( FieldConfig(2,0)>0 ){
      matrix<Type> L1_beta_cf = loadings_matrix( L_beta1_z, n_c, FieldConfig(2,0) );
      matrix<Type> lowercov_uppercor_beta1 = L1_beta_cf * L1_beta_cf.transpose();
      lowercov_uppercor_beta1 = convert_upper_cov_to_cor( lowercov_uppercor_beta1 );
      REPORT( lowercov_uppercor_beta1 );
      ADREPORT( lowercov_uppercor_beta1 );
    }
    if( FieldConfig(0,1)>0 ){
      matrix<Type> L2_omega_cf = loadings_matrix( L_omega2_z, n_c, FieldConfig(0,1) );
      matrix<Type> lowercov_uppercor_omega2 = L2_omega_cf * L2_omega_cf.transpose();
      lowercov_uppercor_omega2 = convert_upper_cov_to_cor( lowercov_uppercor_omega2 );
      REPORT( lowercov_uppercor_omega2 );
      ADREPORT( lowercov_uppercor_omega2 );
    }
    if( FieldConfig(1,1)>0 ){
      matrix<Type> L2_epsilon_cf = loadings_matrix( L_epsilon2_z, n_c, FieldConfig(1,1) );
      matrix<Type> lowercov_uppercor_epsilon2 = L2_epsilon_cf * L2_epsilon_cf.transpose();
      lowercov_uppercor_epsilon2 = convert_upper_cov_to_cor( lowercov_uppercor_epsilon2 );
      REPORT( lowercov_uppercor_epsilon2 );
      ADREPORT( lowercov_uppercor_epsilon2 );
    }
    if( FieldConfig(2,1)>0 ){
      matrix<Type> L1_beta_cf = loadings_matrix( L_beta2_z, n_c, FieldConfig(2,1) );
      matrix<Type> lowercov_uppercor_beta2 = L2_beta_cf * L2_beta_cf.transpose();
      lowercov_uppercor_beta2 = convert_upper_cov_to_cor( lowercov_uppercor_beta2 );
      REPORT( lowercov_uppercor_beta2 );
      ADREPORT( lowercov_uppercor_beta2 );
    }
  }

  // Synchrony
  if( Options(6)==1 ){
    int n_z = yearbounds_zz.rows();
    // Density ("D") or area-expanded total biomass ("B") for each category (use B when summing across sites)
    matrix<Type> D_gy( n_g, n_y );
    matrix<Type> B_cy( n_c, n_y );
    vector<Type> B_y( n_y );
    D_gy.setZero();
    B_cy.setZero();
    B_y.setZero();
    // Sample variance in category-specific density ("D") and biomass ("B")
    array<Type> varD_gcz( n_g, n_c, n_z );
    array<Type> varD_gz( n_g, n_z );
    array<Type> varB_cz( n_c, n_z );
    vector<Type> varB_z( n_z );
    vector<Type> varB_gbar_z( n_z );
    vector<Type> varB_cbar_z( n_z );
    vector<Type> ln_varB_z( n_z );
    vector<Type> ln_varB_gbar_z( n_z );
    vector<Type> ln_varB_cbar_z( n_z );
    array<Type> maxsdD_gz( n_g, n_z );
    array<Type> maxsdB_cz( n_c, n_z );
    vector<Type> maxsdB_z( n_z );
    varD_gcz.setZero();
    varD_gz.setZero();
    varB_cz.setZero();
    varB_z.setZero();
    varB_gbar_z.setZero();
    varB_cbar_z.setZero();
    maxsdD_gz.setZero();
    maxsdB_cz.setZero();
    maxsdB_z.setZero();
    // Proportion of total biomass ("P") for each location or each category
    matrix<Type> propB_gz( n_g, n_z );
    matrix<Type> propB_cz( n_c, n_z );
    propB_gz.setZero();
    propB_cz.setZero();
    // Synchrony indices
    matrix<Type> phi_gz( n_g, n_z );
    matrix<Type> phi_cz( n_c, n_z );
    vector<Type> phi_gbar_z( n_z );
    vector<Type> phi_cbar_z( n_z );
    vector<Type> phi_z( n_z );
    phi_gbar_z.setZero();
    phi_cbar_z.setZero();
    phi_z.setZero();
    // Calculate total biomass for different categories
    for( int y=0; y<n_y; y++ ){
      for( c=0; c<n_c; c++ ){
        for( g=0; g<n_g; g++ ){
          D_gy(g,y) += D_gcy(g,c,y);
          B_cy(c,y) += D_gcy(g,c,y) * a_gl(g,0);
          B_y(y) += D_gcy(g,c,y) * a_gl(g,0);
        }
      }
    }
    // Loop through periods (only using operations available in TMB, i.e., var, mean, row, col, segment)
    Type temp_mean;
    for( int z=0; z<n_z; z++ ){
      for( g=0; g<n_g; g++ ){
        // Variance for biomass in each category, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
        for( c=0; c<n_c; c++ ){
          temp_mean = 0;
          for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) temp_mean += D_gcy(g,c,y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
          for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ){
            varD_gcz(g,c,z) += pow(D_gcy(g,c,y)-temp_mean,2) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0));
          }
        }
        // Variance for combined biomass across categories, use sum(diff^2)/(length(diff)-1) where -1 in denominator is the sample-variance Bessel correction
        temp_mean = 0;
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) temp_mean += D_gy(g,y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
          varD_gz(g,z) += pow(D_gy(g,y)-temp_mean,2) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0));
        }
      }
      for( c=0; c<n_c; c++ ){
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
      for( g=0; g<n_g; g++ ){
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
          propB_gz(g,z) += a_gl(g,0) * D_gy(g,y) / B_y(y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
        }
      }
      // Proportion in each category
      for( c=0; c<n_c; c++ ){
        for( int y=yearbounds_zz(z,0); y<=yearbounds_zz(z,1); y++ ) {
          propB_cz(c,z) += B_cy(c,y) / B_y(y) / float(yearbounds_zz(z,1)-yearbounds_zz(z,0)+1);
        }
      }
      // Species-buffering index (calculate in Density so that areas with zero area are OK)
      for( g=0; g<n_g; g++ ){
        for( c=0; c<n_c; c++ ){
          maxsdD_gz(g,z) += pow(varD_gcz(g,c,z), 0.5);
        }
        phi_gz(g,z) = varD_gz(g,z) / pow( maxsdD_gz(g,z), 2);
        varB_gbar_z(z) += pow(a_gl(g,0),2) * varD_gz(g,z) * propB_gz(g,z);
        phi_gbar_z(z) += phi_gz(g,z) * propB_gz(g,z);
      }
      // Spatial-buffering index
      for( c=0; c<n_c; c++ ){
        for( g=0; g<n_g; g++ ){
          maxsdB_cz(c,z) += a_gl(g,0) * pow(varD_gcz(g,c,z), 0.5);
        }
        phi_cz(c,z) = varB_cz(c,z) / pow( maxsdB_cz(c,z), 2);
        varB_cbar_z(z) += varB_cz(c,z) * propB_cz(c,z);
        phi_cbar_z(z) += phi_cz(c,z) * propB_cz(c,z);
      }
      // Spatial and species-buffering index
      for( c=0; c<n_c; c++ ){
        for( g=0; g<n_g; g++ ){
          maxsdB_z(z) += a_gl(g,0) * pow(varD_gcz(g,c,z), 0.5);
        }
      }
      phi_z(z) = varB_z(z) / pow( maxsdB_z(z), 2);
    }
    ln_varB_gbar_z = log( varB_gbar_z );
    ln_varB_cbar_z = log( varB_cbar_z );
    ln_varB_z = log( varB_z );
    REPORT( B_y );
    REPORT( D_gy );
    REPORT( B_cy );
    REPORT( phi_gz );
    REPORT( phi_gbar_z );
    REPORT( phi_cz );
    REPORT( phi_cbar_z );
    REPORT( phi_z );
    REPORT( propB_gz );
    REPORT( propB_cz );
    REPORT( varD_gcz );
    REPORT( varD_gz );
    REPORT( varB_cz );
    REPORT( varB_z );
    REPORT( varB_gbar_z );
    REPORT( varB_cbar_z );
    REPORT( maxsdB_z );
    REPORT( maxsdD_gz );
    REPORT( maxsdB_cz );
    ADREPORT( varB_gbar_z );
    ADREPORT( varB_cbar_z );
    ADREPORT( varB_z );
    ADREPORT( B_y );
    ADREPORT( ln_varB_gbar_z );
    ADREPORT( ln_varB_cbar_z );
    ADREPORT( ln_varB_z );
    ADREPORT( phi_gbar_z );
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
    if( FieldConfig(1,0)>0 ) CovHat += loadings_matrix(L_epsilon1_z, n_c, FieldConfig(1,0)) * loadings_matrix(L_epsilon1_z, n_c, FieldConfig(1,0)).transpose();
    if( FieldConfig(1,1)>0 ) CovHat += loadings_matrix(L_epsilon2_z, n_c, FieldConfig(1,1)) * loadings_matrix(L_epsilon2_z, n_c, FieldConfig(1,1)).transpose();
    // Coherence ranges from 0 (all factors are equal) to 1 (first factor explains all variance)
    SelfAdjointEigenSolver<Matrix<Type,Dynamic,Dynamic> > es(CovHat);
    vector<Type> eigenvalues_c = es.eigenvalues();       // Ranked from lowest to highest for some reason
    Type psi = 0;
    for(c=0; c<n_c; c++) psi += eigenvalues_c(n_c-c-1) * (n_c - c);
    psi = 2 * ((psi / eigenvalues_c.sum() / n_c) - 0.5);
    // Total variance
    vector<Type> diag_CovHat( n_c );
    vector<Type> log_diag_CovHat( n_c );
    for(c=0; c<n_c; c++) diag_CovHat(c) = CovHat(c,c);
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
      for(c=0; c<n_c; c++){
        sumtemp += Index_cyl(c,y,l);
      }
      for(c=0; c<n_c; c++){
        PropIndex_cyl(c,y,l) = Index_cyl(c,y,l) / sumtemp;
      }
    }}
    ln_PropIndex_cyl = log( PropIndex_cyl );
    REPORT( PropIndex_cyl );
    REPORT( ln_PropIndex_cyl );
    ADREPORT( PropIndex_cyl );
    ADREPORT( ln_PropIndex_cyl );
  }

  ////////////////////////
  // Diagnostic outputs
  ////////////////////////

  REPORT( Q1 );
  REPORT( Q2 );
  REPORT( B_ff );
  REPORT( P1_iz );
  REPORT( P2_iz );
  REPORT( R1_i );
  REPORT( R2_i );
  REPORT( P1_gcy );
  REPORT( P2_gcy );
  REPORT( var_i );
  REPORT( LogProb1_i );
  REPORT( LogProb2_i );
  REPORT( a_average );
  REPORT( eta1_gct );
  REPORT( eta2_gct );
  REPORT( eta1_vc );
  REPORT( eta2_vc );
  REPORT( eta1_vf );
  REPORT( eta2_vf );
  REPORT( zeta1_i );
  REPORT( zeta2_i );
  REPORT( iota_ct );

  REPORT( SigmaM );
  REPORT( sigmaXi1_cp );
  REPORT( sigmaXi2_cp );
  REPORT( Xi1_scp );
  REPORT( Xi2_scp );
  REPORT( Xi1_gcp );
  REPORT( Xi2_gcp );
  REPORT( Beta_rho1_f );
  REPORT( Beta_mean1_c );
  REPORT( Epsilon_rho1_f );
  REPORT( Beta_rho2_f );
  REPORT( Beta_mean2_c );
  REPORT( Epsilon_rho2_f );

  REPORT( Index_cyl );
  REPORT( D_gcy );
  REPORT( R1_gcy );
  REPORT( R2_gcy );
  REPORT( Index_gcyl );
  REPORT( Omega1_sc );
  REPORT( Omega2_sc );
  REPORT( Omegainput1_sf );
  REPORT( Omegainput2_sf );
  REPORT( Omega1_gc );
  REPORT( Omega2_gc );
  REPORT( Epsilon1_sct );
  REPORT( Epsilon2_sct );
  REPORT( Epsiloninput1_sft );
  REPORT( Epsiloninput2_sft );
  REPORT( Epsilon1_gct );
  REPORT( Epsilon2_gct );
  REPORT( H );
  REPORT( Range_raw1 );
  REPORT( Range_raw2 );
  REPORT( beta1_mean_tf );
  REPORT( beta2_mean_tf );
  REPORT( beta1_tc );
  REPORT( beta2_tc );
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( Options );
  REPORT( Options_vec );
  REPORT( yearbounds_zz );
  REPORT( Expansion_cz );

  ADREPORT( Range_raw1 );
  ADREPORT( Range_raw2 );
  ADREPORT( Index_cyl );
  ADREPORT( ln_Index_cyl );

  SIMULATE{
    REPORT( b_i );
  }

  // Additional miscellaneous outputs
  if( Options(0)==1 ){
    ADREPORT( Index_gcyl );
  }
  if( Options(1)==1 ){
    ADREPORT( log(Index_gcyl) );
  }
  if( Options(3)==1 ){
    vector<Type> D_i( n_i );
    D_i = R1_i * R2_i;
    ADREPORT( D_i );
  }
  // Calculate value of vactors at extrapolation-grid cells (e.g., for use when visualizing estimated or rotated factor estimates)
  if( Options(12)==1 ){
    array<Type> Omegainput1_gf( n_g, Omegainput1_sf.cols() );
    array<Type> Epsiloninput1_gft( n_g, Epsiloninput1_sft.col(0).cols(), n_t );
    array<Type> Omegainput2_gf( n_g, Omegainput2_sf.cols() );
    array<Type> Epsiloninput2_gft( n_g, Epsiloninput2_sft.col(0).cols(), n_t );
    Omegainput1_gf = project_knots( n_g, Omegainput1_sf.cols(), int(1), int(0), Omegainput1_sf, Ags_ij, Ags_x );
    Epsiloninput1_gft = project_knots( n_g, Epsiloninput1_sft.col(0).cols(), n_t, int(1), Epsiloninput1_sft, Ags_ij, Ags_x );
    Omegainput2_gf = project_knots( n_g, Omegainput2_sf.cols(), int(1), int(0), Omegainput2_sf, Ags_ij, Ags_x );
    Epsiloninput2_gft = project_knots( n_g, Epsiloninput2_sft.col(0).cols(), n_t, int(1), Epsiloninput2_sft, Ags_ij, Ags_x );
    REPORT( Omegainput1_gf );
    REPORT( Epsiloninput1_gft );
    REPORT( Omegainput2_gf );
    REPORT( Epsiloninput2_gft );
  }


  return jnll;
}
