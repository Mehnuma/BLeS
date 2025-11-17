l_hat_setbb = SETBB_bles(X,Y,statistic,estimator_type,quantile)
% This function implements block length selection for the smooth extended tapered block bootstrap (SETBB) method, 
% provided by Gregory et al. (2018).
% Input:
%        (1) ts- The time series/correlated series of observations 
%        (2) l_init- Initial block length for the moving block bootstrap (MBB) estimator
%        (3) m- Subsample size
%        (4) statistic- Statistic calculated from the data (mean/median)
%        (5) estimator_type- Bootstrapped functional
%        (6) B- Number of bootstrap replications
%        (7) theta_hat-  "True" value used for the bias and distribution estimation
%        (8) threshold- Cut-off point needed for the distribution estimation
%        (9) n_iter- Number of iterations for the HHJ algorithm
% Output:
%        (1) opt_l- Optimal block length selected by the HHJ method
