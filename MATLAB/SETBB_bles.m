function l_hat_setbb = SETBB_bles(X,Y,statistic,estimator_type,quantile)
% This function implements block length selection for the smooth extended tapered block bootstrap (SETBB) method, 
% provided by Gregory et al. (2018).
% Input:
%        (1) X- Input in the data pair
%        (2) Y- Output in the data pair
%        (3) statistic- Statistic calculated from the data ( for modified SETBB is a quantile regression estimator, 
%                       and for the regular SETBB, sample quantiles or other smooth estimators)
%        (4) estimator_type- Bootstrapped functional (for extended SETBB is ‘distribution’, 
%                                                     and for regular SETBB are ‘variance’ and ‘distribution’)
%        (5) quantile- The quantile of interest for the quantile regression estimator
% Output:
%        (1) l_hat_setbb- Optimal block length selected by the SETBB method


end
