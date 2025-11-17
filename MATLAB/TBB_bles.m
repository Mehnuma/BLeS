l_hat = TBB_bles(ts,statistic,estimator_type,method)
% This function selects the optimal block length for the tapered block bootstrap (TBB), proposed by Paparoditis and Politis (2001,2002), 
% and the extended tapered block bootstrap (ETBBB) by Shao (2010).
% Input:
%        (1) ts- The time series/correlated series of observations 
%        (2) statistic- A quantity directly calculated from the data
%        (3) estimator_type- Bootstrapped functional (Options: ‘variance’ and ‘distribution’)
%        (4) method- TBB or ETBB (Options: ‘regular’ and ‘extended’)
% Output:
%        (1) l_hat- Optimal block length selected for the TBB or ETBB method
