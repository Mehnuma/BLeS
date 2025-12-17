function l_hat_setbb = SETBB_bles(X,Y,statistic,estimator_type,quantile,variant)
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
%        (6) variant- Tapered block bootstrap ('TBB') or Moving block
%        bootstrap ('MBB')
% Output:
%        (1) l_hat_setbb- Optimal block length selected by the SETBB method

%% Initial Parameters
[n,~] = size(y);
ell1 = round(n^(1/5));
m = floor(n^(1/3)*ell1^(2/3));

%% Block Bootstrap Estimators

%% Calculate nu-hat
nu_hat = n*ell1^(-r)*var_hat;

%% Calculate B-hat and the Optimal Block Length
if strcmp(variant, 'TBB')
    B_hat = (4/3)*ell1^2*(phi_l1 - phi_2l1);
    l_hat_setbb = round((4*(B_hat^2)/nu_hat)^(1/5)*n^(1/5));
elseif strcmp(variant, 'MBB')
    B_hat = 2*ell1*(phi_l1 - phi_2l1);
    l_hat_setbb = round((2*(B_hat^2)/nu_hat)^(1/3)*n^(1/3));
else
    errordlg('Please choose from the TBB or MBB vaiant')
end
end
