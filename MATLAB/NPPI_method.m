function [opt_l, boot_l] = NPPI_method(ts,statistic,estimator_type,B,theta_hat,threshold,cdf)
% This function implements the nonparametric plugin (NPPI) method by Lahiri, Furukawa, and Lee (2007).
% Input:
%        (1) ts: Dependent series data
%        (2) statistic: Quantity calculated directly from the data
%        (3) estimator_type: Type of bootstrap functional
%        (4) B: Number of Monte-Carlo samples
%        (5) theta_hat: "target" value of the desired statistic
%        (6) threshold: Cut-off point for distribution estimation
%        (7) cdf: Cumulative probability value for bootstrap quantile estimation
% Output:
%        (1) opt_l: Optimal block length from the NPPI method
%        (2) boot_l: Values of the desired bootstrapped functional (e.g., variance of the sample mean)
