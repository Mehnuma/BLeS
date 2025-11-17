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

%% Set Parameter Values
if strcmp(estimator_type, 'bias') || strcmp(estimator_type, 'variance')
    r = 1; c2 = 1;
elseif strcmp(estimator_type, 'one-sided') || strcmp(estimator_type, 'two-sided') ||  strcmp(estimator_type, 'quantile')
    r = 2; c2 = 0.1;
else
    errordlg('Please use the supported estimator types')
end
[n,~] = size(ts);
l = round(n^(1/(r+4)));
m = round(c2*n^(1/3)*l^(2/3));
if isempty(B)
    B = 500;
end
%% Optimal Block Length Calculation
boot_l  = MBB_nppi(ts,statistic,estimator_type,l,  B,theta_hat,threshold,cdf);    % Regular MBB with bl = l
boot_2l = MBB_nppi(ts,statistic,estimator_type,2*l,B,theta_hat,threshold,cdf);    % Regular MBB with bl = 2*l
var_hat = JAB(ts, boot_l,statistic,estimator_type,l,m,B,theta_hat,threshold,cdf); %JAB_var
C1_hat = n*(l^(-r))*var_hat;
C2_hat = 2*l*(boot_l -boot_2l);
opt_l = ((2*(C2_hat^2)/(r*C1_hat))^(1/(r+2)))*n^(1/(r+2));
if opt_l <0.5
    opt_l = 1;
else
    opt_l = round(opt_l);
end
end

%% JAB Variance Estimator
function jab_var = JAB(ts,boot_l,statistic,estimator_type,bl,m,B,theta_hat,threshold,cdf)
[n,~] = size(ts);
N = n-bl+1;    % No. of all possible blocks
M = N-m+1;     % No. of remaining blocks after deleting m blocks
jab_est = nan(M,1);
all_blocks = cell(N,1);
for i=1:N
    all_blocks{i} = ts(i:i+bl-1);
end
parfor i=1:M
    Ii = (setdiff((1:N),(i:i+m-1)))';
    remaining_blocks = all_blocks(Ii);
    jab_est(i) = iJackVal(remaining_blocks,n,bl,statistic,estimator_type,B,theta_hat,threshold,cdf);
end
phi_tilde = ((N*boot_l)-((N-m)*jab_est))/m;
jab_var = (m/(N-m))*mean((phi_tilde-boot_l).^2);
end
