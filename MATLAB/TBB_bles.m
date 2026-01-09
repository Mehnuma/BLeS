function l_hat = TBB_bles(ts,method)
% This function selects the optimal block length for the tapered block bootstrap (TBB), proposed by Paparoditis and Politis (2001,2002), 
% and the extended tapered block bootstrap (ETBBB) by Shao (2010).
% Input:
%        (1) ts- The time series/correlated series of observations 
%        (2) statistic- A quantity directly calculated from the data
%        (3) estimator_type- Bootstrapped functional (Options: ‘variance’ and ‘distribution’)
%        (4) method- TBB or ETBB (Options: ‘regular’ and ‘extended’)
% Output:
%        (1) l_hat- Optimal block length selected for the TBB or ETBB method

[N,~] = size(ts);
K_N = max(5, ceil(sqrt(log10(N))));
m_max = ceil(sqrt(N))+K_N;
c = 2;
m_hat = get_mhat(ts,m_max,c,K_N);
M = 2*m_hat;

% Calculate Gamma_hat and Delta_hat
G_hat = nan(M+1,1);
D_hat = nan(M+1,1);
for k = 0:M
G_hat(k+1) = lambda(k/M).*k^2.*R_hat(ts,k);
D_hat(k+1) = lambda(k/M).*R_hat(ts,k);
end
Gamma_hat = sum(G_hat)*(-5.45);
Delta_hat = 2*((sum(D_hat))^2)*1.1;

if strcmp(method,'TBB') || strcmp(method,'ETBB')
    l_hat = round((4*Gamma_hat^2/Delta_hat)^(1/5)*N^(1/5));
else
    errordlg('Please use TBB or ETBB for the method argument')
end
end

%% Helper Functions
%Calculate m_hat
function m_hat = get_mhat(ts,m_max,c,K_N)
[n,~] = size(ts);
rho_crit = c*sqrt(log10(n)/n);

sample_acf = acf(ts,m_max);
acf_lags = getLags(sample_acf,K_N);
acf_lags = acf_lags(K_N+1:end,:);
indices = abs(acf_lags)<rho_crit;
indices = sum(indices,2);

sig_check = find(indices==K_N,1);
if isempty(sig_check)
    m_hat = find(abs(sample_acf)>rho_crit,1,'last');
else
    m_hat = sig_check;
end
end

%Calculate lambda
function lamb_val = lambda(t)
    lamb_val = (0<=abs(t) && abs(t)<=0.5) +2.*(1-abs(t)).*(0.5<abs(t) && abs(t)<=1);
end

%Calculate R_hat
function R_val= R_hat(ts,k)
    x_bar = mean(ts);
    [N,~] = size(ts);
    term1 = ts(1:(N-abs(k)))-x_bar;
    term2 = ts((1+abs(k)):N)-x_bar;
    R_val = 1/N *sum(term1.*term2);
end

function ta = acf(y,p)
% ACF - Compute Autocorrelations Through p Lags
% >> myacf = acf(y,p) 
%
% Author: Calvin Price (2026). 
% Source: Autocorrelation Function (ACF) (https://www.mathworks.com/matlabcentral/fileexchange/30540-autocorrelation-function-acf), MATLAB Central File Exchange.
% Inputs:
% y - series to compute acf for, nx1 column vector
% p - total number of lags, 1x1 integer
%
% Output:
% myacf - px1 vector containing autocorrelations
%        (First lag computed is lag 1. Lag 0 not computed)
%
%
% A bar graph of the autocorrelations is also produced, with
% rejection region bands for testing individual autocorrelations = 0.
%
% Note that lag 0 autocorelation is not computed, 
% and is not shown on this graph.
%
% Example:
% >> acf(randn(100,1), 10)
%
% --------------------------
% USER INPUT CHECKS
% --------------------------
[n1, n2] = size(y) ;
if n2 ~=1
    error('Input series y must be an nx1 column vector')
end
[a1, a2] = size(p) ;
if ~((a1==1 & a2==1) & (p<n1))
    error('Input number of lags p must be a 1x1 scalar, and must be less than length of series y')
end
% -------------
% BEGIN CODE
% -------------
ta = zeros(p,1) ;
global N 
N = max(size(y)) ;
global ybar 
ybar = mean(y); 
% Collect ACFs at each lag i
for i = 1:p
   ta(i) = acf_k(y,i) ; 
end
end
% ---------------
% SUB FUNCTION
% ---------------
function ta2 = acf_k(y,k)
% ACF_K - Autocorrelation at Lag k
% acf(y,k)
%
% Inputs:
% y - series to compute acf for
% k - which lag to compute acf
% 
global ybar
global N
cross_sum = zeros(N-k,1) ;
% Numerator, unscaled covariance
for i = (k+1):N
    cross_sum(i) = (y(i)-ybar)*(y(i-k)-ybar) ;
end
% Denominator, unscaled variance
yvar = (y-ybar)'*(y-ybar) ;
ta2 = sum(cross_sum) / yvar;
end

function lags = getLags(ts,d)
[row_ts, col_ts] = size(ts);
if col_ts>1
    errordlg('Please use a 1-D series')
end
lags = zeros(row_ts, col_ts*d);

for i = 1:d
    lags(i+1:row_ts,i) = ts(1:row_ts-i,1);
end

end
