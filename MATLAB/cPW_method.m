function b_hat = cPW_method(ts,bb_method)
% This function implements the Politis-White method of block size selection (Politis & White, 2004) with correction (Patton, Politis, & White, 2009)
% Input:
%        (1) ts- Observations from a dependent series
%        (2) bb_method- Selection of the block bootstrap variant ('Circular' and 'Stationary')
% Output:
%        (1) b_hat- Optimal block length selected by the cPW method

[N,~] = size(ts);

K_N = max(5, ceil(sqrt(log10(N))));
m_max = ceil(sqrt(N))+K_N;
b_max = ceil(min(3*sqrt(N),N/3));
c = 2;
m_hat = get_mhat(ts,m_max,c,K_N);
M = 2*m_hat;

G_hat = nan(M+1,1);
g_hat0 = nan(M+1,1);
rhat = nan(M+1,1);

for k = 0:M
rhat(k+1) = R_hat(ts,k);
G_hat(k+1) = lambda(k/M).*k.*R_hat(ts,k);
g_hat0(k+1) = lambda(k/M).*R_hat(ts,k);
end

G_hat = G_hat(1)+2*sum(G_hat(2:end));
g_hat0 = g_hat0(1)+2*sum(g_hat0(2:end));

if strcmp(bb_method,'circularBB')    
    D_CB = 4/3*(g_hat0)^2;
    b_hat = (((2*(G_hat)^2)/D_CB)*N)^(1/3);
    if b_hat < 0.5
        b_hat = 1;
    else 
        b_hat = round(b_hat);
    end
    if b_hat > b_max
        b_hat = b_max;
    end
    % fprintf('Optimal block length for circular block bootstrap is: %f\n', b_hat);

elseif strcmp(bb_method,'stationaryBB')
    D_SB = 2*(g_hat0)^2;
    b_hat = (((2*(G_hat)^2)/D_SB)*N)^(1/3);
    if b_hat < 0.5
        b_hat = 1;
    % else 
    %     b_hat = round(b_hat);
    end
    if b_hat > b_max
        b_hat = b_max;
    end
    % fprintf('Optimal block length for stationary block bootstrap is: %f\n', b_hat);
else
    errordlg('Please choose from circular and stationary bootstrap')
end
end

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

%%
function ta = acf(y,p)
% ACF - Compute Autocorrelations Through p Lags
% Source- Calvin Price (2025). Autocorrelation Function (ACF) (https://www.mathworks.com/matlabcentral/fileexchange/30540-autocorrelation-function-acf), 
% MATLAB Central File Exchange.
% >> myacf = acf(y,p) 
%
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
% Note that lag 0 autocorrelation is not computed, 
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
