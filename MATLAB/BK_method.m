function l_hat = BK_method(ts,statistic)
% This function implements the Bühlmann-Künsch (BK) method of block length selection, provided by Bühlmann & Künsch (1999).
% Input:
%        (1) ts- The time series/correlated series of observations
%        (2) statistic- Statistic calculated from the data (mean/median)
% Output:
%        (1) l_hat- Optimal block length selected by the BK method

[n,~] = size(ts);
new_stat = nan(n,1);
ind = (1:n)';

% Jackknife for IF_hat
if strcmp(statistic, 'mean')
    org_stat = mean(ts);
    for i=1:n
        new_ts = ts(ind~=i);
        new_stat(i) = mean(new_ts);
    end
elseif strcmp(statistic, 'median')
    org_stat = median(ts);
    new_stat = sign(ts-org_stat);
elseif strcmp(statistic, 'variance')
    org_stat = (std(ts))^2;
    for i=1:n
        new_ts = ts(ind~=i);
        new_stat(i) = (std(new_ts))^2;
    end
else
    errordlg('Please supply the any of the supported statistics.')
end
IF_hat = n*(org_stat-new_stat);

% Calculate R-hat
R_hat = zeros(n, 1);
t = 1;
for k=0:n-1
    R_hat(k+1)= (1/n)*sum((IF_hat(t:(n-k))).*IF_hat((t+k):n));
end

% Calculation of b's
b = zeros(5,1);
b(1) = 1/n;
% numerator = 2*sum(R_hat.^2)-((R_hat(1)).^2);
numerator = (R_hat(1))^2 +2*sum((R_hat(2:end).^2));
denom1 = zeros(n,1);
for i = 2:5
    for k = 0:(n-1)
        wsc_val = wSC(k.*b(i-1).*n^(4/21));
        denom1(k+1) = (wsc_val).^2.*k^2.*(R_hat(k+1)).^2;
    end
    Denominator = 6*2*sum(denom1)+1e-3;
    b(i) = n^(-1/3)*(numerator/Denominator)^(1/3);
end

% Calculation of l-hat
num = zeros(n,1);
den = zeros(n,1);
for k = 0:(n-1)
    wth_val = wTH(k*b(5)*n^(4/21));
    wsc_val = wSC(k*b(5)*n^(4/21));
    
    num(k+1) = wth_val*R_hat(k+1); 
    den(k+1) = wsc_val*k*R_hat(k+1);
end
big_num = 2*((2*sum(num(2:end))+num(1))^2);
if (sum(den))^2 == 0
    big_den = 1e-3;
else
    big_den = 3*(2*(sum(den))^2);
end
b_hat = n^(-1/3)*(big_num/big_den)^(1/3);

if (1/b_hat) < 0.5
    l_hat = ceil(1/b_hat);
else
    l_hat = round(1/b_hat);
end
end

%% Window Functions
function w1 = wTH(x)
if abs(x)<=1
    w1 = (1+cos(pi*x))./2;
else
    w1 = 0;
end
end

function w2 = wSC(x)
if abs(x)<0.8
    w2 = 1;
elseif ((0.8 <= abs(x)) && (abs(x) <= 1))
    w2 = (1+cos(5*(x-0.8)*pi))./2;
else
    w2 = 0;
end
end

