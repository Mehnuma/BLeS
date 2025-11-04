function l_hat = BK_method(ts,statistic)
% This function implements the Bühlmann-Künsch (BK) method of block length selection, provided by Bühlmann & Künsch (1999).
% **Input:**
%        (1) ts- The time series/correlated series of observations
%        (2) statistic- Statistic calculated from the data (mean/median)
% **Output:**
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

