function obj = MBB_BK(ts,bl,statistic,B)   
% Input:
% ts = original time-series;
% bl = block length (given);
% B = Number of bootstrap replications;

n_ts = length(ts); % sample size
stat = nan(B,1);   % the statistic to be calculated (from B replications)
rng('default')
for j = 1:B
    ts_star = []; 
    while (length(ts_star)<= n_ts)
        I = randi([1,n_ts-bl-1]);
        ts_star = [ts_star; ts(I:(I+bl-1))];
    end
    if(length(ts_star)> n_ts)
        ts_star = ts_star(1:n_ts, 1);
    end
    if strcmp(statistic, 'mean')
        stat(j) = mean(ts_star);
    elseif strcmp(statistic, 'median')
        stat(j) = median(ts_star);
    elseif strcmp(statistic, 'variance')
        stat(j) = (std(ts_star))^2;
    else
        errordlg('Please supply the any of the supported statistics.')
    end
end
obj = std(stat);
end