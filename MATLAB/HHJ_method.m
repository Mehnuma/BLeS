function opt_l = HHJ_method(ts,l_init,m,statistic,estimator_type,B,theta_hat,threshold,n_iter)
[n,~] = size(ts);
if strcmp(estimator_type, 'bias') || strcmp(estimator_type, 'variance') 
    k = 3;
elseif strcmp(estimator_type, 'one-sided') 
    k = 4;
elseif strcmp(estimator_type, 'two-sided') 
    k = 5;
end
if isempty(l_init)
    l_init = round(n^(1/3)); 
end
if isempty(m) 
    m = round(2*n^(1/2));
end
if m>=n
    errordlg('sub-sample size(m) must be smaller than the original time-series length(n)')
end
if isempty(B)
    B = 100;
end
if isempty(n_iter)
    n_iter = 10;
end
opt_l = [];
l_opt = nan(n_iter,1);
for iter = 1:n_iter
    ts_stat = mbb_hhj(ts,l_init,statistic,estimator_type,B,theta_hat,threshold);
    ts_m = cell(n-m+1,1);
    for i = 1:(n-m+1)
        ts_m{i} = ts(i:(i+m-1));
    end
    boot_result = nan(n-m+1,1);
    mse_stat = nan(m,1);
    for i = 1:round(m/3)
        parfor j = 1:(n-m+1)
            boot_result(j) = mbb_hhj(ts_m{j},i,statistic,estimator_type,B,theta_hat,threshold);
        end
        mse_stat(i) = mean((ts_stat-boot_result).^2);
    end
    [~, lm_hat] = min(mse_stat);
    l_opt(iter) = round((n/m)^(1/k)*lm_hat);
    if iter>1
        if l_opt(iter)==l_opt(iter-1)
            % fprintf('Optimal Block Length: %d\n', l_opt(iter));
            opt_l = l_opt(iter);
            break            
        end
    end
    l_init = l_opt(iter);
end
if isempty(opt_l)
    [~, opt_l] = min(mse_stat);
    % fprintf('Optimal Block Length: %d\n', opt_l);
end
end

function obj_mbb = mbb_hhj(ts,bl,statistic,estimator_type,B,theta_hat,threshold)
% Input:
% ts = Original time-series;
% bl = Block length (given);
% statistic = The statistic to be calculated (from B replications)
% B = Number of bootstrap replications;

[n_ts,~] = size(ts); % sample size
stat = nan(B,1);
tao_star = nan(B,1);
n_block = n_ts-bl+1;
all_blocks = cell(n_block,1);
for i=1:n_block
    all_blocks{i} = ts(i:i+bl-1);
end
rng('default')
parfor j = 1:B
    ts_star = [];
    while (length(ts_star)<= n_ts)
        I = randi([1,n_block]);
        ts_star = [ts_star; cell2mat(all_blocks(I))]; 
    end
    if(length(ts_star)> n_ts)
        ts_star = ts_star(1:n_ts, 1);
    end
    tao_star(j) = getTao(ts_star,bl);
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
if strcmp(estimator_type,'bias')
    obj_mbb = n_ts*(mean(stat)-theta_hat);
elseif strcmp(estimator_type,'variance')
    obj_mbb = n_ts*(std(stat))^2;
elseif strcmp(estimator_type,'one-sided')
    if isempty(threshold)
        disp('Please supply a threshold value')
        return
    else
        scaled_stat = sqrt(n_ts)*(stat-theta_hat)./tao_star;
        obj_mbb = sum(scaled_stat < threshold)/B;
    end
elseif strcmp(estimator_type,'two-sided')
    if isempty(threshold)
        disp('Please supply a threshold value')
        return
    else
        scaled_stat = sqrt(n_ts)*(abs(stat-theta_hat))./tao_star;
        obj_mbb = sum(scaled_stat < threshold)/B;
    end
end
end
function tao_val = getTao(ts,bl)
    [n,~] = size(ts);
    xbar_star = mean(ts);
    b = floor(n/bl);
    big_sum = nan(b,1);
    j=1;
    for i=1:b
        big_sum(i) = (sum(ts(j:(j+bl-1)))-(bl*xbar_star)).^2;
        j=j+bl;
    end
    sigma_n_star = sum(big_sum)/(b*bl);
    tao_val = sqrt(sigma_n_star)+(1/n);
end