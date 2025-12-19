function opt_l = HHJ_method(ts,l_init,m,statistic,estimator_type,B,theta_hat,threshold,n_iter)
% This function presents a faster implementation of the Hall-Horowitz-Jing (HHJ) block length selection technique, provided by Hall, Horowitz, and Jing (1995).
% Input:
%        (1) ts- The time series/correlated series of observations 
%        (2) l_init- Initial block length for the moving block bootstrap (MBB) estimator
%        (3) m- Subsample size
%        (4) statistic- Statistic calculated from the data (mean/median)
%        (5) estimator_type- Bootstrapped functional
%        (6) B- Number of bootstrap replications
%        (7) theta_hat-  "True" value used for the bias and distribution estimation
%        (8) threshold- Cut-off point needed for the distribution estimation
%        (9) n_iter- Number of iterations for the HHJ algorithm
% Output:
%        (1) opt_l- Optimal block length selected by the HHJ method

[n,~] = size(ts);
% --- convert estimator_type and statistic strings to numeric flags once ---
switch lower(estimator_type)
    case 'bias',       est_flag = 1;
    case 'variance',   est_flag = 2;
    case 'one-sided',  est_flag = 3;
    case 'two-sided',  est_flag = 4;
    otherwise, error('Unsupported estimator_type');
end
statistic = lower(statistic);
if ~ismember(statistic, {'mean','median','variance'})
    error('Please supply one of: ''mean'', ''median'', ''variance''.');
end

% defaults
if isempty(l_init), l_init = round(n^(1/3)); end
if isempty(m),      m = round(2*n^(1/2)); end
if m>=n, error('sub-sample size(m) must be smaller than the original time-series length(n)'); end
if isempty(B),      B = 100; end
if isempty(n_iter), n_iter = 10; end

% rng('default');
opt_l = [];
l_opt = nan(n_iter,1);

% create sliding windows matrix once (m x (n-m+1)) using indexing
winCount = n - m + 1;
idx = bsxfun(@plus, (1:m)', 0:(winCount-1));  % m x winCount
windows = ts(idx);                              % m x winCount
for iter = 1:n_iter
    % compute ts_stat on full ts with current l_init
    ts_stat = mbb_hhj_fast(ts, l_init, statistic, est_flag, B, theta_hat, threshold);    
    max_lm = round(m/3);
    mse_stat = nan(max_lm,1);
    for lm = 1:max_lm
        % For each sliding window (columns of 'windows'), compute boot result in parallel
        boot_result = nan(winCount,1);
        parfor j = 1:winCount
            col_ts = windows(:,j);         % m x 1 vector
            boot_result(j) = mbb_hhj_fast(col_ts, lm, statistic, est_flag, B, theta_hat, threshold);
        end
        % mean squared error across windows
        mse_stat(lm) = mean((ts_stat - boot_result).^2);
    end
    [~, lm_hat] = min(mse_stat);
    k = estimator_k_from_flag(est_flag); % 3/4/5 depending on estimator_type
    l_opt(iter) = round((n/m)^(1/k) * lm_hat);
    if iter > 1 && l_opt(iter) == l_opt(iter-1)
        opt_l = l_opt(iter);
        break
    end
    l_init = l_opt(iter);
end
% fallback if no convergence
if isempty(opt_l)
    [~, opt_l] = min(mse_stat);
end
end

function k = estimator_k_from_flag(est_flag)
    switch est_flag
        case {1,2}, k = 3;
        case 3,     k = 4;
        case 4,     k = 5;
        otherwise, error('Unknown estimator flag');
    end
end

function obj_mbb = mbb_hhj_fast(ts, bl, statistic, est_flag, B, theta_hat, threshold)
% Faster MBB-based functional used by HHJ_method.
% Inputs: ts (column vector), bl (block length), statistic (string), est_flag numeric.

n_ts = numel(ts);
% number of overlapping blocks
n_block = n_ts - bl + 1;
if n_block < 1
    error('Block length bl must be <= length(ts)');
end

blockIdx = bsxfun(@plus, (1:bl)', 0:(n_block-1));  % bl x n_block
blocks = ts(blockIdx);                              % bl x n_block
L = ceil(n_ts / bl);
stat_vals = nan(B,1);
tao_star = nan(B,1);
parfor boot = 1:B
    sidx = randi(n_block, L, 1);
    ts_star_full = reshape(blocks(:, sidx), [], 1);
    ts_star = ts_star_full(1:n_ts);  % trim to required length
    tao_star(boot) = getTao_fast(ts_star, bl);
    switch statistic
        case 'mean'
            stat_vals(boot) = mean(ts_star);
        case 'median'
            stat_vals(boot) = median(ts_star);
        case 'variance'
            stat_vals(boot) = var(ts_star, 1); % population variance (consistent with (std)^2)
        otherwise
            error('Unsupported statistic');
    end
end
switch est_flag
    case 1  % bias
        obj_mbb = n_ts * (mean(stat_vals) - theta_hat);
    case 2  % variance
        obj_mbb = n_ts * (std(stat_vals,1))^2;
    case 3  % one-sided
        if isempty(threshold), error('Please supply a threshold value'); end
        scaled_stat = sqrt(n_ts) * (stat_vals - theta_hat) ./ tao_star;
        obj_mbb = mean(scaled_stat < threshold);
    case 4  % two-sided
        if isempty(threshold), error('Please supply a threshold value'); end
        scaled_stat = sqrt(n_ts) * (abs(stat_vals - theta_hat) ./ tao_star);
        obj_mbb = mean(scaled_stat < threshold);
    otherwise
        error('Unknown estimator flag');
end
end

function tao_val = getTao_fast(ts, bl)
n = numel(ts);
xbar_star = mean(ts);
b = floor(n / bl);
if b < 1
    sigma_n_star = var(ts,1);
else
    trimmed = ts(1:(b*bl));
    resh = reshape(trimmed, bl, b);              % bl x b
    col_sums = sum(resh, 1);                     % 1 x b
    big_sum = (col_sums - bl * xbar_star).^2;    % 1 x b
    sigma_n_star = sum(big_sum) / (b * bl);
end
tao_val = sqrt(sigma_n_star) + (1 / n);
end
