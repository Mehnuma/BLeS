function l_hat = BK_faster(ts, statistic, windowType)
% BK_METHOD_FFTPAR: Extremely fast Bühlmann–Künsch block length selection
% using FFT-based autocovariance estimation and optional parallel execution.
%
% Inputs:
%   ts          - vector or matrix (each column is a separate time series)
%   statistic   - 'mean', 'median', or 'variance'
%   windowType  - 'TH' (Tukey-Hanning) or 'SC' (Spectral-Cosine)
%
% Output:
%   l_hat       - Optimal block length(s)
%
% Reference: Bühlmann & Künsch (1999), Block length selection for the bootstrap
% Note: ChatGPT was used for speeding up the code.

if isrow(ts)
    ts = ts(:); % ensure column
end
[n, m] = size(ts);
l_hat = zeros(1, m);

%% Run each time series in parallel
parfor j = 1:m
    x = ts(:, j);
    l_hat(j) = computeBK_single(x, statistic, windowType);
end
end

%% ------------------------------------------------------------------------
function l_hat = computeBK_single(ts, statistic, windowType)
n = numel(ts);

%% 1. Influence function IF_hat
switch lower(statistic)
    case 'mean'
        org_stat = mean(ts);
        new_stat = (n * org_stat - ts) / (n - 1);

    case 'median'
        org_stat = median(ts);
        new_stat = sign(ts - org_stat);

    case 'variance'
        org_stat = var(ts, 1);
        ts2 = ts.^2;
        sum_x = sum(ts);
        sum_x2 = sum(ts2);
        new_stat = ((sum_x2 - ts2) - (sum_x - ts).^2 / (n - 1)) / (n - 2);

    otherwise
        error('Unsupported statistic. Use "mean", "median", or "variance".');
end

IF_hat = n * (org_stat - new_stat);
IF_hat = IF_hat - mean(IF_hat);

%% 2. FFT-based autocovariance estimation
% Zero-pad to next power of 2 for efficiency
NFFT = 2^nextpow2(2*n);
F = fft(IF_hat - mean(IF_hat), NFFT);
acf = real(ifft(abs(F).^2)) / n;
acf = acf(1:n); % keep first n lags

maxLag = floor(n / 4);
autoCov = acf(1:maxLag + 1);

%% 3. Window weights
l_vec = (0:maxLag) / maxLag;
switch upper(windowType)
    case 'TH' % Tukey–Hanning
        w = (abs(l_vec) <= 1) .* (1 + cos(pi * l_vec)) / 2;
    case 'SC' % Spectral–Cosine
        w = zeros(size(l_vec));
        idx1 = abs(l_vec) < 0.8;
        idx2 = abs(l_vec) >= 0.8 & abs(l_vec) <= 1;
        w(idx1) = 1;
        w(idx2) = (1 + cos(5 * (abs(l_vec(idx2)) - 0.8) * pi)) / 2;
    otherwise
        error('Unknown window type. Use ''TH'' or ''SC''.');
end

autoCov_win = autoCov .* w(:);

%% 4. Spectral density estimate at frequency zero
f0_hat = autoCov_win(1) + 2 * sum(autoCov_win(2:end));

%% 5. Optimal block length
var_IF = var(IF_hat, 1);
c_star = (2 * f0_hat / var_IF)^2;
l_hat = c_star^(1/3) * n^(1/3);

% Keep it within reasonable range
l_hat = max(2, min(round(l_hat), n/2));
end
