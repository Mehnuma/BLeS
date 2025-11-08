function b_hat = cPW_method(ts,bb_method)
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
