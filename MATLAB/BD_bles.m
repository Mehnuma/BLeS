function bopt_BD = BD_bles(ts,p,statistic,method)
% This function implements block length selection for the Bertail_Dudek (BD) method, 
% provided by Bertail and Dudek (2024).
% Input:
%        (1) ts- Periodic time-series
%        (2) p- Period
%        (3) statistic- Parameter to be estimated with options: 'overall mean' and 'seasonal mean'
%        (4) method- 
%            'GSBB' : Generalized seasonal block bootstrap
%            'CGSBB': Circular Generalized seasonal block bootstrap
%            'EMBB': Extended moving block bootstrap
%            'CEMBB': Circular extended moving block bootstrap%        
% Output:
%        (1) bopt_BD- Optimal block length selected by the BD method
% Example:
%	pData = readmatrix('pData.csv');
% 	y_bd = pData(1:200,1);
%	bd_lhat_overall = BD_bles(y_bd,6,'overall mean','EMBB');
%	bd_lhat_seasonal = BD_bles(y_bd,6,'seasonal mean','GSBB');

if strcmp(statistic,'overall mean')
    bopt_BD = BD_overallMean(ts,p,method);
elseif strcmp(statistic,'seasonal mean')
    bopt_BD = BD_seasonalMean(ts,p,method);
else
    errordlg('Please use either overall mean or seasonal mean as the statistic')
end
end
function bopt_overallMean = BD_overallMean(ts,p,method)   
    n = length(ts);
    Kn_min  = floor(n^(1/4));
    Kn_max = floor(n^(1/2));
    Ln = round(4*n^(1/4));
    G = 0;
    sumG = nan(Kn_max-Kn_min+1,1);
    for Kn = Kn_min:Kn_max
        Gacf = getSums(ts,Kn,p,'G');
        G = G+Gacf;
        sumG(Kn-Kn_min+1) = G/(Kn-Kn_min+1);
    end 
    nrowG = floor((size(sumG,1))/2);
    Gn_hat = mean(sumG(nrowG:end));
    fk_hat = nan(p,1);
    for k = 0:(p-1)
        fkCov1 = sum(seasonalACF(ts,1:Ln,p));
        fkCov2 = sum(seasonalACF(ts,0,p));
        fkSum = 2*fkCov1*exp((-2*pi*1i*k*(1:p))/p)'+fkCov2 ;
        fk_hat(k+1) = (1/(2*pi*p))*fkSum;
    end
    Dn_hat = (4*(2*pi*p)^2*sum(abs(fk_hat).^2))/3;

    if strcmp(method, 'EMBB') || strcmp(method, 'CEMBB')
        bopt_overallMean = ceil((2*Gn_hat^2/Dn_hat)^(1/3) * n^(1/3));
    elseif strcmp(method, 'GSBB') || strcmp(method, 'CGSBB')
        bopt_overallMean = ceil((2*Gn_hat^2/Dn_hat)^(1/3) * n^(1/3));
        if mod(bopt_overallMean, p) < p/2
            bopt_overallMean = floor(bopt_overallMean/p)*p + 1;
        else
            bopt_overallMean = (floor(bopt_overallMean/p) + 1)*p - 1;
        end
    else
        errodlg('Please use EMBB, CEMBB, GSBB, or CGSBB as the method')
    end
end
function bopt_seasonalMean = BD_seasonalMean(ts,p,method)
    n = length(ts);
    Kn_min  = floor(n^(1/3));
    Kn_max = floor(n^(1/2));
    Gs = nan(p,Kn_max-Kn_min);
    for Kn = Kn_min:Kn_max
         Gs(:,Kn-Kn_min+1) = getSums(ts,Kn,p,'Gs');
    end
    ncolGs = (size(Gs, 2))/2;
    Gs = Gs(:, ncolGs:end);
    Gs_mean = mean(Gs,2);
    Gn_s_hat = sum(Gs_mean.^2);
    Ds = nan(p,Kn_max-Kn_min);
    for Kn = Kn_min:Kn_max
         Ds(:,Kn-Kn_min+1) = getSums(ts,Kn,p,'Ds');
    end
    ncolDs = (size(Ds, 2))/2;
    Ds = Ds(:, ncolDs:end);
    Ds_mean = mean(Ds,2)*p;
    Dn_s_hat = (4/3)*sum(Ds_mean.^2);

    if strcmp(method, 'EMBB') || strcmp(method, 'CEMBB')
        bopt_seasonalMean = ceil((2*Gn_s_hat/Dn_s_hat)^(1/3)*n^(1/3));
    elseif strcmp(method, 'GSBB') || strcmp(method, 'CGSBB')
        bopt_seasonalMean = ceil((2*Gn_s_hat/Dn_s_hat)^(1/3)*n^(1/3));
        if mod(bopt_seasonalMean, p) < p/2
            bopt_seasonalMean = floor(bopt_seasonalMean/p)*p+1;
        else
            bopt_seasonalMean = (floor(bopt_seasonalMean/p)+1)*p-1;
        end
    else
         errodlg('Please use EMBB, CEMBB, GSBB, or CGSBB as the method')
    end
end
function sum_val = getSums(ts,Kn,p,term)
cov1 = seasonalACF(ts,1:Kn,p);
cov2 = seasonalACF(ts,(Kn-1):-1:0,p);
if strcmp(term, 'G')
    sum_val = sum(cov1'*((1:Kn)/p)')+sum(cov2'*((1:Kn)/p)');
elseif strcmp(term, 'Gs')
    sum_val = cov1'*((1:Kn)/p)'+cov2'*(floor((1:Kn)/p))';
elseif strcmp(term, 'Ds')
    sum_val = cov1'*((1:Kn)/p)'+cov2'*((1:Kn)/p)';
end
end

function acf = seasonalACF(ts,lag,p)
% ChatGPT was used to speed up this code.
    mu_hat = seasonalMean(ts,p); 
    ts = ts(:);
    n = length(ts);
    w  = floor(n/p);
    lagh = length(lag);
    acf  = zeros(lagh, p);
    base = (0:(w-1))*p;
    for i = 1:lagh
        h = lag(i);
        mxi = w*p + h;
        if mxi > n
            ts = [ts; zeros(mxi-n,1)];
        end
        for s = 1:p
            mus  = mu_hat(s);
            mush = mu_hat(mod(s+h-1,p)+1);
            acf(i,s) = mean((ts(base+s)- mus).*(ts(base+s+h)- mush));
        end
    end
    if lagh == 1
        acf = acf(:);
    end
end
function mu = seasonalMean(ts,p)
% ChatGPT was used to speed up this code.
    ts = ts(:);
    n = length(ts);
    if n < p
        p = n;
    end
    mu = zeros(p,1);
    for s = 1:p
        idx = s:p:n;
        mu(s) = mean(ts(idx));
    end
end
