function l_hat = SETBB_bles(X,Y,quantile,boot_variant,B)
% This function implements block length selection for the smooth extended tapered block bootstrap (SETBB) method,
% provided by Gregory et al. (2018).
% Input:
%        (1) X- Input in the data pair
%        (2) Y- Output in the data pair
%        (3) quantile- The quantile of interest for the quantile regression estimator
%        (4) boot_variant- Smooth Extended Tapered block bootstrap ('SETBB') or Moving block
%        bootstrap ('MBB')
%        (5) B- Number of Monte-Carlo replications
% Output:
%        (1) l_hat- Optimal block length selected for the quantile regression (QR) estimator by the SETBB or MBB method

%% Initial Parameters
[n,~] = size(Y);
d = size(X,2)-1;
ell1 = round(n^(1/5));
m = floor(n^(1/3)*ell1^(2/3));
ell2 = 2*ell1;
if isempty(B); B=1000; end

model = rq(X,Y,quantile);
beta_hat = model.coefficients;
[~,~,bw] = kde(model.residuals,Bandwidth='plug-in');

%% Calculate Optimal Block Length
if strcmp(boot_variant, 'SETBB')
    piTildeSETBB_ell1 = getPiTilde(n,ell1,'SETBB');
    piTildeSETBB_ell2 = getPiTilde(n,ell2,'SETBB');
    betaTildeSETBB_ell1 = fminunc(@(b) calcStilde(X,Y,b,piTildeSETBB_ell1,quantile,bw),beta_hat);
    betaTildeSETBB_ell2 = fminunc(@(b) calcStilde(X,Y,b,piTildeSETBB_ell2,quantile,bw),beta_hat);
    [~,piStarSETBB_ell1,Points1,~] = getPiStar(n,ell1,B);
    [~,piStarSETBB_ell2,~,~] = getPiStar(n,ell2,B);

    X_pert = X+normrnd(0,bw,[n,d]);
    X_pert(:,1) = X(:,1);
    Y_pert = Y+normrnd(0,bw,[n,1]);
    DnSETBB_ell1 = calcDn(X_pert,Y_pert,betaTildeSETBB_ell1,quantile,piStarSETBB_ell1,B,d);
    DnSETBB_ell2 = calcDn(X_pert,Y_pert,betaTildeSETBB_ell2,quantile,piStarSETBB_ell2,B,d);
    phiSETBB_ell1 = sum(diag(cov(DnSETBB_ell1)));
    phiSETBB_ell2 = sum(diag(cov(DnSETBB_ell2)));

    N = n-ell1+1; M=N-m+1; jabPhiMBB_ell1=nan(M,1);
    for k = 1:M
        blockIndex = k:(k+m-1);
        kIndex = ismember(Points1(k,:),blockIndex);
        jabPhiMBB_ell1(k) = sum(diag(cov(DnSETBB_ell1(kIndex))));
    end
    jabPhiMBB_ell1 = jabPhiMBB_ell1(~isnan(jabPhiMBB_ell1));
    jabPhiTildeSETBB = (1/m)*(N*phiSETBB_ell1-(N-m)*jabPhiMBB_ell1);
    var_hat_SETBB = (m/(N-m))*mean((jabPhiTildeSETBB-phiSETBB_ell1).^2);
    nu_hat_SETBB =  (n/ell1)*var_hat_SETBB;
    B_hat_SETBB = (4/3)*ell1^2*(phiSETBB_ell1 - phiSETBB_ell2);
    l_hat = round((4*(B_hat_SETBB^2)/nu_hat_SETBB)^(1/5)*n^(1/5));
    if l_hat<0.5
        l_hat = 1;
    end
elseif strcmp(boot_variant, 'MBB')
    piTildeMBB_ell1 = getPiTilde(n,ell1,'MBB');
    piTildeMBB_ell2 = getPiTilde(n,ell2,'MBB');
    model1 = rq(X,Y,quantile,'Weights',piTildeMBB_ell1);
    betaTildeMBB_ell1 = model1.coefficients;
    model2= rq(X,Y,quantile,'Weights',piTildeMBB_ell2);
    betaTildeMBB_ell2 = model2.coefficients;
    [piStarMBB_ell1,~,Points1,~] = getPiStar(n,ell1,B);
    [piStarMBB_ell2,~,~,~] = getPiStar(n,ell2,B);

    X_pert = X+normrnd(0,bw,[n,d]);
    X_pert(:,1) = X(:,1);
    Y_pert = Y+normrnd(0,bw,[n,1]);
    DnMBB_ell1 = calcDn(X_pert,Y_pert,betaTildeMBB_ell1,quantile,piStarMBB_ell1,B,d);
    DnMBB_ell2 = calcDn(X_pert,Y_pert,betaTildeMBB_ell2,quantile,piStarMBB_ell2,B,d);
    phiMBB_ell1 = sum(diag(cov(DnMBB_ell1)));
    phiMBB_ell2 = sum(diag(cov(DnMBB_ell2)));

    N = n-ell1+1; M=N-m+1; jabPhiMBB_ell1=nan(M,1);
    for k = 1:M
        blockIndex = k:(k+m-1);
        kIndex = ismember(Points1(k,:),blockIndex);
        jabPhiMBB_ell1(k) = sum(diag(cov(DnMBB_ell1(kIndex))));
    end
    jabPhiMBB_ell1 = jabPhiMBB_ell1(~isnan(jabPhiMBB_ell1));
    jabPhiTildeMBB = (1/m)*(N*phiMBB_ell1-(N-m)*jabPhiMBB_ell1);
    var_hat_MBB = (m/(N-m))*mean((jabPhiTildeMBB-phiMBB_ell1).^2);
    nu_hat_MBB =  n/ell1*var_hat_MBB;
    B_hat_MBB = 2*ell1*(phiMBB_ell1 - phiMBB_ell2);
    l_hat = round((2*(B_hat_MBB^2)/nu_hat_MBB)^(1/3)*n^(1/3));
    if l_hat<0.5
        l_hat = 1;
    end
else
    errordlg('Please choose from the SETBB or MBB vaiant')
end
end
function pi_tilde = getPiTilde(n,ell,boot_variant)
if strcmp(boot_variant, 'MBB')
    if ell == 1
        pi_tilde = ones(1,n)/n;
    else
        pi_tilde = zeros(1, n);
        for t = 1:(ell-1)
            pi_tilde(t) = (t/ell)/(n-ell+1);
            pi_tilde(n-t+1)= pi_tilde(t);
        end
        pi_tilde(ell:(n-ell+1)) = 1/(n-ell+1);
    end
elseif strcmp(boot_variant, 'SETBB')
    if ell==1, pi_tilde = ones(n,1)/n; return; end
    u=((1:ell)-0.5)/ell;
    omega_vec = omega(u,0.43);
    omega_vec = omega_vec/sum(omega_vec);
    pi_tilde=zeros(n,1);
    for t=1:ell-1
        pi_tilde(t)=sum(omega_vec)/(n-ell+1);
        pi_tilde(n-t+1)=pi_tilde(t);
    end
    pi_tilde(ell:n-ell+1)=1/(n-ell+1);
else
    errordlg('Please choose from the SETBB or MBB vaiant')
end
end
function S_tilde = calcStilde(X,Y,beta,piTilde,tau,h)
var_eta = h^2*(1 + sum(beta(2:end).^2));
e = Y - X*beta;
Phi = normcdf(e,0,sqrt(var_eta));
phi = normpdf(e,0,sqrt(var_eta));
S_tilde = sum(piTilde.*(e.*(tau-(1-Phi)) + var_eta.*phi));
end
function Dn_star = calcDn(X,Y,betaTilde,tau,piStar,B,d)
Dn_star = zeros(B,d+1);
piStar = piStar';
for i=1:B
    Dn_star(i,:) = sum(piStar(i,1)*X.*signFunc(tau,(Y-X*betaTilde)));
end
end
function val = signFunc(tau,diff)
val = tau*(diff>0)-(1-tau)*(diff<=0);
end
function fit = rq(X, Y, tau, varargin)
% ChatGPT was used to speed up this code.
p = inputParser;
addParameter(p, 'Weights', []);
parse(p, varargin{:});
weights = p.Results.Weights;
[n, pX] = size(X);
tau = unique(sort(tau(:)'));
eps_val = eps^(2/3);
tau(tau == 0) = eps_val;
tau(tau == 1) = 1 - eps_val;
if any(tau < 0 | tau > 1)
    error('Invalid tau: must be in [0,1]');
end
Rho = @(u,t)u.*(t-(u<0));
coef = zeros(pX, numel(tau));
resid = zeros(n, numel(tau));
fitted = zeros(n, numel(tau));
rho = zeros(1, numel(tau));
for i = 1:numel(tau)
    t = tau(i);
    beta = rq_lp(X, Y, t, weights);
    coef(:, i) = beta;
    fitted(:, i) = X * beta;
    resid(:, i) = Y - fitted(:, i);
    rho(i) = sum(Rho(resid(:, i), t));
end
fit.coefficients = coef;
fit.residuals = resid;
fit.fitted_values = fitted;
fit.tau = tau;
fit.rho = rho;
fit.method = 'br (LP)';
fit.weights = weights;
fit.call = 'rq';
end
function beta = rq_lp(X, Y, tau, weights)
% ChatGPT was used to speed up this code.
[n, p] = size(X);
if isempty(weights)
    W = eye(n);
else
    W = diag(weights);
end
f = [zeros(p,1);
    tau * ones(n,1);
    (1 - tau) * ones(n,1)];
Aeq = [X, eye(n), -eye(n)];
beq = Y;
lb = [-inf(p,1); zeros(2*n,1)];
ub = [];
opts = optimoptions('linprog','Display','none');
sol = linprog(f, [], [], Aeq, beq, lb, ub, opts);
beta = sol(1:p);
end
function [weight_MBB, weight_SETBB, blockpts, m_l] = getPiStar(n,l,B)
rng('default')
c=0.43;
maxv = n-l+1;
blocks = floor(n/l);
weights = zeros(l,1);
tmp = ((1:l)'-0.5)/l;

idx1 = tmp <= c; idx2 = tmp > c & tmp < (1-c); idx3 = tmp >= (1-c);
weights(idx1) = tmp(idx1)/c; weights(idx2) = 1; weights(idx3) = (1-tmp(idx3))/c;
weightnorm1 = sum(weights);
weightnorm2 = sum(weights.^2);
scale = 1 / (blocks*weightnorm1);
m_l   = weightnorm1^2/(l*weightnorm2);

blockpts = floor(rand(B,blocks)*maxv);
weight_MBB  = zeros(n,B);
weight_SETBB = zeros(n,B);
for i = 1:B
    blk = blockpts(i,:);
    for t = 1:n
        k = (t-1)-blk;
        valid = (k>=0) & (k<l);
        accumMBB = sum(valid);
        if accumMBB > 0
            accumETBB = sum(weights(k(valid)+1));
        else
            accumETBB = 0;
        end
        weight_MBB(t, i)  = accumMBB/(blocks*l);
        weight_SETBB(t, i) = accumETBB*scale;
    end
end
end
function y = omega(t,c)
if nargin < 2; c = 0.43; end
y = (t<=c).*(t/c)+((t>c) & (t< (1-c)))+(t>=(1-c)).*((1-t)/c);
end
