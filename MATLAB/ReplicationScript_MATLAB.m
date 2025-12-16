% This script contains all the code snippets included in the manuscript titled "BLeS: A MATLAB and Octave Toolbox for Block Length Selection in Block Bootstrap."
% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu


%% Corrected Politis-White (cPW) Method (From Politis and White (2004))
nsims_cpw = 1000;                        % No. of simulations
N1 = 200; N2 = 800;                      % Sample sizes
data_N1R1 = readmatrix('cpw_N1R1.csv');
data_N1R2 = readmatrix('cpw_N1R2.csv');
data_N1R3 = readmatrix('cpw_N1R3.csv');

data_N2R1 = readmatrix('cpw_N2R1.csv');
data_N2R2 = readmatrix('cpw_N2R2.csv');
data_N2R3 = readmatrix('cpw_N2R3.csv');

opt_N1R1_CB = nan(nsims_cpw,1);
opt_N1R2_CB = nan(nsims_cpw,1);
opt_N1R3_CB = nan(nsims_cpw,1);
opt_N2R1_CB = nan(nsims_cpw,1);
opt_N2R2_CB = nan(nsims_cpw,1);
opt_N2R3_CB = nan(nsims_cpw,1);

opt_N1R1_SB = nan(nsims_cpw,1);
opt_N1R2_SB = nan(nsims_cpw,1);
opt_N1R3_SB = nan(nsims_cpw,1);
opt_N2R1_SB = nan(nsims_cpw,1);
opt_N2R2_SB = nan(nsims_cpw,1);
opt_N2R3_SB = nan(nsims_cpw,1);
cpw_index1 =1; cpw_index2 =1;

for sim =1:nsims_cpw
Xt1 = data_N1R1((cpw_index1:(cpw_index1+N1-1))');
Xt2 = data_N1R2((cpw_index1:(cpw_index1+N1-1))');
Xt3 = data_N1R3((cpw_index1:(cpw_index1+N1-1))');

opt_N1R1_CB(sim) = cPW_method(Xt1, 'circularBB');
opt_N1R2_CB(sim) = cPW_method(Xt2, 'circularBB');
opt_N1R3_CB(sim) = cPW_method(Xt3, 'circularBB');

opt_N1R1_SB(sim) = cPW_method(Xt1, 'stationaryBB');
opt_N1R2_SB(sim) = cPW_method(Xt2, 'stationaryBB');
opt_N1R3_SB(sim) = cPW_method(Xt3, 'stationaryBB');
cpw_index1 = cpw_index1+N1;
end

for sim =1:nsims_cpw
Xt1 = data_N2R1((cpw_index2:(cpw_index2+N2-1))');
Xt2 = data_N2R2((cpw_index2:(cpw_index2+N2-1))');
Xt3 = data_N2R3((cpw_index2:(cpw_index2+N2-1))');

opt_N2R1_CB(sim) = cPW_method(Xt1, 'circularBB');
opt_N2R2_CB(sim) = cPW_method(Xt2, 'circularBB');
opt_N2R3_CB(sim) = cPW_method(Xt3, 'circularBB');

opt_N2R1_SB(sim) = cPW_method(Xt1, 'stationaryBB');
opt_N2R2_SB(sim) = cPW_method(Xt2, 'stationaryBB');
opt_N2R3_SB(sim) = cPW_method(Xt3, 'stationaryBB');
cpw_index2 = cpw_index2+N2;
end

%% Calculation of the Table Results
% Theoretical Optimal Block Lengths
bopt_N1R1_CB = 19;
bopt_N1R2_CB = 2;
bopt_N1R3_CB = 6;
bopt_N1R1_SB = 12.0043; 
bopt_N1R2_SB = 1.3106; 
bopt_N1R3_SB = 2.7991; 

bopt_N2R1_CB = 29;
bopt_N2R2_CB = 4;
bopt_N2R3_CB = 9;
bopt_N2R1_SB = 19.0557; 
bopt_N2R2_SB = 2.0805; 
bopt_N2R3_SB = 4.4432; 

cb_ratio_N1R1 = opt_N1R1_CB/bopt_N1R1_CB;
cb_ratio_N1R2 = opt_N1R2_CB/bopt_N1R2_CB;
cb_ratio_N1R3 = opt_N1R3_CB/bopt_N1R3_CB;
sb_ratio_N1R1 = opt_N1R1_SB/bopt_N1R1_SB;
sb_ratio_N1R2 = opt_N1R2_SB/bopt_N1R2_SB;
sb_ratio_N1R3 = opt_N1R3_SB/bopt_N1R3_SB;

cb_ratio_N2R1 = opt_N2R1_CB/bopt_N2R1_CB;
cb_ratio_N2R2 = opt_N2R2_CB/bopt_N2R2_CB;
cb_ratio_N2R3 = opt_N2R3_CB/bopt_N2R3_CB;
sb_ratio_N2R1 = opt_N2R1_SB/bopt_N2R1_SB;
sb_ratio_N2R2 = opt_N2R2_SB/bopt_N2R2_SB;
sb_ratio_N2R3 = opt_N2R3_SB/bopt_N2R3_SB;

cb_mean_N1R1 = mean(cb_ratio_N1R1); cb_std_N1R1 = std(cb_ratio_N1R1);
cb_mean_N1R2 = mean(cb_ratio_N1R2); cb_std_N1R2 = std(cb_ratio_N1R2);
cb_mean_N1R3 = mean(cb_ratio_N1R3); cb_std_N1R3 = std(cb_ratio_N1R3);

sb_mean_N1R1 = mean(sb_ratio_N1R1); sb_std_N1R1 = std(sb_ratio_N1R1);
sb_mean_N1R2 = mean(sb_ratio_N1R2); sb_std_N1R2 = std(sb_ratio_N1R2);
sb_mean_N1R3 = mean(sb_ratio_N1R3); sb_std_N1R3 = std(sb_ratio_N1R3);

cb_mean_N2R1 = mean(cb_ratio_N2R1); cb_std_N2R1 = std(cb_ratio_N2R1);
cb_mean_N2R2 = mean(cb_ratio_N2R2); cb_std_N2R2 = std(cb_ratio_N2R2);
cb_mean_N2R3 = mean(cb_ratio_N2R3); cb_std_N2R3 = std(cb_ratio_N2R3);

sb_mean_N2R1 = mean(sb_ratio_N2R1); sb_std_N2R1 = std(sb_ratio_N2R1);
sb_mean_N2R2 = mean(sb_ratio_N2R2); sb_std_N2R2 = std(sb_ratio_N2R2);
sb_mean_N2R3 = mean(sb_ratio_N2R3); sb_std_N2R3 = std(sb_ratio_N2R3);

disp('***************************************** cPW Results *****************************************')

disp('===================================== Column 3 Results =====================================')
fprintf("Mean for optimal block length ratio,CB (rho = 0.7, N = 200): %f\n", cb_mean_N1R1)
fprintf("Mean for optimal block length ratio,CB (rho = 0.7, N = 800): %f\n", cb_mean_N2R1)
fprintf("Mean for optimal block length ratio,CB (rho = 0.1, N = 200): %f\n", cb_mean_N1R2)
fprintf("Mean for optimal block length ratio,CB (rho = 0.1, N = 800): %f\n", cb_mean_N2R2)
fprintf("Mean for optimal block length ratio,CB (rho = -0.4, N = 200): %f\n", cb_mean_N1R3)
fprintf("Mean for optimal block length ratio,CB (rho = -0.4, N = 800): %f\n", cb_mean_N2R3)
fprintf("\n")

disp('===================================== Column 4 Results =====================================')
fprintf("S.D. for optimal block length ratio,CB (rho = 0.7, N = 200): %f\n", cb_std_N1R1)
fprintf("S.D. for optimal block length ratio,CB (rho = 0.7, N = 800): %f\n", cb_std_N2R1)
fprintf("S.D. for optimal block length ratio,CB (rho = 0.1, N = 200): %f\n", cb_std_N1R2)
fprintf("S.D. for optimal block length ratio,CB (rho = 0.1, N = 800): %f\n", cb_std_N2R2)
fprintf("S.D. for optimal block length ratio,CB (rho = -0.4, N = 200): %f\n", cb_std_N1R3)
fprintf("S.D. for optimal block length ratio,CB (rho = -0.4, N = 800): %f\n", cb_std_N2R3)
fprintf("\n")

disp('===================================== Column 6 Results =====================================')
fprintf("Mean for optimal block length ratio,SB (rho = 0.7, N = 200): %f\n", sb_mean_N1R1)
fprintf("Mean for optimal block length ratio,SB (rho = 0.7, N = 800): %f\n", sb_mean_N2R1)
fprintf("Mean for optimal block length ratio,SB (rho = 0.1, N = 200): %f\n", sb_mean_N1R2)
fprintf("Mean for optimal block length ratio,SB (rho = 0.1, N = 800): %f\n", sb_mean_N2R2)
fprintf("Mean for optimal block length ratio,SB (rho = -0.4, N = 200): %f\n", sb_mean_N1R3)
fprintf("Mean for optimal block length ratio,SB (rho = -0.4, N = 800): %f\n", sb_mean_N2R3)
fprintf("\n")

disp('===================================== Column 7 Results =====================================')
fprintf("S.D. for optimal block length ratio,SB (rho = 0.7, N = 200): %f\n", sb_std_N1R1)
fprintf("S.D. for optimal block length ratio,SB (rho = 0.7, N = 800): %f\n", sb_std_N2R1)
fprintf("S.D. for optimal block length ratio,SB (rho = 0.1, N = 200): %f\n", sb_std_N1R2)
fprintf("S.D. for optimal block length ratio,SB (rho = 0.1, N = 800): %f\n", sb_std_N2R2)
fprintf("S.D. for optimal block length ratio,SB (rho = -0.4, N = 200): %f\n", sb_std_N1R3)
fprintf("S.D. for optimal block length ratio,SB (rho = -0.4, N = 800): %f\n", sb_std_N2R3)

disp('*******************************************************************************************')

%% Plot for cPW Simulations
h1 = histogram(opt_N1R1_CB);
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
xlabel('$\hat{b}_{opt, CB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N1R2_CB);
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
xlabel('$\hat{b}_{opt, CB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N1R3_CB);
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
xlabel('$\hat{b}_{opt, CB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N2R1_CB);
h1.FaceColor = 'r';
h1.FaceAlpha = 0.8;
xlabel('$\hat{b}_{opt, CB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N2R2_CB);
h1.FaceColor = 'r';
h1.FaceAlpha = 0.8;
xlabel('$\hat{b}_{opt, CB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N2R3_CB);
h1.FaceColor = 'r';
h1.FaceAlpha = 0.8;
xlabel('$\hat{b}_{opt, CB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N1R1_SB);
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
xlabel('$\hat{b}_{opt, SB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N1R2_SB);
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
xlabel('$\hat{b}_{opt, SB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N1R3_SB);
h1.FaceColor = 'b';
h1.FaceAlpha = 0.5;
xlabel('$\hat{b}_{opt, SB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N2R1_SB);
h1.FaceColor = 'r';
h1.FaceAlpha = 0.8;
xlabel('$\hat{b}_{opt, SB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N2R2_SB);
h1.FaceColor = 'r';
h1.FaceAlpha = 0.8;
xlabel('$\hat{b}_{opt, SB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

h1 = histogram(opt_N2R3_SB);
h1.FaceColor = 'r';
h1.FaceAlpha = 0.8;
xlabel('$\hat{b}_{opt, SB}$','interpreter','latex')
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')
