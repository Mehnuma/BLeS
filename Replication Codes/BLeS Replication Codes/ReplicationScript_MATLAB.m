% This script contains all the code snippets included in the manuscript titled "BLeS: A MATLAB and Octave Toolbox for Block Length Selection in Block Bootstrap."
% Authors: Mehnuma Tabassum, Kris De Brabanter
% Email: mehnuma@iastate.edu


%% Bühlmann-Künsch (BK) Method (From Bühlmann and Künsch (1999))
nsims_bk = 1000;
data_bk128 = readmatrix('bk_bilinear128_new.csv');
data_bk512 = readmatrix('bk_bilinear512_new.csv');
bk_n1 = 128; bk_n2 = 512;
B_bk = 5000;

bk_lhat_mean128 = nan(nsims_bk,1);
ssboot_mean128 = nan(nsims_bk,1);
ssboot_l13_mean128 = nan(nsims_bk,1);
ssboot_lopt_mean128 = nan(nsims_bk,1);

bk_lhat_mean512 = nan(nsims_bk,1);
ssboot_mean512 = nan(nsims_bk,1);
ssboot_l13_mean512 = nan(nsims_bk,1);
ssboot_lopt_mean512 = nan(nsims_bk,1);

bk_lhat_median128 = nan(nsims_bk,1);
ssboot_median128 = nan(nsims_bk,1);
ssboot_l13_median128 = nan(nsims_bk,1);
ssboot_lopt_median128 = nan(nsims_bk,1);

bk_lhat_median512 = nan(nsims_bk,1);
ssboot_median512 = nan(nsims_bk,1);
ssboot_l13_median512 = nan(nsims_bk,1);
ssboot_lopt_median512 = nan(nsims_bk,1);

bk_index1 = 1; bk_index2=1;
for sim =1:nsims_bk
    y128 = data_bk128((bk_index1:(bk_index1+bk_n1-1))');

    bk_lhat_mean128(sim) = BK_method(y128,'mean');    
    ssboot_mean128(sim) = sqrt(bk_n1)*MBB_BK(y128,bk_lhat_mean128(sim),'mean',B_bk);
    ssboot_l13_mean128(sim) = sqrt(bk_n1)*MBB_BK(y128,5,'mean',B_bk);
    ssboot_lopt_mean128(sim) = sqrt(bk_n1)*MBB_BK(y128,6,'mean',B_bk);

    bk_lhat_median128(sim) = BK_method(y128,'median');
    ssboot_median128(sim) = sqrt(bk_n1)*MBB_BK(y128,bk_lhat_median128(sim),'median',B_bk);
    ssboot_l13_median128(sim) = sqrt(bk_n1)*MBB_BK(y128,5,'median',B_bk);
    ssboot_lopt_median128(sim) = sqrt(bk_n1)*MBB_BK(y128,3,'median',B_bk);
    bk_index1 = bk_index1+bk_n1;
end
for sim =1:nsims_bk
    y512 = data_bk512((bk_index2:(bk_index2+bk_n2-1))');

    bk_lhat_mean512(sim) = BK_method(y512,'mean');    
    ssboot_mean512(sim) = sqrt(bk_n2)*MBB_BK(y512,bk_lhat_mean512(sim),'mean',B_bk);
    ssboot_l13_mean512(sim) = sqrt(bk_n2)*MBB_BK(y512,8,'mean',B_bk);
    ssboot_lopt_mean512(sim) = sqrt(bk_n2)*MBB_BK(y512,9,'mean',B_bk);

    bk_lhat_median512(sim) = BK_method(y512,'median');
    ssboot_median512(sim) = sqrt(bk_n2)*MBB_BK(y512,bk_lhat_median512(sim),'median',B_bk);
    ssboot_l13_median512(sim) = sqrt(bk_n2)*MBB_BK(y512,8,'median',B_bk);
    ssboot_lopt_median512(sim) = sqrt(bk_n2)*MBB_BK(y512,5,'median',B_bk);
    bk_index2 = bk_index2+bk_n2;
end

% Calculate "True" Statistics
nsims_bk_true = 1000;
data_bk128_true = readmatrix('bk_bilinear128_truestat_new.csv');
data_bk512_true = readmatrix('bk_bilinear512_truestat_new.csv');

true_stat_mean128 = nan(nsims_bk_true,1); true_stat_median128 = nan(nsims_bk_true,1);
true_stat_mean512 = nan(nsims_bk_true,1); true_stat_median512 = nan(nsims_bk_true,1);

bk_index_true1 = 1; bk_index_true2=1;
for sim =1:nsims_bk_true
    y128_true = data_bk128_true((bk_index_true1:(bk_index_true1+bk_n1-1))');
    true_stat_mean128(sim) = mean(y128_true);
    true_stat_median128(sim) = median(y128_true);
    bk_index_true1 = bk_index_true1+bk_n1;
end

for sim =1:nsims_bk_true
    y512_true = data_bk512_true((bk_index_true2:(bk_index_true2+bk_n2-1))');
    true_stat_mean512(sim) = mean(y512_true);
    true_stat_median512(sim) = median(y512_true);
    bk_index_true2 = bk_index_true2+bk_n2;
end

True_Stat_mean128 = sqrt(bk_n1)*std(true_stat_mean128);
True_Stat_median128 = sqrt(bk_n1)*std(true_stat_median128);
True_Stat_mean512 = sqrt(bk_n2)*std(true_stat_mean512);
True_Stat_median512 = sqrt(bk_n2)*std(true_stat_median512);

% Calculate Bias, Standard Deviation, and RMSE
% Adaptive, Mean, n=128
bias_adaptive_mean128 = mean(ssboot_mean128)-True_Stat_mean128;
sd_ssboot_mean128 = std(ssboot_mean128);
rmse_mean128 = ((ssboot_mean128-True_Stat_mean128).^2)/(True_Stat_mean128)^2;
mean_rmse_mean128 = mean(rmse_mean128);
sd_rmse_mean128 = std(rmse_mean128);
% l^1/3, Mean, n=128
bias_l13_mean128 = mean(ssboot_l13_mean128)-True_Stat_mean128;
sd_ssboot_l13_mean128 = std(ssboot_l13_mean128);
rmse_l13_mean128 = ((ssboot_l13_mean128-True_Stat_mean128).^2)/(True_Stat_mean128)^2;
mean_rmse_l13_mean128 = mean(rmse_l13_mean128);
sd_rmse_l13_mean128 = std(rmse_l13_mean128);
% l_opt, Mean, n=128
bias_lopt_mean128 = mean(ssboot_lopt_mean128)-True_Stat_mean128;
sd_ssboot_lopt_mean128 = std(ssboot_lopt_mean128);
rmse_lopt_mean128 = ((ssboot_lopt_mean128-True_Stat_mean128).^2)/(True_Stat_mean128)^2;
mean_rmse_lopt_mean128 = mean(rmse_lopt_mean128);
sd_rmse_lopt_mean128 = std(rmse_lopt_mean128);

% Adaptive, Median, n=128
bias_adaptive_median128 = mean(ssboot_median128)-True_Stat_median128;
sd_ssboot_median128 = std(ssboot_median128);
rmse_median128 = ((ssboot_median128-True_Stat_median128).^2)/(True_Stat_median128)^2;
mean_rmse_median128 = mean(rmse_median128);
sd_rmse_median128 = std(rmse_median128);
% l^1/3, Median, n=128
bias_l13_median128 = mean(ssboot_l13_median128)-True_Stat_median128;
sd_ssboot_l13_median128 = std(ssboot_l13_median128);
rmse_l13_median128 = ((ssboot_l13_median128-True_Stat_median128).^2)/(True_Stat_median128)^2;
mean_rmse_l13_median128 = mean(rmse_l13_median128);
sd_rmse_l13_median128 = std(rmse_l13_median128);
% l_opt, Median, n=128
bias_lopt_median128 = mean(ssboot_lopt_median128)-True_Stat_median128;
sd_ssboot_lopt_median128 = std(ssboot_lopt_median128);
rmse_lopt_median128 = ((ssboot_lopt_median128-True_Stat_median128).^2)/(True_Stat_median128)^2;
mean_rmse_lopt_median128 = mean(rmse_lopt_median128);
sd_rmse_lopt_median128 = std(rmse_lopt_median128);

% Adaptive, Mean, n=512
bias_adaptive_mean512 = mean(ssboot_mean512)-True_Stat_mean512;
sd_ssboot_mean512 = std(ssboot_mean512);
rmse_mean512 = ((ssboot_mean512-True_Stat_mean512).^2)/(True_Stat_mean512)^2;
mean_rmse_mean512 = mean(rmse_mean512);
sd_rmse_mean512 = std(rmse_mean512);
% l^1/3, Mean, n=512
bias_l13_mean512 = mean(ssboot_l13_mean512)-True_Stat_mean512;
sd_ssboot_l13_mean512 = std(ssboot_l13_mean512);
rmse_l13_mean512 = ((ssboot_l13_mean512-True_Stat_mean512).^2)/(True_Stat_mean512)^2;
mean_rmse_l13_mean512 = mean(rmse_l13_mean512);
sd_rmse_l13_mean512 = std(rmse_l13_mean512);
% l_opt, Mean, n=512
bias_lopt_mean512 = mean(ssboot_lopt_mean512)-True_Stat_mean512;
sd_ssboot_lopt_mean512 = std(ssboot_lopt_mean512);
rmse_lopt_mean512 = ((ssboot_lopt_mean512-True_Stat_mean512).^2)/(True_Stat_mean512)^2;
mean_rmse_lopt_mean512 = mean(rmse_lopt_mean512);
sd_rmse_lopt_mean512 = std(rmse_lopt_mean512);

% Adaptive, Median, n=512
bias_adaptive_median512 = mean(ssboot_median512)-True_Stat_median512;
sd_ssboot_median512 = std(ssboot_median512);
rmse_median512 = ((ssboot_median512-True_Stat_median512).^2)/(True_Stat_median512)^2;
mean_rmse_median512 = mean(rmse_median512);
sd_rmse_median512 = std(rmse_median512);
% l^1/3, Median, n=512
bias_l13_median512 = mean(ssboot_l13_median512)-True_Stat_median512;
sd_ssboot_l13_median512 = std(ssboot_l13_median512);
rmse_l13_median512 = ((ssboot_l13_median512-True_Stat_median512).^2)/(True_Stat_median512)^2;
mean_rmse_l13_median512 = mean(rmse_l13_median512);
sd_rmse_l13_median512 = std(rmse_l13_median512);
% l_opt, Median, n=512
bias_lopt_median512 = mean(ssboot_lopt_median512)-True_Stat_median512;
sd_ssboot_lopt_median512 = std(ssboot_lopt_median512);
rmse_lopt_median512 = ((ssboot_lopt_median512-True_Stat_median512).^2)/(True_Stat_median512)^2;
mean_rmse_lopt_median512 = mean(rmse_lopt_median512);
sd_rmse_lopt_median512 = std(rmse_lopt_median512);

% Table
disp('***************************************** cPW Results *****************************************')

disp('===================================== Column 3 Results =====================================')
fprintf("True Sigma_n for Mean (n=128): %f\n", True_Stat_mean128)
fprintf("True Sigma_n for Mean (n=512): %f\n", True_Stat_mean512)
fprintf("True Sigma_n for Median (n=128): %f\n", True_Stat_median128)
fprintf("True Sigma_n for Median (n=512): %f\n", True_Stat_median512)
fprintf("\n")

disp('===================================== Column 4 Results =====================================')
fprintf("Bias for l_BK (Mean, n=128): %f\n", bias_adaptive_mean128)
fprintf("Bias for l^1/3 (Mean, n=128): %f\n", bias_l13_mean128)
fprintf("Bias for l_opt (Mean, n=128): %f\n", bias_lopt_mean128)
fprintf("Bias for l_BK (Mean, n=512): %f\n", bias_adaptive_mean512)
fprintf("Bias for l^1/3 (Mean, n=512): %f\n", bias_l13_mean512)
fprintf("Bias for l_opt (Mean, n=512): %f\n", bias_lopt_mean512)

fprintf("Bias for l_BK (Median, n=128): %f\n", bias_adaptive_median128)
fprintf("Bias for l^1/3 (Median, n=128): %f\n", bias_l13_median128)
fprintf("Bias for l_opt (Median, n=128): %f\n", bias_lopt_median128)
fprintf("Bias for l_BK (Median, n=512): %f\n", bias_adaptive_median512)
fprintf("Bias for l^1/3 (Median, n=512): %f\n", bias_l13_median512)
fprintf("Bias for l_opt (Median, n=512): %f\n", bias_lopt_median512)
fprintf("\n")

disp('===================================== Column 5 Results =====================================')
fprintf("SD for l_BK (Mean, n=128): %f\n", sd_ssboot_mean128)
fprintf("SD for l^1/3 (Mean, n=128): %f\n", sd_ssboot_l13_mean128)
fprintf("SD for l_opt (Mean, n=128): %f\n", sd_ssboot_lopt_mean128)
fprintf("SD for l_BK (Mean, n=512): %f\n", sd_ssboot_mean512)
fprintf("SD for l^1/3 (Mean, n=512): %f\n", sd_ssboot_l13_mean512)
fprintf("SD for l_opt (Mean, n=512): %f\n", sd_ssboot_lopt_mean512)

fprintf("SD for l_BK (Median, n=128): %f\n", sd_ssboot_median128)
fprintf("SD for l^1/3 (Median, n=128): %f\n", sd_ssboot_l13_median128)
fprintf("SD for l_opt (Median, n=128): %f\n", sd_ssboot_lopt_median128)
fprintf("SD for l_BK (Median, n=512): %f\n", sd_ssboot_median512)
fprintf("SD for l^1/3 (Median, n=512): %f\n", sd_ssboot_l13_median512)
fprintf("SD for l_opt (Median, n=512): %f\n", sd_ssboot_lopt_median512)
fprintf("\n")

disp('===================================== Column 6 Results =====================================')
fprintf("RMSE for l_BK (Mean, n=128): %f(%f)\n", mean_rmse_mean128, sd_rmse_mean128)
fprintf("RMSE for l^1/3 (Mean, n=128): %f(%f)\n", mean_rmse_l13_mean128, sd_rmse_l13_mean128)
fprintf("RMSE for l_opt (Mean, n=128): %f(%f)\n", mean_rmse_lopt_mean128, sd_rmse_lopt_mean128)
fprintf("RMSE for l_BK (Mean, n=512): %f(%f)\n", mean_rmse_mean512, sd_rmse_mean512)
fprintf("RMSE for l^1/3 (Mean, n=512): %f(%f)\n", mean_rmse_l13_mean512, sd_rmse_l13_mean512)
fprintf("RMSE for l_opt (Mean, n=512): %f(%f)\n", mean_rmse_lopt_mean512, sd_rmse_lopt_mean512)

fprintf("RMSE for l_BK (Median, n=128): %f(%f)\n", mean_rmse_median128, sd_rmse_median128)
fprintf("RMSE for l^1/3 (Median, n=128): %f(%f)\n", mean_rmse_l13_median128, sd_rmse_l13_median128)
fprintf("RMSE for l_opt (Median, n=128): %f(%f)\n", mean_rmse_lopt_median128, sd_rmse_lopt_median128)
fprintf("RMSE for l_BK (Median, n=512): %f(%f)\n", mean_rmse_median512, sd_rmse_median512)
fprintf("RMSE for l^1/3 (Median, n=512): %f(%f)\n", mean_rmse_l13_median512, sd_rmse_l13_median512)
fprintf("RMSE for l_opt (Median, n=512): %f(%f)\n", mean_rmse_lopt_median512, sd_rmse_lopt_median512)
fprintf("\n")
disp('*******************************************************************************************')


% Plots
% Mean, n=128
finalStat_mean128 = ssboot_mean128/True_Stat_mean128;
finalStat_l13_mean128  = ssboot_l13_mean128/True_Stat_mean128;
finalStat_lopt_mean128  = ssboot_lopt_mean128/True_Stat_mean128;

% Median, n=128
finalStat_median128 = ssboot_median128/True_Stat_median128;
finalStat_l13_median128  = ssboot_l13_median128/True_Stat_median128;
finalStat_lopt_median128  = ssboot_lopt_median128/True_Stat_median128;

% Mean, n=512
finalStat_mean512 = ssboot_mean512/True_Stat_mean512;
finalStat_l13_mean512  = ssboot_l13_mean512/True_Stat_mean512;
finalStat_lopt_mean512  = ssboot_lopt_mean512/True_Stat_mean512;

% Median, n=512
finalStat_median512 = ssboot_median512/True_Stat_median512;
finalStat_l13_median512  = ssboot_l13_median512/True_Stat_median512;
finalStat_lopt_median512  = ssboot_lopt_median512/True_Stat_median512;


b_colors = ['#eb8410'; '#9a1fae'; '#8beb10'];
b_markerColor = ['#eb8410'; '#9a1fae'; '#8beb10'];
% b_colors = ['#2D2926'; '#FFCD00'; '#C8102E'];
% b_markerColor = ['#2D2926'; '#FFCD00'; '#C8102E'];
x = 1:3;

b1_boxdata = [finalStat_mean128 finalStat_l13_mean128 finalStat_lopt_mean128];
figure();
ax = axes();
hold(ax);
for i=1:numel(x)
    boxchart(x(i)*ones(size(b1_boxdata(:,i))), b1_boxdata(:,i), 'BoxFaceColor', b_colors(i,:), 'BoxFaceAlpha',0.7, 'MarkerColor',b_markerColor(i,:), 'MarkerStyle','.')
end
hold off
box on
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

b2_boxdata = [finalStat_median128 finalStat_l13_median128 finalStat_lopt_median128];
figure();
ax = axes();
hold(ax);
for i=1:numel(x)
    boxchart(x(i)*ones(size(b2_boxdata(:,i))), b2_boxdata(:,i), 'BoxFaceColor', b_colors(i,:), 'BoxFaceAlpha',0.7, 'MarkerColor',b_markerColor(i,:), 'MarkerStyle','.')
end
hold off
box on
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')


b3_boxdata = [finalStat_mean512 finalStat_l13_mean512 finalStat_lopt_mean512];
figure();
ax = axes();
hold(ax);
for i=1:numel(x)
    boxchart(x(i)*ones(size(b3_boxdata(:,i))), b3_boxdata(:,i), 'BoxFaceColor', b_colors(i,:), 'BoxFaceAlpha',0.7, 'MarkerColor',b_markerColor(i,:), 'MarkerStyle','.')
end
hold off
box on
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

b4_boxdata = [finalStat_median512 finalStat_l13_median512 finalStat_lopt_median512];
figure();
ax = axes();
hold(ax);
for i=1:numel(x)
    boxchart(x(i)*ones(size(b4_boxdata(:,i))), b4_boxdata(:,i), 'BoxFaceColor', b_colors(i,:), 'BoxFaceAlpha',0.7, 'MarkerColor',b_markerColor(i,:), 'MarkerStyle','.')
end
hold off
box on
set(gca, 'FontName','Arial','FontSize',24, 'LineWidth',3)
set(gcf,'RendererMode','manual')

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
