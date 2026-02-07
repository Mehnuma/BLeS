
# HHJ Example
% Example:
% 	y_hhj = readmatrix('example_dataset.csv'); 
% 	hhj_lhat_var = HHJ_bles(y_hhj,5,30,'mean','variance',500,0,0,[]);
% 	hhj_lhat_dist = HHJ_bles(y_hhj,5,30,'mean','one-sided',500,7.8960e-04,0,[]);


# BK
% Example:
% 	y_bk = readmatrix('example_dataset.csv');
%	bk_lhat_mean = BK_bles(y_bk,'mean');
%	bk_lhat_median = BK_bles(y_bk,'median');


# cPW
% Example:
%	y_cpw = readmatrix('example_dataset.csv');
%	cpw_lhat_cbb = cPW_bles(y_cpw, 'circularBB');
%	cpw_lhat_sbb = cPW_bles(y_cpw, 'stationaryBB');

# NPPI
% Example:
%	y_nppi = readmatrix('example_dataset.csv');
%	[nppi_lhat_var, nppi_stat_val] = NPPI_bles(y_nppi,'mean','variance',500,[],[],[]);
%	[nppi_lhat_var, nppi_stat_val] = NPPI_bles(y_nppi,'mean','quantile',500,7.8960e-04,[],0.35);

# TBB/ETBB
% Example:
%	y_tbb = readmatrix('example_dataset.csv');
%	tbb_lhat = TBB_bles(y_tbb,'TBB');
%	etbb_lhat = TBB_bles(y_tbb,'ETBB'); 

# SETBB
% Example:
%	qData = readmatrix('qregData.csv');
%	X = qdata(:,1:2); Y = qData(:,3);
%	setbb_lhat = SETBB_bles(X,Y,0.35,'SETBB',500); 
%	setbb_lhat = SETBB_bles(X,Y,0.85,'SETBB',750); 


# BD
% Example:
%	pData = readmatrix('pData.csv');
% 	y_bd = pData(1:200,1);
%	bd_lhat_overall = BD_bles(y_bd,6,'overall mean','EMBB');
%	bd_lhat_seasonal = BD_bles(y_bd,6,'seasonal mean','GSBB');
