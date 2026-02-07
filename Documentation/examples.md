# BLeS Examples
Below are code snippets for the functions in the `BLeS` toolbox. These examples use the following datasets: `example_dataset.csv`, `qregData.csv`, and `pData.csv`, which are available in the [Replication Codes Folder](Replication Codes/BLeS Replication Codes).

## HHJ Example
For the demonstration of the HHJ method, we use the following examples. The first example below estimates the block length for the bootstrap variance estimation of the sample mean using the HHJ method. The second example selects the block length for the one-sided distribution estimation of the sample mean. In both cases, a pilot block length of size 5 and a subsample size of 30 is used, in accordance with Lahiri (2003). The "theta_hat" value was obtained via additional simulations for model (6.1) in Lahiri et al. (2007).
```matlab
% Example:
 	y_hhj = readmatrix('example_dataset.csv'); 
 	hhj_lhat_var = HHJ_bles(y_hhj,5,30,'mean','variance',500,[],[],[]);
 	hhj_lhat_dist = HHJ_bles(y_hhj,5,30,'mean','one-sided',500,7.8960e-04,0,[]);
```

## BK Example
The BK method uses the sample mean or median as the statistic and the variance functional to approximate the optimal block length. The usage is pretty straightforward, as follows:
```matlab
% Example:
  y_bk = readmatrix('example_dataset.csv');
  bk_lhat_mean = BK_bles(y_bk,'mean');
  bk_lhat_median = BK_bles(y_bk,'median');
```

## cPW Example
Like the BK method, the cPW method uses the variance functional only for the sample mean. The input arguments are the data and the block bootstrap variant.
```matlab
% Example:
	y_cpw = readmatrix('example_dataset.csv');
	cpw_lhat_cbb = cPW_bles(y_cpw, 'circularBB');
	cpw_lhat_sbb = cPW_bles(y_cpw, 'stationaryBB');
```

## NPPI Example
The NPPI method can be used to estimate bias, variance, distribution, and quantiles for linear statistics. Below, the first example selects the optimal block length for variance estimation of the sample mean, while the second example concerns the 35% quantile of a studentized statistic; see Lahiri (2007, Eq. 4.5).
```matlab
% Example:
	y_nppi = readmatrix('example_dataset.csv');
	[nppi_lhat_var, nppi_stat_val] = NPPI_bles(y_nppi,'mean','variance',500,[],[],[]);
	[nppi_lhat_var, nppi_stat_val] = NPPI_bles(y_nppi,'mean','quantile',500,7.8960e-04,[],0.35);
```

## TBB/ETBB Example
The following example estimates for the TBB/ETBB method for the example dataset.
```matlab
% Example:
	y_tbb = readmatrix('example_dataset.csv');
	tbb_lhat = TBB_bles(y_tbb,'TBB');
	etbb_lhat = TBB_bles(y_tbb,'ETBB'); 
```

## SETBB Example
The block length selection for the quantile regression estimator requires the input-output pair and the desire quantile. Additionally, the `boot_variant` needs to be specified to `SETBB` or `MBB` to specify how the bootstrap estimates are calculated. 
```matlab
% Example:
	qData = readmatrix('qregData.csv');
	X = qdata(:,1:2); Y = qData(:,3);
	setbb_lhat = SETBB_bles(X,Y,0.35,'SETBB',500); 
	setbb_lhat = SETBB_bles(X,Y,0.85,'SETBB',750); 
```

## BD Example
The following examples implement the BD method by Bertail and Dudek (2024), for the `overall mean` using the `EMBB` method, and for the `seasonal mean` using the `GSBB` method, respectively. 
```matlab
% Example:
	pData = readmatrix('pData.csv');
 	y_bd = pData(1:200,1);
	bd_lhat_overall = BD_bles(y_bd,6,'overall mean','EMBB');
	bd_lhat_seasonal = BD_bles(y_bd,6,'seasonal mean','GSBB');
```
