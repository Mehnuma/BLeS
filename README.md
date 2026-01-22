# BLeS: A MATLAB and Octave Toolbox for Block Length Selection in Block Bootstrap.

We present a MATLAB and Octave toolbox, `BLeS`, for selecting the block length, the main tuning parameter in the block bootstrap (a resampling method for estimating the empirical distribution of dependent data). This package enables users to select block lengths using various methods provided by different scholars for several well-known block bootstrap approaches. This toolbox offers the following implementations for block length selection:
1. Hall-Horowitz-Jing (HHJ) method
2. Bühlmann-Künsch (BK) method
3. Corrected Politis-White (cPW) method
4. Nonparametric plugin (NPPI) method by Lahiri et al.
5. Methods for tapered block bootstrap variants:
     (a) Original tapered block bootstrap (TBB) by Paparoditis-Politis
     (b) Extended tapered block bootstrap (ETBB) by Shao
     (c) Smooth extended TBB (SETBB) by Gregory et al.
6. Bertail-Dudek (BD) Method

This toolbox has been tested using different examples from published articles, with the replication scripts available in the `Replication Codes` folder. BLeS is also supported in Octave with the codes placed in the `Octave` folder. Additionally, this repository also contains the codes for the Simulation section in Tabassum and De Brabanter (2026) in the `Replication Codes` folder.


# References
[1] Hall, P., Horowitz, J. L., & Jing, B. Y. (1995). On blocking rules for the bootstrap with dependent data. Biometrika, 82(3), 561-574.

[2] Bühlmann, P., & Künsch, H. R. (1999). Block length selection in the bootstrap for time series. Computational Statistics & Data Analysis, 31(3), 295-310.

[3] Politis, D. N., & White, H. (2004). Automatic block-length selection for the dependent bootstrap. Econometric reviews, 23(1), 53-70.

[4] Patton, A., Politis, D. N., & White, H. (2009). Correction to “Automatic block-length selection for the dependent bootstrap” by D. Politis and H. White. Econometric Reviews, 28(4), 372-375.

[5] Lahiri, S. N., Furukawa, K., & Lee, Y. D. (2007). A nonparametric plug-in rule for selecting optimal block lengths for block bootstrap methods. Statistical methodology, 4(3), 292-321.

[6] Paparoditis, E., & Politis, D. (2002). The tapered block bootstrap for general statistics from stationary sequences. The Econometrics Journal, 5(1), 131-148.

[7] Shao, X. (2010). Extended tapered block bootstrap. Statistica Sinica, 807-821.

[8] Gregory, K. B., Lahiri, S. N., & Nordman, D. J. (2018). A smooth block bootstrap for quantile regression with time series. The Annals of Statistics, 46(3), 1138-1166.

[9] Bertail, P., & Dudek, A. E. (2024). Optimal choice of bootstrap block length for periodically correlated time series. Bernoulli, 30(3), 2521-2545.

[10] Tabassum, M., & De Brabanter, K. (2025).  A State-of-the-Science Overview of Block Length Selection Methods in Block Bootstrap. Statistical Science, In review.

