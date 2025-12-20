# BLeS: A MATLAB and Octave Toolbox for Block Length Selection in Block Bootstrap.

We present a MATLAB toolbox called BLeS for selecting block length, the main tuning parameter in block bootstrap, which is a resampling method for estimating an empirical distribution for dependent data. This package enables users to select block lengths using various methods provided by different scholars for several well-known block bootstrap approaches. This toolbox offers the following implementations for block length selection:
1. Hall-Horowitz-Jing (HHJ) method
2. Bühlmann-Künsch (BK) method
3. Corrected Politis-White (cPW) method
4. Nonparametric plugin (NPPI) method by Lahiri et al.
5. Methods for tapered block bootstrap variants:
     (a) Original tapered block bootstrap (TBB) by Paparoditis-Politis
     (b) Extended tapered block bootstrap (ETBB) by Shao
     (c) Smooth extended TBB (SETBB) by Gregory et al. 

This toolbox has been tested using different examples from published articles, and a simulation study illustrating its range and functionality has been presented. BLeS is also supported in Octave.
