---
title: 'Block Length Selection for Block Bootstrap in MATLAB and Octave using BLeS'
tags:
  - Block length selection
  - Block bootstrap
  - MATLAB
  - Octave
  - BLeS
authors:
  - name: Mehnuma Tabassum
    orcid: 0000-0001-6057-6763
    affiliation: 1
  - name:
      given: Kris
      non-dropping-particle: De
      family: Brabanter
    orcid: 0000-0002-7697-1079
    affiliation: "2, 1"
affiliations:
 - name: Department of Industrial and Manufacturing Systems Engineering, Iowa State University, Ames, IA 50011
   index: 1
 - name: Department of Statistics, Iowa State University, Ames, IA 50011
   index: 2
date: 9 January, 2026
bibliography: paper.bib

---

# Summary
We present a MATLAB toolbox called BLeS for selecting block length, the primary tuning parameter in block bootstrap, a resampling method for estimating an empirical distribution of dependent data. This package enables users to select block lengths using various methods provided by different scholars for several well-known block bootstrap approaches. The implementations include the Hall-Horowitz-Jing method, Bühlmann-Künsch method, corrected Politis-White method, and the nonparametric plug-in method. This toolbox has been tested using different examples from published articles, and a simulation study illustrating its range and functionality has been presented. BLeS is also supported in Octave.

# Statement of Need
Bootstrapping is a resampling method used to approximate properties of an estimator, such as bias, variance, and its distribution. Standard bootstrap assumes independent and identically distributed (IID) data [@efron1993], and hence it is unsuitable for correlated observations [@kunsch1989]. To handle dependence, block bootstrap resamples blocks of consecutive observations from a correlated series instead of individual IID points[@hall1995]. Block bootstrap has several variants—fixed, moving/overlapping, circular, stationary, etc.—but in all cases the key issue is choosing the block length[@lahiri1999], which largely determines performance. Thus, block length is the main tuning parameter and its selection is an important research problem.

There exist several mean-squared error (MSE)-optimal block length in some cases; however, they cannot be obtained in practice due to including population parameters. As a result, the available literature focuses on empirical block length choices for a certain block bootstrap variant. Pioneering articles include @hall1995, @buhlmann1999, @paparoditis2002, @politis2004, @lahiri2007,  @shao2010, @gregory2018, @bertail2024. For a detailed overview and their comparative capabilities, we refer the readers to @tabassum2026. 

Because these methods have been developed for practical purposes, their implementations and the ease-of-use is of prime importance. Although some free-source codes for block bootstrap can be found (see [here](https://github.com/lamferzon/bboot) for bboot in R, and [here](https://www.mathworks.com/matlabcentral/fileexchange/53701-bootstrapping-time-series) for some MATLAB codes), there are only a few sources that implement block length selection algorithms. We have identified some sources that, in some capacity, have implemented some of the methods available in the literature. The most frequently implemented ones are by @politis2004 and @hall1995, while there are almost nonexistent sources for the other methods, to the best of our knowledge


# State of the Field 
The available implementations are presented below:

| Package / Function(s)        | Author(s)         | Year | Capabilities                                                                             | Environment       | Source                                                                                                             |
| ---------------------------- | ----------------- | ---- | ---------------------------------------------------------------------------------------- | ----------------- | ------------------------------------------------------------------------------------------------------------------ |
| `opt_block_length_REV_dec07` | @patton2007       | 2007 | Block length selector for CBB and SBB                                                    | MATLAB            | [Andrew Patton's Matlab Page](https://public.econ.duke.edu/~ap172/code.html)                                       |
| `b.star`                     | @hayfield2008     | 2008 | Implements block length selection for CBB and SBB                                        | R                 | [np package](https://CRAN.R-project.org/package=np)                                                                |
| `getNPPIblksizesQR`          | @gregory2022      | 2022 | Implements block length selection for MBB, SMBB, ETBB, and SETBB for quantile regression | R                 | [QregBB package](https://cran.r-project.org/package=QregBB)                                                        |
| `optimal_block_length`       | @nowotny2019      | 2019 | Block length selectors for CBB and SBB                                                   | Python            | [recombinator package](https://github.com/InvestmentSystems/recombinator)                                          |
| `optimal_block_length`       | @sheppard2021     | 2021 | Block length selectors for CBB and SBB                                                   | Python            | [arch package](https://arch.readthedocs.io/en/latest/bootstrap/generated/arch.bootstrap.optimal_block_length.html) |
| `OBL`                        | @james2022        | 2022 | Calculates optimal block length for NBB, MBB, CBB, TMBB, and TCBB                        | R                 | [OBL package](https://CRAN.R-project.org/package=OBL)                                                              |
| `blocklength`                | @stashevsky2025   | 2025 | Implements HHJ, cPW, and NPPI                                                            | R                 | [blocklength package](https://cran.r-project.org/package=blocklength)                                              |
| `boodd`                      | @bertail2025      | 2025 | Calculates optimal block length for SB and BD                                            | R                 | [boodd package](https://cran.r-project.org/web/packages/boodd/index.html)                                          |
| `BLeS`                       | @tabassum2026     | 2026 | Implements HHJ, BK, cPW, NPPI, TBB/ETBB, SETBB/modified SETBB, and BD                    | MATLAB and Octave | [BLeS Toolbox](https://github.com/Mehnuma/BLeS)                                                                    |

Notes:
BB: Block bootstrap; CBB: Circular BB; SBB: Stationary BB; MBB: Moving BB; SMBB: Smooth MBB; ETBB: Extended tapered BB; SETBB: Smooth ETBB; HHJ: Hall–Horowitz–Jing; cPW: corrected Politis–White; NPPI: Nonparametric plug-in; NBB: Non-overlapping BB; TMBB: Tapered MBB; TCBB: Tapered CBB; TBB: Tapered BB.

# Software Design
The toolbox is designed as follows. It includes the implements for the HHJ method (\autoref{fig:hhj}), BK method (\autoref{fig:bk}), cPW method (\autoref{fig:cpw}), NPPI method (\autoref{fig:nppi}), the block length selectors for the tapered block bootstrap variants, and the BD method. 

# Research Impact Statement
Given that the resources for block length selection methods are not available in many cases, and the available ones are scattered in a manner, we propose a MATLAB and Octave toolbox called BLeS that contains the implementations of all the aforementioned methods in one place. The unique characteristics of this toolbox are that it:

(1) Provides a complete implementation scheme for almost every block length selection method available in the literature.

(2) Offers ease-of-use through only one line of code and the opportunity to customize the inputs.

(3) Contains both MATLAB and Octave implementations to provide more flexibility to the users.

# References
<div id="refs"></div>

# Appendix
In this section, we present the graphical representations of the block length selectors in BLeS. For details of these methods, the readers are referred to the source articles, and to @tabassum2026 for a comprehensive overview. 

![Hall-Horowitz-Jing (HHJ) Method of Block Length Selection.\label{fig:hhj}](bles_hhj.pdf)

![Bühlmann-Künsch (BK) Method of Block Length Selection.\label{fig:bk}](bles_buhlmann.pdf)

![Corrected Politis-White (cPW) Method of Block Length Selection.\label{fig:cpw}](bles_politis.pdf)

![Nonparametric plug-in (NPPI) Method of Block Length Selection.\label{fig:nppi}](bles_lahiri.pdf)
