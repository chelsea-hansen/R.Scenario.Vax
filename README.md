
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R.Scenario.Vax

<!-- badges: start -->
<!-- badges: end -->

The goal of R.Scenario.Vax is to fit a deterministic, compartmental
transmission model to RSV hospitalization data and provide scenario
projections for the number of hospitalizations averted due to new
immunization products.

## Installation

You can install the development version of R.Scenario.Vax from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
remotes::install_github("chelsea-hansen/R.Scenario.Vax")
```

## Example

This is a basic example of how to prepare the population dynamics data
and initial values for compartments for New York state:

``` r
library(R.Scenario.Vax)


ny_data = get_data(state_or_county="state",state_abbr="NY",county_name=NULL)

#Parameter Values 
print(ny_data[[1]])
#> $PerCapitaBirthsYear
#>  [1] 4065.833    0.000    0.000    0.000    0.000    0.000    0.000    0.000
#>  [9]    0.000    0.000    0.000    0.000    0.000
#> 
#> $WidthAgeClassMonth
#>  [1]   2   2   2   2   2   2  12  36  60 120 240 240 240
#> 
#> $DurationMatImmunityDays
#> [1] 112
#> 
#> $RRHm
#> [1] 1
#> 
#> $recover1
#> [1] 182.625
#> 
#> $recover2
#> [1] 182.625
#> 
#> $recover3
#> [1] 358.9
#> 
#> $recover4
#> [1] 358.9
#> 
#> $um
#> [1] 0.0002148846
#> 
#> $rho1
#> [1] 0.75
#> 
#> $rho2
#> [1] 0.51
#> 
#> $dur.days1
#> [1] 10
#> 
#> $dur.days2
#> [1] 7
#> 
#> $dur.days3
#> [1] 5
#> 
#> $yinit.matrix
#>                  M0           S0          I1           R1           S1
#> <2m    2.327601e+04 1.219042e+04  5.20694384 2.235766e+02 3.078539e+02
#> 2-3m   1.493964e+04 1.887370e+04  6.76176381 8.467768e+02 1.295791e+03
#> 4-5m   9.476456e+03 2.159630e+04  7.08179310 1.823025e+03 2.971670e+03
#> 6-7m   5.933561e+03 2.189950e+04  6.86524825 2.907862e+03 4.977038e+03
#> 8-9m   3.681500e+03 2.104421e+04  6.38078909 3.851422e+03 6.910233e+03
#> 10-11m 2.279123e+03 1.978613e+04  5.89654601 4.520882e+03 8.558731e+03
#> 1Y     2.915675e+03 7.499557e+04 22.56024170 2.466448e+04 7.320475e+04
#> 2-4Y   7.848828e+02 7.581870e+04 25.22740145 3.229712e+04 1.731002e+05
#> 5-9Y   6.559234e+01 1.346010e+04 11.34628797 1.367530e+04 6.108295e+04
#> 10-19Y 3.237104e+00 1.128243e+03  1.11624320 2.574224e+03 1.053315e+04
#> 20-39Y 6.247554e-02 4.479497e+01  0.06247554 2.057320e+02 8.273637e+02
#> 40-59Y 0.000000e+00 7.211316e-01  0.00000000 5.547166e+00 2.579432e+01
#> 60Y+   0.000000e+00 0.000000e+00  0.00000000 1.582355e-01 1.107648e+00
#>                 I2           R2           S2          I3           R3
#> <2m     0.05367985     4.562787 5.582704e+00  0.00000000 5.367985e-02
#> 2-3m    0.21465922    21.948904 2.801303e+01  0.00000000 3.756537e-01
#> 4-5m    0.42919952    58.263847 7.827526e+01  0.00000000 1.233949e+00
#> 6-7m    0.75088649   116.280139 1.661605e+02  0.00000000 3.432624e+00
#> 8-9m    1.01878101   199.734732 3.059561e+02  0.00000000 8.043008e+00
#> 10-11m  1.23291439   315.840530 5.178778e+02  0.05360498 1.693917e+01
#> 1Y     11.17519269  8159.254473 2.139298e+04  1.83629927 2.788709e+03
#> 2-4Y   29.80430030 34238.107844 1.514069e+05 14.65695926 4.190637e+04
#> 5-9Y   27.23108671 32676.010521 1.125825e+05 28.31168451 9.394428e+04
#> 10-19Y  5.52540319 10125.385160 3.507755e+04 10.32525014 5.390071e+04
#> 20-39Y  0.49980440  1233.704782 4.448509e+03  1.56188886 1.152630e+04
#> 40-59Y  0.00000000    48.093928 2.098493e+02  0.05547166 7.543037e+02
#> 60Y+    0.00000000     1.793335 1.292256e+01  0.00000000 4.646848e+01
#>                  S3          I4           R4 Mn Mv N Si Vs1 Vs2
#> <2m    5.367985e-02   0.0000000 0.000000e+00  0  0 0  0   0   0
#> 2-3m   1.609944e-01   0.0000000 0.000000e+00  0  0 0  0   0   0
#> 4-5m   6.437993e-01   0.0000000 0.000000e+00  0  0 0  0   0   0
#> 6-7m   1.930851e+00   0.0000000 0.000000e+00  0  0 0  0   0   0
#> 8-9m   4.825807e+00   0.0000000 5.362007e-02  0  0 0  0   0   0
#> 10-11m 1.050658e+01   0.0000000 1.608149e-01  0  0 0  0   0   0
#> 1Y     3.381728e+03   0.2098628 3.044584e+02  0  0 0  0   0   0
#> 2-4Y   1.299883e+05   8.3909713 2.373174e+04  0  0 0  0   0   0
#> 5-9Y   5.073945e+05  85.2592329 2.705507e+05  0  0 0  0   0   0
#> 10-19Y 1.286831e+06 251.6570429 9.226240e+05  0  0 0  0   0   0
#> 20-39Y 2.929426e+06 670.8622703 2.432961e+06  0  0 0  0   0   0
#> 40-59Y 2.801750e+06 553.1078703 2.140549e+06  0  0 0  0   0   0
#> 60Y+   3.103731e+06 414.2076745 1.759833e+06  0  0 0  0   0   0
#> 
#> $q
#> [1] 1
#> 
#> $c2
#>            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,] 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.3840000
#>  [2,] 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.3840000
#>  [3,] 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.3840000
#>  [4,] 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.3840000
#>  [5,] 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.3840000
#>  [6,] 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.0640000 0.3840000
#>  [7,] 0.3840000 0.3840000 0.3840000 0.3840000 0.3840000 0.3840000 0.3840000
#>  [8,] 1.1520000 1.1520000 1.1520000 1.1520000 1.1520000 1.1520000 1.1520000
#>  [9,] 0.4066667 0.4066667 0.4066667 0.4066667 0.4066667 0.4066667 0.4066667
#> [10,] 0.3710833 0.3710833 0.3710833 0.3710833 0.3710833 0.3710833 0.3710833
#> [11,] 1.5046667 1.5046667 1.5046667 1.5046667 1.5046667 1.5046667 1.5046667
#> [12,] 0.6354167 0.6354167 0.6354167 0.6354167 0.6354167 0.6354167 0.6354167
#> [13,] 0.3431250 0.3431250 0.3431250 0.3431250 0.3431250 0.3431250 0.3431250
#>            [,8]      [,9]     [,10]    [,11]     [,12]    [,13]
#>  [1,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [2,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [3,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [4,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [5,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [6,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [7,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [8,] 1.1520000 0.4066667 0.3710833 1.504667 0.6354167 0.343125
#>  [9,] 0.4066667 6.6400000 1.7350000 3.355000 1.8550000 1.020000
#> [10,] 0.3710833 1.7350000 8.0550000 2.727500 2.5675000 1.162500
#> [11,] 1.5046667 3.3550000 2.7275000 4.907500 3.1775000 1.505000
#> [12,] 0.6354167 1.8550000 2.5675000 3.177500 3.8225000 2.300000
#> [13,] 0.3431250 1.0200000 1.1625000 1.505000 2.3000000 3.300000
#> 
#> $sigma1
#> [1] 0.89
#> 
#> $sigma2
#> [1] 0.7209
#> 
#> $sigma3
#> [1] 0.237897
#> 
#> $length.step
#> [1] 7
#> 
#> $time.step
#> [1] "week"
#> 
#> $seed
#> [1] 196.7715
#> 
#> $waningN
#> [1] 150
#> 
#> $waningV
#> [1] 180
#> 
#> $waningS
#> [1] 730.5
#> 
#> $RRIn
#> [1] 1
#> 
#> $RRIv
#> [1] 1
#> 
#> $RRIs
#> [1] 1
#> 
#> $RRHn
#> [1] 0.2
#> 
#> $RRHv
#> [1] 0.45
#> 
#> $RRHs
#> [1] 0.2
```
