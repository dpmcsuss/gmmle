Readme
================
Daniel Sussman
5/13/2020

## GMMLE Package

A package for running simulations and creating plots found in [Maximum
Likelihood Estimation and Graph Matching in Errorfully Observed
Networks](https://arxiv.org/abs/1812.10519).

``` r
devtools::install_github("dpmcsuss/gmmle")
```

``` r
library(gmmle)
```

``` r
plot_prob_bound_norm_disag()
```

![](Readme_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plot_matching_errors()
```

![](Readme_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
plot_cum_acc_by_norm_deg()
```

    ## Computing cumaltive errors for each match. This can take a few minutes

![](Readme_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
make_er_plot()
```

![](Readme_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
make_ws_plot()
```

![](Readme_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
make_model_compare_plot()
```

![](Readme_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->
