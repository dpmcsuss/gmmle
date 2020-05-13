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

Load the
    library.

``` r
library(gmmle)
```

### Figure 1

Left:

``` r
make_er_plot()
```

![](Readme_files/figure-gfm/er_plot-1.png)<!-- -->

Right:

``` r
make_ws_plot()
```

![](Readme_files/figure-gfm/ws_plot-1.png)<!-- -->

### Figure 2

``` r
make_model_compare_plot()
```

![](Readme_files/figure-gfm/model_compare_plot-1.png)<!-- -->

### Figure 3

``` r
plot_prob_bound_norm_disag()
```

![](Readme_files/figure-gfm/prob_bound_norm_disag-1.png)<!-- -->

### Figure 4

``` r
plot_matching_errors()
```

![](Readme_files/figure-gfm/matching_errors-1.png)<!-- -->

### Figure 5

``` r
plot_cum_acc_by_norm_deg()
```

    ## Computing cumaltive errors for each match. This can take a few minutes

![](Readme_files/figure-gfm/cum_acc_by_norm_deg-1.png)<!-- -->
