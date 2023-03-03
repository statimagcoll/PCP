---
title: Supplementary Material for "Small sample evaluation of resampling-based inference
  for topological features of neuroimages"
author: Simon Vandekar, Kaidi Kang, Neil Woodward, Anna Huang, Maureen McHugo, Shawn
  Garbett, Jeremy Stephens, Russell T. Shinohara, Armin Schwartzman, and Jeffrey Blume
date: "2/7/2020"
output:
  pdf_document: default
  html_document: default
---







# Simulation functions
















## Simulation results

Details of the simulation analyses and evaluation metrics are given in Section 4 of the paper. This section presents results for all simulations for parametric and robust test statistics evaluating the marginal and global CDFs.

## Sex covariate

The effect of sex was removed from the imaging data using equation (13), and then bootstrap samples of the residuals were modeled as a function of sex and tested on one degree of freedom.
Type 1 error rates and KL-divergence are given below.



```
## MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001
```

![Actual versus target type 1 error rates for the three inference procedures considered for testing the marginal distribution of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-2-1.pdf)

![Actual versus target type 1 error rates for the three inference procedures considered for testing the marginal distribution of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-3-1.pdf)


![KL-divergence for the three inference procedures considered for the marginal distribution of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-4-1.pdf)

![KL-divergence for the three inference procedures considered for the marginal distribution of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-5-1.pdf)





### Global distributions

![Actual versus target type 1 error rates for the three inference procedures considered for testing the distribution of the global maximum of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-6-1.pdf)

![Actual versus target type 1 error rates for the three inference procedures considered for testing the distribution of the global maximum of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-7-1.pdf)



![KL-divergence for the three inference procedures considered for the distribution of the global maximum of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-8-1.pdf)

![KL-divergence for the three inference procedures considered for the distribution of the global maximum of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-9-1.pdf)




## Fake group variable

Group was simulated independently of the imaging data and the bootstrap samples of the imaging data were modeled and tested on 3 degrees of freedom.
Type 1 error rates and KL-divergence are given below.


### Marginal distribution 


```
## MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001
```

![Actual versus target type 1 error rates for the three inference procedures considered for testing the marginal distribution of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-10-1.pdf)

![Actual versus target type 1 error rates for the three inference procedures considered for testing the marginal distribution of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-11-1.pdf)


![KL-divergence for the three inference procedures considered for the marginal distribution of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-12-1.pdf)

![KL-divergence for the three inference procedures considered for the marginal distribution of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-13-1.pdf)





### Global distributions

![Actual versus target type 1 error rates for the three inference procedures considered for testing the distribution of the global maximum of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-14-1.pdf)

![Actual versus target type 1 error rates for the three inference procedures considered for testing the distribution of the global maximum of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-15-1.pdf)



![KL-divergence for the three inference procedures considered for the distribution of the global maximum of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-16-1.pdf)

![KL-divergence for the three inference procedures considered for the distribution of the global maximum of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-17-1.pdf)








## Age continuous covariate fit with splines

The imaging data were residualized to age and other covariates using equation (13) in the main paper, then bootstrap samples of the residuals were modeled with age using splines on 4 degrees of freedom.
the test was for the nonlinear effect of age over the linear age effect on 3 degrees of freedom.


### Marginal distribution 


```
## MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001MaximaExtent; cft = 0.01Extent; cft = 0.001Mass; cft = 0.01Mass; cft = 0.001Global MaximaGlobal Extent; cft = 0.01Global Extent; cft = 0.001Global Mass; cft = 0.01Global Mass; cft = 0.001
```

![Actual versus target type 1 error rates for the three inference procedures considered for testing the marginal distribution of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-18-1.pdf)

![Actual versus target type 1 error rates for the three inference procedures considered for testing the marginal distribution of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-19-1.pdf)


![KL-divergence for the three inference procedures considered for the marginal distribution of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-20-1.pdf)

![KL-divergence for the three inference procedures considered for the marginal distribution of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-21-1.pdf)





### Global distributions

![Actual versus target type 1 error rates for the three inference procedures considered for testing the distribution of the global maximum of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-22-1.pdf)

![Actual versus target type 1 error rates for the three inference procedures considered for testing the distribution of the global maximum of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-23-1.pdf)



![KL-divergence for the three inference procedures considered for the distribution of the global maximum of each topological feature (TF) of the parametric test statistics image.](figure/unnamed-chunk-24-1.pdf)

![KL-divergence for the three inference procedures considered for the distribution of the global maximum of each topological feature (TF) of the robust test statistics image.](figure/unnamed-chunk-25-1.pdf)
