---
title: Supplementary Material for "Small sample evaluation of resampling-based longitudinal inference
  for topological features of neuroimages"
author: Simon Vandekar
date: "8/3/2022"
output:
  pdf_document: default
  html_document: default
---





[1] 51925







# Simulation functions







```
## Error: <text>:60:1: unexpected '}'
## 59:   unlink(simdirs)
## 60: }
##     ^
```









## Simulation results

Details of the simulation analyses and evaluation metrics are given in Section 4 of the paper. This section presents results for all simulations for parametric and robust test statistics evaluating the marginal and global CDFs.

## Sex covariate

The effect of sex was removed from the imaging data using equation (13), and then bootstrap samples of the residuals were modeled as a function of sex and tested on one degree of freedom.
Type 1 error rates and QQ-plot are given below.















### Global distributions














## Fake group variable

Group was simulated independently of the imaging data and the bootstrap samples of the imaging data were modeled and tested on 3 degrees of freedom.
Type 1 error rates and QQ-plot are given below.


### Marginal distribution 


```
## Error in load(rdata): object 'fakeGroupSimConfig' not found
```

```
## Error in eval(expr, envir, enclos): object 'fakeGroupSimConfig' not found
```

```
## Error in errorPlots(pd, nsim = nsim, tStat = TRUE, global = FALSE): object 'pd' not found
```


```
## Error in errorPlots(pd, nsim = nsim, tStat = FALSE, global = FALSE): object 'pd' not found
```



```
## Error in qqPlots(pd, nsim = nsim, tStat = TRUE, global = FALSE): object 'pd' not found
```


```
## Error in qqPlots(pd, nsim = nsim, tStat = FALSE, global = FALSE): object 'pd' not found
```





### Global distributions


```
## Error in errorPlots(pd, nsim = nsim, tStat = TRUE, global = TRUE): object 'pd' not found
```


```
## Error in errorPlots(pd, nsim = nsim, tStat = FALSE, global = TRUE): object 'pd' not found
```




```
## Error in qqPlots(pd, nsim = nsim, tStat = TRUE, global = TRUE): object 'pd' not found
```


```
## Error in qqPlots(pd, nsim = nsim, tStat = FALSE, global = TRUE): object 'pd' not found
```








## Age continuous covariate fit with splines

The imaging data were residualized to age and other covariates using equation (13) in the main paper, then bootstrap samples of the residuals were modeled with age using splines on 4 degrees of freedom.
the test was for the nonlinear effect of age over the linear age effect on 3 degrees of freedom.


### Marginal distribution 


```
## Error in load(rdata): object 'ageSplineSimConfig' not found
```

```
## Error in eval(expr, envir, enclos): object 'ageSplineSimConfig' not found
```

```
## Error in errorPlots(pd, nsim = nsim, tStat = TRUE, global = FALSE): object 'pd' not found
```


```
## Error in errorPlots(pd, nsim = nsim, tStat = FALSE, global = FALSE): object 'pd' not found
```



```
## Error in qqPlots(pd, nsim = nsim, tStat = TRUE, global = FALSE): object 'pd' not found
```


```
## Error in qqPlots(pd, nsim = nsim, tStat = FALSE, global = FALSE): object 'pd' not found
```





### Global distributions


```
## Error in errorPlots(pd, nsim = nsim, tStat = TRUE, global = TRUE): object 'pd' not found
```


```
## Error in errorPlots(pd, nsim = nsim, tStat = FALSE, global = TRUE): object 'pd' not found
```




```
## Error in qqPlots(pd, nsim = nsim, tStat = TRUE, global = TRUE): object 'pd' not found
```


```
## Error in qqPlots(pd, nsim = nsim, tStat = FALSE, global = TRUE): object 'pd' not found
```




