
# A Mixed Effect Similarity Matrix Regression Model (SMRmix) for Integrating Multiple Microbiome Datasets at Community Level


## Overview

This package provides implementation of SMRmix, which integrates heterogeneous microbiome datasets and tests the association between microbial composition and an outcome of interest.

## Required packages 

Required packages for functions in MIDASim include: stats, harmonicmeanp, lme4, nlme, MASS, expm. These packages are all available on CRAN.

## Installation

A github installation is available via

```{r}
install_github("mengyu-he/SMRmix")
```

Install from the local file (the .tar.gz file), which you can download from https://github.com/mengyu-he/SMRmix .

```{r}
devtools::install_local("SMRmix_0.1.0.tar.gz", dependencies = TRUE)
```


## Usage

### Example data description

The example data can be directly loaded through 

```{r}
data("example_SMRmix")
attach(example_SMRmix)
```

The example dataset is built upon two microbiome studies [1,2]. The sample sizes are 23, 20, respectively. The datasets consists of two objects:

1. The kernel matrices calculated based on 4 different kernel metrics: weighted UniFrac, unweighted UniFrac, Jaccard, Bray-Curtis. 

2. The metadata including HIV and age information, where HIV is a binary status (1-positive, 0-negative), and age is a continuous variable.

In the following, with the example dataset, we demonstrate the usage of main function *SMRmix()* for continuous outcomes and binary outcomes.

### Binary outcome

Suppose the outcome of interest is HIV status, and the goal is to test the association between microbiome composition and HIV status adjusting for age,
$$logit(P(HIV_{ik}=1)) = \beta \cdot age_{ik}+ f(\boldsymbol{M_{ik}})+h_k,$$
where $\boldsymbol{M_{ik}}$ represents the microbiome profile for $i$-th subject in $k$-th study, and $f(\boldsymbol{M_{ik}})$ depicts the effect of microbiome composition on HIV status, $h_k$ is a study-specific random effect. Then SMRmix can be applied by

```{r}
SMRmix( formula =  HIV ~ 1 + Age + (1 | Study.id), Kernels = Kernels, outcome = "binary")
```

### Continous outcome

Suppose the outcome of interest is age, and the goal is to test the association between microbiome composition and age adjusting for HIV status,
$$age_{ik}= \beta \cdot HIV_{ik}+ f(\boldsymbol{M_{ik}})+h_k,$$
where $\boldsymbol{M_{ik}}$ represents the microbiome profile for $i$-th subject in $k$-th study, and $f(\boldsymbol{M_{ik}})$ depicts the effect of microbiome composition on HIV status, $h_k$ is a study-specific random effect. Then SMRmix can be applied by

```{r}
SMRmix( formula =  Age ~ 1 + HIV + (1 | Study.id), Kernels = Kernels, outcome = "continuous")
```

## References

[1] Dillon, S. M., Lee, E. J., Kotter, C. V., Austin, G. L., Dong, Z., Hecht, D. K., Gianella, S., Siewe, B., Smith, D. M., Landay, A. L., Robertson, C. E., Frank, D. N., & Wilson, C. C. (2014). An altered intestinal mucosal microbiome in HIV-1 infection is associated with mucosal and systemic immune activation and endotoxemia. Mucosal immunology, 7(4), 983–994. https://doi.org/10.1038/mi.2013.116

[2] Dinh, D. M., Volpe, G. E., Duffalo, C., Bhalchandra, S., Tai, A. K., Kane, A. V., Wanke, C. A., & Ward, H. D. (2015). Intestinal microbiota, microbial translocation, and systemic inflammation in chronic HIV infection. The Journal of infectious diseases, 211(1), 19–27. https://doi.org/10.1093/infdis/jiu409

