
# Nearest-neighbor Projected-Distance Regression (NPDR)

NPDR is a nearest-neighbor feature selection algorithm that fits a generalized linear model for projected distances of a given attribute over all pairs of instances in a neighborhood. In the NPDR model, the predictor is the attribute distance between neighbors projected onto the attribute dimension, and the outcome is the projected phenotype distance (for quantitative traits) or hit/miss (for case/control) between all pairs of nearest neighbor instances. NPDR can fit any combination of predictor data types (catgorical or numeric) and outcome data types (case-control or quantitative) as well as adjust for covariates that may be confounding. As with STIR (STatistical Inference Relief), NDPR allows for the calculation of statistical significance of importance scores and adjustment for multiple testing.   

#### Websites

[NPDR Github Page](https://insilico.github.io/npdr/)

[STIR Github Page](https://insilico.github.io/stir/)

[insilico Github Organization](https://github.com/insilico)

[insilico McKinney Lab](http://insilico.utulsa.edu/)

#### Related References. 

[2017 STIR paper in Bioinformatics](https://doi.org/10.1093/bioinformatics/bty788)

[2013 Gene-Wise Adaptive-Neighbors paper in PLoS One](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081527)

### To install:

    >library(devtools)
    
    >install_github("insilico/npdr")  # todo (optional build_vignettes = TRUE)
    >library(npdr)
    >data(package="npdr")
    
    # >vignette(" ") # todo (if you build_vignettes)

### Dependencies
To install the `tidyverse` collection of R packages:

```
install.packages('tidyverse')
```
To set `fast.reg = TRUE` or `fast.dist = TRUE` or `use.glmnet = TRUE`, please install the `speedglm` and `glmnet` packages:

```
install.packages(c('speedglm', 'wordspace', 'glmnet'))
```

If an issue arises with updating `openssl`, try updating it on your own system, e.g. for MacOS:

```brew install openssl@1.1```

### Examples


### Abstract

Relief-based methods are nearest-neighbor machine learning feature selection algorithms that compute the importance of attributes that may involve interactions in high-dimensional data. Previously we introduced STIR, which extended Relief-based methods to compute statistical significance of attributes in case-control data by reformulating the Relief score as a pseudo t-test. Here we extend the statistical formalism of STIR to a generalized linear model (glm) formalism to handle quantitative and case-control outcome variables, any predictor data type (continuous or categorical), and  adjust  for  covariates  while  computing statistical significance of attributes.


#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
