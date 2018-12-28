
# glmSTIR: generalized linear model (GLM) based STatistical Inference Relief (STIR)

glmSTIR is a nearest-neighbor feature selection algorithm that fits generalized linear models of pairwise distances between instances. glmSTIR can fit any combination of predictor data types (catgorical or numeric) and outcome data types (case-control or quantitative) as well as control for covariates. As with STIR, glmSTIR allows for the calulation of statistical significance of importance scores.   

Information about the original STIR algorithm for case-control outcome data can be found here. 
Trang T. Le, Ryan J. Urbanowicz, Jason H. Moore, B. A McKinney. “STatistical Inference Relief (STIR) feature selection,” Bioinformatics. 18 September 2018. [https://doi.org/10.1093/bioinformatics/bty788](https://doi.org/10.1093/bioinformatics/bty788)

[http://insilico.utulsa.edu/software/STIR](http://insilico.utulsa.edu/software/STIR)

### To install:

    >library(devtools)
    
    >install_github("insilico/glmSTIR")  # todo (optional build_vignettes = TRUE)
    >library(glmSTIR)
    >data(package="glmSTIR")
    
    # >vignette(" ") # todo (if you build_vignettes)
    
### Examples


### Abstract

Relief-based methods are nearest-neighbor machine learning feature selection algorithms that compute the importance of attributes that may involve interactions in high-dimensional data. Previously we introduced STIR, which extended Relief-based methods to compute statistical significance of attributes in case-control  data  by reformulating the Relief score as a pseudo t-test. Here we extend the statistical formalism of STIR to a generalized linear model (glm) formalism to handle quantitative and case-control outcome variables, any predictor data type (continuous or categorical), and  correct  for  covariates  while  computing statistical significance of attributes.

#### Websites
[https://github.com/insilico](https://github.com/insilico)

#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
