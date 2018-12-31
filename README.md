
# glmSTIR: generalized linear model (GLM) based STatistical Inference Relief (STIR)

glmSTIR is a nearest-neighbor feature selection algorithm that fits a generalized linear model for each attribute over all pairs of instances in a neighborhood. In the glmSTIR model, the predictor is the attribute distance between neighbors and the outcome is the phenotype distance (for quantitative traits) or hit/miss (for case/control) between all pairs of nearest neighbor instances. glmSTIR can fit any combination of predictor data types (catgorical or numeric) and outcome data types (case-control or quantitative) as well as adjust for covariates that may be confounding. As with STIR, glmSTIR allows for the calculation of statistical significance of importance scores and adjustment for multiple testing.   

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

Relief-based methods are nearest-neighbor machine learning feature selection algorithms that compute the importance of attributes that may involve interactions in high-dimensional data. Previously we introduced STIR, which extended Relief-based methods to compute statistical significance of attributes in case-control data by reformulating the Relief score as a pseudo t-test. Here we extend the statistical formalism of STIR to a generalized linear model (glm) formalism to handle quantitative and case-control outcome variables, any predictor data type (continuous or categorical), and  adjust  for  covariates  while  computing statistical significance of attributes.

#### Websites
[https://github.com/insilico](https://github.com/insilico)

#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
