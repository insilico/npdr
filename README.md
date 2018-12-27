
# generalized linear model (GLM) based STatistical Inference Relief (STIR)

glmSTIR is a nearest-neighbor feature selection algorithm fits generalized linear models of pairwise distances. glmSTIR can fit any combination of predictor (catgorical or numeric) and outcome (case-control or quantitative) data types and correct for covariates. As with STIR, glmSTIR allows for the calulation of statistical significance of importance scores.   

Information about the original STIR algorithm for case-control outcome data can be found here. 
Trang T. Le, Ryan J. Urbanowicz, Jason H. Moore, B. A McKinney. “STatistical Inference Relief (STIR) feature selection,” Bioinformatics. 18 September 2018. [https://doi.org/10.1093/bioinformatics/bty788](https://doi.org/10.1093/bioinformatics/bty788)

[http://insilico.utulsa.edu/software/STIR](http://insilico.utulsa.edu/software/STIR)

Example STIR usage and output: [STIRexample.md](https://github.com/insilico/glmSTIR/blob/master/inst/example/STIRexample.md).

### To install:

    >library(devtools)
    
    >install_github("insilico/glmSTIR")  # you can use build_vignettes = TRUE but slows down install

    >library(glmSTIR)
    
    >data(package="glmSTIR")
    
    # >vignette("STIRsimulated") # if you build_vignettes
    
    # >vignette("STIRmdd")

    
 ### Examples

[Simulated Data Example with privateEC Simulation](https://github.com/insilico/glmSTIR/blob/master/inst/example/STIRexample.md) 

[RNA-Seq Example](https://github.com/insilico/glmSTIR/blob/master/vignettes/STIRmdd.Rmd) 

### Abstract

Relief-based methods are nearest-neighbor machine learning feature selection algorithms that compute the importance of attributes that may involve interactions in high-dimensional data. Previously we introduced STIR, which extended Relief-based methods to compute statistical significance of attributes  in  case-control  data  by reformulating the Relief score as a pseudo t-testHere we extend the statistical formalism of STIR to regression to handle continuous and case-control outcome variables, any predictor data type (continuous  or  categorical),  and  adjust  for  covariates  while  computing statistical significance of attributes.

#### Websites
[https://github.com/insilico](https://github.com/insilico)

#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
