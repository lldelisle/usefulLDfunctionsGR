# usefulLDfunctionsGR
Functions That I Use Regularly Which depends on GenomicRanges

## Installation
This package depends on GenomicRanges, rtracklayer, usefulLDfunctions, stats, grDevices, utils, combinat.
If usefulLDfunctions, and/or GenomicRanges and/or rtracklayer and/or combinat are not installed on your R see Installation of dependencies section.

The easiest way to install usefulLDfunctions is using devtools::install_github() from R:
```
if (!"devtools" %in% installed.packages()){
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctionsGR")
```

## Installation of dependencies
usefulLDfunctions is another R package I wrote which is on github too. So to install it:
```
if (!"devtools" %in% installed.packages()){
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
```
As the installation of Bioconductor package depends on the R version you have, I recommend you to use the function `safelyLoadAPackageInCRANorBioconductor` from the usefulLDfunctions package:
```
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
safelyLoadAPackageInCRANorBioconductor("combinat")
```

## Issues
If you have issues, use the Issues in github or send an email to lucille.delisle\@epfl.ch
