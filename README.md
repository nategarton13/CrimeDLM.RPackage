# CrimeDLM_RPackage

In order to download and install this package, first install the devtools R package. Then run...

library("devtools")

devtools::install_github("nategarton13/CrimeDLM.RPackage" , build_vignettes = TRUE)

This may take a few minutes, as building the vignettes takes a bit of time.

For some users, the above code will fail with an RCMD.exe error. In that case, running 

devtools::install_github("nategarton13/CrimeDLM.RPackage" , build_vignettes = FALSE)

will install the package, but will not build the vignettes. After installing the package, you can build the vignettes manually by downloading or cloning this repo and knitting the .rmd vignette files. This is especially easy in RStudio.
