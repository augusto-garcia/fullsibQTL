# fullsibQTL

[![Build Status](https://travis-ci.org/augusto-garcia/fullsibQTL.svg?branch=master)](https://travis-ci.org/augusto-garcia/fullsibQTL) [![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)

<!-- [![Build Status](https://travis-ci.org/augusto-garcia/fullsibQTL.svg?branch=master)](https://travis-ci.org/augusto-garcia/fullsibQTL) -->

# fullsibQTL

`fullsibQTL` is a software for for QTL mapping in outbred (outcrossing) species, considering a full-sib cross (or F1 population) as a mapping population. It is assumed a previous genetic map obtained with `onemap` package. The method is based on maximum likelihood approach using mixture models and EM algorithm. The probabilities for QTL genotypes are inferred using multipoint approach based on Hidden Markov Models. We implemented functions to perform interval mapping (IM) and composite interval mapping (CIM). In a first step one need to scan the genome for QTLs and after locating them, one can characterize these loci, *i.e.*, it is possible to estimate their effects, segregation pattern and the linkage phase between QTL and markers.

It has been available on CRAN (http://cran.r-project.org/web/packages/fullsibQTL/index.html). Its last version was updated on YYYY-MM-DD. CRAN has `fullsibQTL`'s stable version, which is recommended for most users.

This github page has its version under development. New functions will be added (test phase) and, once it is done, we will synchronizethe repositories and add it to CRAN.

`fullsibQTL` comprises a set of functions that allows users to QTL mapping. Some functions are used internally by the package, and should not be used directly.

# How to install

## From CRAN (stable version)

It is easy, just type (within R):

```R
install.packages( "fullsibQTL" )
```

The `onemap` package is a dependency, if you are having trouble installing it please check its git page (https://github.com/augusto-garcia/onemap)

## From github (version under development)

Within R, you need to install and load the package `devtools` and install `fullsibQTL` from this very repository:

```R
install.packages( "devtools" )
library( "devtools" )

install_github( "augusto-garcia/fullsibQTL" )
```

This will allow you to automatically build and install packages from github. If you use Windows and the above script did not work, try to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. If you are facing problems with Rtools installation, try to do it by selecting *Run as Admnistrator* option with right mouse button. On a Mac, you will need Xcode (available on the App Store). On Linux, you may need to install `r-cran-tkrplot`.

Then, to install `fullsibQTL` from github (this very repo):

```R
install_github( "augusto-garcia/fullsibQTL" )
```

<div>
<iframe src='https://rramadeu.github.io/fullsibQTL/vignettes_highres/graphics/cim_im.html', height="600", width = "650"></iframe>
</div>

# Tutorials

You can read `fullsibQTL` tutorial going to the vignettes of the installed package, or clicking below.

[fullsibQTL Tutorial](https://rramadeu.github.io/fullsibQTL/vignettes_highres/fullsibQTL_Tutorial.html)

# Additional material
[Graphical options with plot_fullsibQTL](https://rramadeu.github.io/fullsibQTL/vignettes_highres/Graphical_options_with_plot_fullsibQTL.html)

[Graphical options with plot](https://rramadeu.github.io/fullsibQTL/vignettes_highres/Graphical_options_with_plot.html)

[QTL mapping with partially informative markers](https://rramadeu.github.io/fullsibQTL/vignettes_highres/QTL_mapping_with_partially_informative_markers.html)
