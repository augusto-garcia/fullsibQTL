# fullsibQTL

[![Build Status](https://travis-ci.org/augusto-garcia/fullsibQTL.svg?branch=master)](https://travis-ci.org/augusto-garcia/fullsibQTL) [![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)

<!-- [![Build Status](https://travis-ci.org/augusto-garcia/fullsibQTL.svg?branch=master)](https://travis-ci.org/augusto-garcia/fullsibQTL) -->

# fullsibQTL

`fullsibQTL` is a software for QTL mapping in outbred (outcrossing) species, considering a full-sib cross (or F1 population) as a mapping population, using the model develop by Gazaffi et al 2014 (https://link.springer.com/article/10.1007/s11295-013-0664-2). It assumes that a previous genetic map was obtained with R package `onemap`.

The method is based on a maximum likelihood approach, using mixture models and the EM algorithm. The probabilities of QTL genotypes are inferred using a multipoint approach based on Hidden Markov Models. We implemented functions to perform interval mapping (IM) and composite interval mapping (CIM) for F1 populations. In a first step, one needs to scan the genome for QTLs and, after locating them, characterize these loci, estimating their effects, segregation pattern and the linkage phase between markers and QTL.

The package on its stable version is available on CRAN (http://cran.r-project.org/web/packages/fullsibQTL/index.html). This version is recommended for most users. The version on github has new experimental features and new developments, and should be useful for some audiences. Once they are tested, we will synchronize the repositories.

`fullsibQTL` comprises a set of functions that allows users to perform QTL mapping. Some functions are used internally by the package, and should not be used directly.

`fullsibQTL` also has interactive graphical features to help in the analysis, that will facilitate most of the tasks related to a QTL analysis:

![http://i.imgur.com/Macz6ph.gif](http://i.imgur.com/Macz6ph.gif)

# How to install

## From CRAN (stable version)

It is easy, just type (within R):

```R
install.packages( "fullsibQTL" )
```

The `onemap` package is a dependency, so if you are having trouble installing `fullsibQTL` due to this dependency, please check `onempa` github page (https://github.com/augusto-garcia/onemap) for information about how to install it.

## From github (version under development)

Within R, you need to install and load the package `devtools` and install `fullsibQTL` from this very repository:

```R
install.packages( "devtools" )
library( "devtools" )
install_github( "augusto-garcia/fullsibQTL" )
```

This will allow you to automatically build and install packages from github. If you use MS Windows and the above script did not work, try to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. If you are facing problems with Rtools installation, try to do it by selecting *Run as Admnistrator* option with right mouse button. On a Mac, you will need _Xcode_ (available on the App Store). On Linux, you may need to install `r-cran-tkrplot` in the command line, since this is a `onemap` dependency that sometimes cause some problems.


# Tutorials

You can read `fullsibQTL` tutorial going to the vignettes of the installed package, or by clicking below.

Comments and suggestions are welcome!

[fullsibQTL Tutorial](https://rramadeu.github.io/fullsibQTL/vignettes_highres/fullsibQTL_Tutorial.html)

# Additional material
[Graphical options with plot_fullsibQTL](https://rramadeu.github.io/fullsibQTL/vignettes_highres/Graphical_options_with_plot_fullsibQTL.html)

[Graphical options with plot](https://rramadeu.github.io/fullsibQTL/vignettes_highres/Graphical_options_with_plot.html)

[QTL mapping with partially informative markers](https://rramadeu.github.io/fullsibQTL/vignettes_highres/QTL_mapping_with_partially_informative_markers.html)
