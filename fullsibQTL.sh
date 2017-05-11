rm fullsibQTL/src/fullsibQTL.o
rm fullsibQTL/src/fullsibQTL.so
#R CMD check fullsibQTL
#R CMD check --as-cran --timings fullsibQTL
#R CMD check --timings fullsibQTL
R CMD build fullsibQTL
R CMD INSTALL fullsibQTL_1.0-0.tar.gz
