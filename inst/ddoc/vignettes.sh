R CMD Sweave Tutorial_FullsibQTL.Rnw
%pdflatex Tutorial_FullsibQTL.tex
%pdflatex Tutorial_FullsibQTL.tex
latex Tutorial_FullsibQTL.tex
latex Tutorial_FullsibQTL.tex
dvips -Z Tutorial_FullsibQTL.dvi -o Tutorial_FullsibQTL.ps
ps2pdf Tutorial_FullsibQTL.ps Tutorial_FullsibQTL.pdf

rm Tutorial_FullsibQTL.dvi
rm Tutorial_FullsibQTL.ps
rm Tutorial_FullsibQTL.tex
rm Tutorial_FullsibQTL.log
rm Tutorial_FullsibQTL.aux
evince Tutorial_FullsibQTL.pdf