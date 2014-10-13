#!/bin/sh

rm top.aux top.bbl 
pdflatex top 
bibtex top 
perl -i -p -e 's|^(\\begin{thebibliography}{.*})$|$1\n\\normalsize|' top.bbl 
pdflatex top 
pdflatex top 
