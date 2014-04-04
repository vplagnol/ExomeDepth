
both: run_code
	pdflatex vignette.tex

run_code: vignette.Rnw
	R CMD Sweave vignette.Rnw

clean :
	rm -f vignette.aux vignette.log vignette.pdf vignette.tex vignette.toc fig/*.pdf

