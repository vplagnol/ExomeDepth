full: createPackage build

createPackage:
	Rscript  --default-packages=methods,utils,graphics util/create_lib.R --no-save

build: createPackage
	R CMD build --resave-data working/ExomeDepth

checks: 
	R CMD check --as-cran  ExomeDepth_1.1.12.tar.gz

make_vignette:
	cd vignette && Rscript -e 'knitr::knit("vignette.Rnw")' && pdflatex vignette.tex

html: build
	Rscript -e "setwd('working/ExomeDepth'); devtools::load_all(); pkgdown::build_site()"

