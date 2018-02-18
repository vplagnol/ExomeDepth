createPackage:
	Rscript  --default-packages=methods,utils,graphics create_lib.R --no-save

build: createPackage
	R CMD build --resave-data working/ExomeDepth

checks: build
	R CMD check --as-cran  ExomeDepth_1.1.11.tar.gz

html: build
	Rscript -e "setwd('working/ExomeDepth'); devtools::load_all(); pkgdown::build_site()"

