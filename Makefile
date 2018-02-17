createPackage:
	Rscript  --default-packages=methods,utils,graphics create_lib.R --no-save

build: createPackage
	R CMD build --resave-data working/ExomeDepth



