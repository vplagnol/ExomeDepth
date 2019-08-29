base <- 'working'
if (!file.exists(base)) dir.create(base)

system(paste('rm ', base, '/ExomeDepth/R/*', sep = ''))
package.skeleton(name="ExomeDepth",
                 code_files = c('R/class_definition.R', 
		 'R/optimize_reference_set.R', 
		 'R/tools.R',		       
		 'R/countBamInGranges.R', 
		 'R/plot_CNVs_method.R',
		 'R/annotate_extra.R'),
                 path= base,
                 force=TRUE)

file.copy(from = "doc/_pkgdown.yml", to = paste0(base, "/ExomeDepth/"), overwrite = TRUE)

for (folder in c('data', 'src', 'vignettes', 'inst', 'inst/doc')) {
  my.folder <- paste(base, '/ExomeDepth/', folder, sep = '')
  if (!file.exists(my.folder)) dir.create(path = my.folder)
}


methods.file <- list.files(path = paste(base, '/ExomeDepth/man', pattern = '*methods*', full.names = TRUE))
file.remove(methods.file)
print(list.files(path = paste(base, '/ExomeDepth/man', sep = ''), pattern = '*methods*', full.names = TRUE))

system(paste("cp src/*.cpp src/*.c src/*.h ", base, '/ExomeDepth/src', sep = ''))
file.copy(from = "README.md", to = paste0(base, "/ExomeDepth/"), overwrite = TRUE)
                        
file.copy(from = 'R/zzz.R', to= paste(base, '/ExomeDepth/R/zzz.R', sep = ''), overwrite = TRUE)


file.copy(from = 'NAMESPACE', to = paste(base, '/ExomeDepth/NAMESPACE', sep = ''), overwrite = TRUE)
file.copy(from = 'DESCRIPTION', to = paste(base, '/ExomeDepth/DESCRIPTION', sep= ''), overwrite = TRUE)

file.copy(from = "vignette/vignette.Rnw", to = paste(base, '/ExomeDepth/vignettes/ExomeDepth-vignette.Rnw', sep =''), overwrite = TRUE)

file.copy(from = 'data/ExomeCount.RData', to = paste(base, '/ExomeDepth/data/ExomeCount.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/exons.hg19.RData', to = paste(base, '/ExomeDepth/data/exons.hg19.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/genes.hg19.RData', to = paste(base, '/ExomeDepth/data/genes.hg19.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/exons.hg19.X.RData', to = paste(base, '/ExomeDepth/data/exons.hg19.X.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/Conrad.hg19.RData', to = paste(base, '/ExomeDepth/data/Conrad.hg19.RData', sep = ''), overwrite = TRUE)

clean.Rd.files <- file.remove(list.files(paste0(base, "/ExomeDepth/man"), 
	       full.names = TRUE)) ## now useful to let roxygen write the doc
file.remove(paste0(base, "/ExomeDepth/Read-and-delete-me"))

## now oxygenize
pkgbuild::compile_dll(paste0(base, "/ExomeDepth"))
devtools::document(paste0(base, "/ExomeDepth"))
