#Rbin <- "/cluster/project8/vyp/vincent/Software/R-3.2.0/bin/R"
#Rbin <- "/usr/bin/R"
Rbin <- "/share/apps/R-3.2.2/bin/R"
#Rbin <- " ../R-devel/bin/R"


base <- 'working'
system(paste('rm ', base, '/ExomeDepth/R/*', sep = ''))
package.skeleton(name="ExomeDepth",
                 code_files = c('R/class_definition.R', 'R/optimize_reference_set.R', 'R/tools.R', 'R/countBamInGranges.R', 'R/plot_CNVs_method.R'),
                 path= base,
                 force=TRUE)


for (folder in c('data', 'src', 'vignettes', 'inst', 'inst/doc')) {
  my.folder <- paste(base, '/ExomeDepth/', folder, sep = '')
  if (!file.exists(my.folder)) dir.create(path = my.folder)
}

file.remove (paste(base, '/ExomeDepth/man/CallCNVs-methods.Rd', sep = ''))
file.remove (paste(base, '/ExomeDepth/man/TestCNV-methods.Rd', sep = ''))
file.remove (paste(base, '/ExomeDepth/man/plot-methods.Rd', sep = ''))

methods.file <- list.files(path = paste(base, '/ExomeDepth/man', pattern = '*methods*', full.names = TRUE))
file.remove(methods.file)
print(list.files(path = paste(base, '/ExomeDepth/man', sep = ''), pattern = '*methods*', full.names = TRUE))

system(paste("cp src/*.cpp src/*.c src/*.h ", base, '/ExomeDepth/src', sep = ''))
                        
file.copy(from = 'R/zzz.R', to= paste(base, '/ExomeDepth/R/zzz.R', sep = ''), overwrite = TRUE)


file.copy(from = 'doc/NAMESPACE', to = paste(base, '/ExomeDepth/NAMESPACE', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/DESCRIPTION', to = paste(base, '/ExomeDepth/DESCRIPTION', sep= ''), overwrite = TRUE)
file.copy(from = 'data/ExomeCount.RData', to = paste(base, '/ExomeDepth/data/ExomeCount.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/exons.hg19.RData', to = paste(base, '/ExomeDepth/data/exons.hg19.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/genes.hg19.RData', to = paste(base, '/ExomeDepth/data/genes.hg19.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/exons.hg19.X.RData', to = paste(base, '/ExomeDepth/data/exons.hg19.X.RData', sep = ''), overwrite = TRUE)
file.copy(from = 'data/Conrad.hg19.RData', to = paste(base, '/ExomeDepth/data/Conrad.hg19.RData', sep = ''), overwrite = TRUE)

file.copy(from = 'doc/exons.hg19.Rd', to= paste(base, '/ExomeDepth/man/exons.hg19.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/exons.hg19.X.Rd', to = paste(base, '/ExomeDepth/man/exons.hg19.X.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/Conrad.hg19.common.CNVs.Rd', to = paste(base, '/ExomeDepth/man/Conrad.hg19.common.CNVs.Rd', sep = ''), overwrite = TRUE)

file.copy(from = 'doc/ExomeCount.Rd', to= paste(base, '/ExomeDepth/man/ExomeCount.Rd', sep = ''), overwrite = TRUE)

file.copy(from = 'doc/ExomeDepth-class.Rd', to= paste(base, '/ExomeDepth/man/ExomeDepth-class.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/countBamInGRanges.exomeDepth.Rd', to= paste(base, '/ExomeDepth/man/countBamInGRanges.exomeDepth.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/getBamCounts.Rd', to= paste(base, '/ExomeDepth/man/getBamCounts.Rd', sep = ''),overwrite = TRUE)

#file.copy(from = 'doc/AddAnnotations.Rd', to= '/ugi/home/shared/vincent/libraries/R/working/ExomeDepth/man/AddAnnotations.Rd', overwrite = TRUE)
file.copy(from = 'doc/CallCNVs.Rd', to= paste(base, '/ExomeDepth/man/CallCNVs.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/TestCNV.Rd', to= paste(base, '/ExomeDepth/man/TestCNV.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/AnnotateExtra-methods.Rd', to= paste(base, '/ExomeDepth/man/AnnotateExtra-methods.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/AnnotateExtra.Rd', to= paste(base, '/ExomeDepth/man/AnnotateExtra.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/count.everted.reads.Rd', to= paste(base, '/ExomeDepth/man/count.everted.reads.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/countBam.everted.Rd', to= paste(base, '/ExomeDepth/man/countBam.everted.Rd', sep = ''), overwrite = TRUE)

file.copy(from = 'doc/genes.hg19.Rd', to= paste(base, '/ExomeDepth/man/genes.hg19.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/exons.hg19.Rd', to= paste(base, '/ExomeDepth/man/exons.hg19.Rd', sep = ''), overwrite = TRUE)

file.copy(from = 'doc/plot.ExomeDepth.Rd', to= paste(base, '/ExomeDepth/man/plot.ExomeDepth.Rd', sep = ''), overwrite = TRUE)

file.copy(from = 'doc/somatic.CNV.call.Rd', to= paste(base, '/ExomeDepth/man/somatic.CNV.call.Rd', sep = ''), overwrite = TRUE)


file.copy(from = 'doc/select.reference.set.Rd', to= paste(base, '/ExomeDepth/man/select.reference.set.Rd',sep = ''), overwrite = TRUE)
file.copy(from = 'doc/qbetabinom.Rd', to= paste(base, '/ExomeDepth/man/qbetabinom.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/get.power.betabinom.Rd', to= paste(base, '/ExomeDepth/man/get.power.betabinom.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/qbetabinom.ab.Rd', to= paste(base, '/ExomeDepth/man/qbetabinom.ab.Rd', sep = ''), overwrite = TRUE)
file.copy(from = 'doc/viterbi.hmm.Rd', to= paste(base, '/ExomeDepth/man/viterbi.hmm.Rd', sep = ''), overwrite = TRUE)

file.copy(from = 'doc/ExomeDepth-package.Rd', to= paste(base, '/ExomeDepth/man/ExomeDepth-package.Rd', sep = ''), overwrite = TRUE)
file.copy(from = "vignette/vignette.Rnw", to = paste(base, '/ExomeDepth/vignettes/ExomeDepth-vignette.Rnw', sep =''), overwrite = TRUE)



#system(paste(Rbin, " CMD build ", base, "/ExomeDepth", sep = ''))
system(paste(Rbin, " CMD build --resave-data ", base, "/ExomeDepth", sep = ''))
system(paste(Rbin, " CMD INSTALL  ExomeDepth_1.1.10.tar.gz"))
system(paste('cp /home/ucbtvyp/vyp/vincent/libraries/R/installed/ExomeDepth/doc/ExomeDepth-vignette.pdf .'))

##system("/cluster/project8/vyp/vincent/Software/R-devel/bin/R CMD check --as-cran ExomeDepth_1.1.10.tar.gz")









