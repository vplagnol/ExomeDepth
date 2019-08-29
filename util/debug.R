source("scripts/ExoDepth/R/class_definition.R")
source("scripts/ExoDepth/R/optimize_reference_set.R")
source('scripts/ExoDepth/R/tools.R')
dyn.load('/ugi/home/shared/vincent/libraries/R/installed/ExomeDepth/libs/ExomeDepth.so' )

load('debug.RData')


print(.Object@phi)
print(.Object@expected[1:100])
print(as.integer(.Object@reference + .Object@test)[1:100])
print(as.integer(.Object@test)[1:100])

.Object@likelihood <- .Call("get_loglike_matrix",
                              phi = .Object@phi,
                              expected = .Object@expected[1:10],
                              total = as.integer(.Object@reference + .Object@test)[1:10],
                              observed = as.integer(.Object@test)[1:10])
