# For a fresh workstation. Set RNASEQ_LIBRARY to an independent analysis library.
# Current Cloud dependencies are already prepared; this is not a renv lockfile update.
lib <- Sys.getenv('RNASEQ_LIBRARY', file.path(getwd(),'analysis-library'))
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
cran <- c('BiocManager','ggplot2','cowplot')
need <- cran[!vapply(cran,requireNamespace,logical(1),quietly=TRUE)]
if(length(need)) install.packages(need,repos='https://cloud.r-project.org',lib=lib)
bio <- c('DESeq2','AnnotationDbi','GO.db','org.Hs.eg.db','org.Mm.eg.db')
need <- bio[!vapply(bio,requireNamespace,logical(1),quietly=TRUE)]
if(length(need)) BiocManager::install(need,ask=FALSE,update=FALSE,lib=lib)
stopifnot(all(vapply(c(cran,bio),requireNamespace,logical(1),quietly=TRUE)))
