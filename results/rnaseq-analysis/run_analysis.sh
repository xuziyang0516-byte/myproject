#!/bin/bash
set -euo pipefail
export RNASEQ_PROJECT="$(cd "$(dirname "$0")" && pwd)"
if [ -x /workspace/.cloud-setup/bin/Rscript ]; then export PATH=/workspace/.cloud-setup/bin:$PATH; fi
if [ -d "$RNASEQ_PROJECT/dependencies/usr/lib/R/site-library" ]; then
 export R_LIBS="$RNASEQ_PROJECT/dependencies/usr/lib/R/site-library:$RNASEQ_PROJECT/dependencies/usr/lib/R/library${R_LIBS:+:$R_LIBS}"
fi
export RENV_CONFIG_AUTOLOADER_ENABLED=false
export XDG_CACHE_HOME="$RNASEQ_PROJECT/fontcache"
mkdir -p "$RNASEQ_PROJECT/results/qc" "$XDG_CACHE_HOME"
python3 "$RNASEQ_PROJECT/scripts/extract_selection.py"
Rscript --vanilla "$RNASEQ_PROJECT/scripts/qc.R"
Rscript --vanilla "$RNASEQ_PROJECT/scripts/paired_de.R"
Rscript --vanilla "$RNASEQ_PROJECT/scripts/module_figures.R"
Rscript --vanilla -e 'invisible(lapply(c("DESeq2","AnnotationDbi","GO.db","org.Hs.eg.db","org.Mm.eg.db","ggplot2","cowplot"),requireNamespace,quietly=TRUE));capture.output(sessionInfo(),file=file.path(Sys.getenv("RNASEQ_PROJECT"),"results","sessionInfo.txt"))'
if [ "${1:-}" = "--full-ora" ]; then
 Rscript --vanilla "$RNASEQ_PROJECT/scripts/full_pathway_ora.R"
 Rscript --vanilla "$RNASEQ_PROJECT/scripts/plot_full_ora.R"
 Rscript --vanilla "$RNASEQ_PROJECT/scripts/paper_figures.R"
fi
