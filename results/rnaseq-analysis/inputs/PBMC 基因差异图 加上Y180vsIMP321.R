# =============================================================================
# PBMC: Sheet5-selected integrated GO/KEGG GSEA
#
# WHY GSEA HERE?
# --------------
# Y180 vs IMP321 has too few threshold-defined DEGs for ORA.
# GSEA avoids the Up/Down DEG cutoff and uses the complete ranked gene list.
#
# Ranking:
#   DESeq2 Wald statistic ("stat") from each complete all_genes.tsv table.
#
# Comparisons:
#   IMP vs Ctrl
#   Y180 vs Ctrl
#   Y180 vs IMP321
#
# Pathway selector:
#   Sheet5 = biologically relevant pathways to display
#
# IMPORTANT:
# 1. The Excel sheets select which pathway IDs are displayed only.
# 2. GSEA is performed against the COMPLETE GO/KEGG collections first.
# 3. BH adjustment is therefore performed before pathway selection.
# 4. The Excel GeneName column is not used as the GSEA gene list.
# 5. Plot aesthetics requested for this figure:
#      color = GSEA BH-adjusted P value
#      size  = number of threshold-defined UP DEGs overlapping the pathway
# 6. NES is retained in the exported result table for direction interpretation:
#      NES > 0: numerator group is enriched
#      NES < 0: denominator group is enriched
# =============================================================================


# =============================================================================
# 1. User settings
# =============================================================================

dataset_name <- "PBMC"

all_gene_files <- c(
  "IMP vs Ctrl" = paste0(
    "D:/个人记录/科研/Y180/PBMC Paired_GO_KEGG_results/",
    "IMP_vs_Ctrl/IMP_vs_Ctrl_paired_all_genes.tsv"
  ),
  "Y180 vs Ctrl" = paste0(
    "D:/个人记录/科研/Y180/PBMC Paired_GO_KEGG_results/",
    "Y180_vs_Ctrl/Y180_vs_Ctrl_paired_all_genes.tsv"
  ),
  "Y180 vs IMP321" = paste0(
    "D:/个人记录/科研/Y180/PBMC Paired_GO_KEGG_results/",
    "Y180_vs_IMP321/Y180_vs_IMP321_paired_all_genes.tsv"
  )
)

selected_pathway_excel <- paste0(
  "D:/个人记录/科研/Y180/PBMC Paired_GO_KEGG_results/",
  "final 10 selected.xlsx"
)

selected_sheet_direction <- c(
  Sheet5 = "Up"
)

comparison_order <- c(
  "IMP vs Ctrl",
  "Y180 vs Ctrl",
  "Y180 vs IMP321"
)

# DEG Count uses the same threshold definition as the prior ORA analysis.
# Because Sheet5 is the UP selector, the plotted count is the number of
# UP DEGs from each comparison that overlap the pathway.
deg_padj_cutoff <- 0.05
deg_abs_shrunken_log2fc_cutoff <- 1

minimum_gene_set_size <- 10L
maximum_gene_set_size <- 500L

# Keep all GSEA results first; selected pathways are filtered afterward.
gsea_pvalue_cutoff <- 1

# Setting eps = 0 requests more accurate fgsea P-value estimation.
gsea_eps <- 0

# Reproducibility for fgsea-based calculations.
gsea_seed <- 20260825L

# FALSE gives labels such as "GO:0006935, chemotaxis".
use_excel_display_label <- FALSE

output_dir <- file.path(
  dirname(selected_pathway_excel),
  "selected_pathways_GSEA_Sheet5_AdjustedP_DEGCount_KEGG_integrated"
)

complete_gsea_dir <- file.path(
  output_dir,
  "01_complete_GSEA"
)

plot_dir <- file.path(
  output_dir,
  "02_selected_pathway_plots"
)

dir.create(
  complete_gsea_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 2. Packages
# =============================================================================

cran_packages <- c(
  "readxl",
  "ggplot2",
  "msigdbr",
  "scales"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_core_packages <- c(
  "clusterProfiler",
  "AnnotationDbi"
)

missing_bioc_core <- bioc_core_packages[
  !vapply(
    bioc_core_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc_core) > 0) {
  BiocManager::install(
    missing_bioc_core,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(msigdbr)
  library(scales)
  library(clusterProfiler)
  library(AnnotationDbi)
})


# =============================================================================
# 3. General helpers
# =============================================================================

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("　", " ", x, fixed = TRUE)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "NULL", "-")] <- NA_character_
  x
}


write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}


remove_blank_columns <- function(x) {
  keep <- vapply(
    x,
    function(column) {
      any(
        !is.na(
          clean_text(column)
        )
      )
    },
    FUN.VALUE = logical(1)
  )
  
  x[
    ,
    keep,
    drop = FALSE
  ]
}


normalize_database <- function(
    database,
    source_file,
    panel,
    id
) {
  database <- clean_text(database)
  source_file <- clean_text(source_file)
  panel <- clean_text(panel)
  id <- clean_text(id)
  
  context <- toupper(
    paste(
      database,
      source_file,
      panel,
      sep = " | "
    )
  )
  
  context <- gsub(
    "[[:space:].-]+",
    "_",
    context
  )
  
  answer <- rep(
    NA_character_,
    length(id)
  )
  
  answer[
    grepl(
      "^(HSA|MMU|KO)[0-9]+",
      toupper(id)
    )
  ] <- "KEGG"
  
  answer[
    is.na(answer) &
      grepl(
        "KEGG",
        context
      )
  ] <- "KEGG"
  
  answer[
    is.na(answer) &
      grepl(
        "GO_?BP|BIOLOGICAL_PROCESS",
        context
      )
  ] <- "GO_BP"
  
  answer[
    is.na(answer) &
      grepl(
        "GO_?MF|MOLECULAR_FUNCTION",
        context
      )
  ] <- "GO_MF"
  
  answer[
    is.na(answer) &
      grepl(
        "GO_?CC|CELLULAR_COMPONENT",
        context
      )
  ] <- "GO_CC"
  
  answer
}


save_plot <- function(
    plot_object,
    stem,
    width,
    height
) {
  ggplot2::ggsave(
    filename = file.path(
      plot_dir,
      paste0(
        stem,
        ".png"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )
  
  tryCatch(
    ggplot2::ggsave(
      filename = file.path(
        plot_dir,
        paste0(
          stem,
          ".pdf"
        )
      ),
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    ),
    error = function(e) {
      warning(
        "PNG saved, but PDF export failed: ",
        conditionMessage(e)
      )
    }
  )
}


# =============================================================================
# 4. Detect organism from the complete DE tables
# =============================================================================

read_ids_for_species <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(
      "Cannot find complete DE table:\n",
      file_path
    )
  }
  
  x <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  
  colnames(x) <- make.unique(
    clean_text(
      colnames(x)
    ),
    sep = "_"
  )
  
  id_column <- if (
    "EnsemblID_clean" %in%
    colnames(x)
  ) {
    "EnsemblID_clean"
  } else if (
    "Geneid" %in%
    colnames(x)
  ) {
    "Geneid"
  } else {
    stop(
      "The complete DE table has no EnsemblID_clean/Geneid column:\n",
      file_path
    )
  }
  
  ids <- sub(
    "\\.[0-9]+$",
    "",
    clean_text(
      x[
        ,
        id_column,
        drop = TRUE
      ]
    )
  )
  
  unique(
    ids[
      !is.na(ids) &
        grepl(
          "^(ENSG|ENSMUSG)[0-9]+$",
          ids
        )
    ]
  )
}


species_ids <- unlist(
  lapply(
    all_gene_files,
    read_ids_for_species
  ),
  use.names = FALSE
)

if (length(species_ids) == 0) {
  stop(
    "No human or mouse Ensembl IDs were recovered."
  )
}

human_fraction <- mean(
  grepl(
    "^ENSG[0-9]+$",
    species_ids
  )
)

mouse_fraction <- mean(
  grepl(
    "^ENSMUSG[0-9]+$",
    species_ids
  )
)

if (mouse_fraction > 0.95) {
  analysis_species <- "mouse"
  orgdb_package <- "org.Mm.eg.db"
  kegg_prefix <- "mmu"
  valid_ensembl_pattern <- "^ENSMUSG[0-9]+$"
  
} else if (human_fraction > 0.95) {
  analysis_species <- "human"
  orgdb_package <- "org.Hs.eg.db"
  kegg_prefix <- "hsa"
  valid_ensembl_pattern <- "^ENSG[0-9]+$"
  
} else {
  stop(
    paste0(
      "Mixed or unsupported Ensembl IDs.\n",
      "Human fraction: ",
      round(human_fraction, 3),
      "\nMouse fraction: ",
      round(mouse_fraction, 3)
    )
  )
}

if (
  !requireNamespace(
    orgdb_package,
    quietly = TRUE
  )
) {
  BiocManager::install(
    orgdb_package,
    ask = FALSE,
    update = FALSE
  )
}

OrgDb <- get(
  orgdb_package,
  envir = asNamespace(
    orgdb_package
  )
)

message(
  dataset_name,
  " species detected: ",
  analysis_species,
  "; OrgDb = ",
  orgdb_package,
  "; KEGG prefix = ",
  kegg_prefix
)


# =============================================================================
# 5. Species-aware annotation helpers
# =============================================================================

normalize_kegg_pathway_id <- function(x) {
  x <- clean_text(x)
  
  x <- sub(
    "^path:",
    "",
    x,
    ignore.case = TRUE
  )
  
  x <- sub(
    "^(hsa|mmu|ko)",
    "",
    x,
    ignore.case = TRUE
  )
  
  answer <- ifelse(
    !is.na(x) &
      grepl(
        "^[0-9]{5}$",
        x
      ),
    paste0(
      kegg_prefix,
      x
    ),
    x
  )
  
  clean_text(answer)
}


infer_go_database_from_id <- function(go_ids) {
  go_ids <- clean_text(go_ids)
  
  answer <- rep(
    NA_character_,
    length(go_ids)
  )
  
  go_positions <- which(
    !is.na(go_ids) &
      grepl(
        "^GO:[0-9]{7}$",
        go_ids
      )
  )
  
  if (length(go_positions) == 0) {
    return(answer)
  }
  
  candidate_ids <- unique(
    go_ids[
      go_positions
    ]
  )
  
  valid_ids <- intersect(
    candidate_ids,
    AnnotationDbi::keys(
      OrgDb,
      keytype = "GOALL"
    )
  )
  
  if (length(valid_ids) == 0) {
    return(answer)
  }
  
  ontology_table <- suppressMessages(
    AnnotationDbi::select(
      OrgDb,
      keys = valid_ids,
      keytype = "GOALL",
      columns = "ONTOLOGYALL"
    )
  )
  
  ontology_table <- ontology_table[
    !is.na(
      ontology_table$GOALL
    ) &
      !is.na(
        ontology_table$ONTOLOGYALL
      ),
    c(
      "GOALL",
      "ONTOLOGYALL"
    ),
    drop = FALSE
  ]
  
  ontology_table <- unique(
    ontology_table
  )
  
  ontology_by_id <- tapply(
    ontology_table$ONTOLOGYALL,
    ontology_table$GOALL,
    function(x) {
      unique(x)[1]
    }
  )
  
  matched_ontology <- unname(
    ontology_by_id[
      go_ids[
        go_positions
      ]
    ]
  )
  
  answer[
    go_positions
  ] <- ifelse(
    is.na(matched_ontology),
    NA_character_,
    paste0(
      "GO_",
      matched_ontology
    )
  )
  
  answer
}


# =============================================================================
# 6. Read complete all-gene tables and build ranked vectors
# =============================================================================

read_rank_table <- function(
    file_path,
    comparison_name
) {
  x <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  
  colnames(x) <- make.unique(
    clean_text(
      colnames(x)
    ),
    sep = "_"
  )
  
  ensembl_column <- if (
    "EnsemblID_clean" %in%
    colnames(x)
  ) {
    "EnsemblID_clean"
  } else if (
    "Geneid" %in%
    colnames(x)
  ) {
    "Geneid"
  } else {
    stop(
      comparison_name,
      " lacks EnsemblID_clean/Geneid."
    )
  }
  
  x$EnsemblIDForGSEA <- sub(
    "\\.[0-9]+$",
    "",
    clean_text(
      x[
        ,
        ensembl_column,
        drop = TRUE
      ]
    )
  )
  
  valid_id <- (
    !is.na(
      x$EnsemblIDForGSEA
    ) &
      grepl(
        valid_ensembl_pattern,
        x$EnsemblIDForGSEA
      )
  )
  
  x <- x[
    valid_id,
    ,
    drop = FALSE
  ]
  
  # Reconstruct one uniform Direction definition when necessary.
  if (
    !"Direction" %in%
    colnames(x)
  ) {
    required_direction_columns <- c(
      "padj",
      "log2FoldChange_shrunken"
    )
    
    missing_direction_columns <- setdiff(
      required_direction_columns,
      colnames(x)
    )
    
    if (
      length(
        missing_direction_columns
      ) > 0
    ) {
      stop(
        comparison_name,
        " lacks Direction and the columns needed to reconstruct it: ",
        paste(
          missing_direction_columns,
          collapse = ", "
        )
      )
    }
    
    padj_for_direction <- suppressWarnings(
      as.numeric(
        x[
          ,
          "padj",
          drop = TRUE
        ]
      )
    )
    
    shrunken_lfc_for_direction <- suppressWarnings(
      as.numeric(
        x[
          ,
          "log2FoldChange_shrunken",
          drop = TRUE
        ]
      )
    )
    
    x$Direction <- "Not_significant"
    
    valid_deg <- (
      !is.na(
        padj_for_direction
      ) &
        padj_for_direction <
        deg_padj_cutoff &
        !is.na(
          shrunken_lfc_for_direction
        )
    )
    
    x$Direction[
      valid_deg &
        shrunken_lfc_for_direction >=
        deg_abs_shrunken_log2fc_cutoff
    ] <- "Up"
    
    x$Direction[
      valid_deg &
        shrunken_lfc_for_direction <=
        -deg_abs_shrunken_log2fc_cutoff
    ] <- "Down"
  }
  
  direction_upper <- toupper(
    trimws(
      as.character(
        x$Direction
      )
    )
  )
  
  x$DirectionForCount <- "Not_significant"
  x$DirectionForCount[
    direction_upper == "UP"
  ] <- "Up"
  x$DirectionForCount[
    direction_upper == "DOWN"
  ] <- "Down"
  
  if ("stat" %in% colnames(x)) {
    rank_score <- suppressWarnings(
      as.numeric(
        x[
          ,
          "stat",
          drop = TRUE
        ]
      )
    )
    
    rank_method <- "DESeq2 Wald statistic"
    
  } else if (
    all(
      c(
        "log2FoldChange",
        "pvalue"
      ) %in%
      colnames(x)
    )
  ) {
    lfc <- suppressWarnings(
      as.numeric(
        x[
          ,
          "log2FoldChange",
          drop = TRUE
        ]
      )
    )
    
    p_value <- suppressWarnings(
      as.numeric(
        x[
          ,
          "pvalue",
          drop = TRUE
        ]
      )
    )
    
    rank_score <- sign(lfc) *
      -log10(
        pmax(
          p_value,
          .Machine$double.xmin
        )
      )
    
    rank_method <- "signed -log10(P value) fallback"
    
    warning(
      comparison_name,
      ": stat column missing; using signed -log10(P value)."
    )
    
  } else {
    stop(
      comparison_name,
      " has neither stat nor log2FoldChange+pvalue columns."
    )
  }
  
  x$RankScore <- rank_score
  
  x <- x[
    is.finite(
      x$RankScore
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(x) == 0) {
    stop(
      comparison_name,
      ": no finite ranking values."
    )
  }
  
  # Retain one row per Ensembl ID using the strongest absolute ranking value.
  x <- x[
    order(
      -abs(
        x$RankScore
      )
    ),
    ,
    drop = FALSE
  ]
  
  x <- x[
    !duplicated(
      x$EnsemblIDForGSEA
    ),
    ,
    drop = FALSE
  ]
  
  entrez_map <- suppressMessages(
    AnnotationDbi::mapIds(
      OrgDb,
      keys = unique(
        x$EnsemblIDForGSEA
      ),
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  
  x$ENTREZIDForGSEA <- clean_text(
    unname(
      entrez_map[
        x$EnsemblIDForGSEA
      ]
    )
  )
  
  x$Comparison <- comparison_name
  x$RankMethod <- rank_method
  
  x
}


all_rank_tables <- lapply(
  names(all_gene_files),
  function(comparison_name) {
    read_rank_table(
      all_gene_files[[comparison_name]],
      comparison_name
    )
  }
)

names(all_rank_tables) <- names(all_gene_files)


make_named_rank <- function(
    values,
    ids
) {
  keep <- (
    is.finite(values) &
      !is.na(ids) &
      ids != ""
  )
  
  values <- values[
    keep
  ]
  
  ids <- ids[
    keep
  ]
  
  ranking_table <- data.frame(
    ID = ids,
    Score = values,
    stringsAsFactors = FALSE
  )
  
  ranking_table <- ranking_table[
    order(
      -abs(
        ranking_table$Score
      )
    ),
    ,
    drop = FALSE
  ]
  
  ranking_table <- ranking_table[
    !duplicated(
      ranking_table$ID
    ),
    ,
    drop = FALSE
  ]
  
  ranking <- ranking_table$Score
  names(ranking) <- ranking_table$ID
  
  ranking <- sort(
    ranking,
    decreasing = TRUE
  )
  
  ranking
}


go_rank_lists <- lapply(
  all_rank_tables,
  function(x) {
    make_named_rank(
      values = x$RankScore,
      ids = x$EnsemblIDForGSEA
    )
  }
)


# KEGG MSigDB gene sets below are represented with target-species
# Ensembl IDs, so the same complete Ensembl ranking is used for KEGG GSEA.
kegg_rank_lists <- go_rank_lists


rank_summary <- do.call(
  rbind,
  lapply(
    names(all_rank_tables),
    function(comparison_name) {
      current_table <- all_rank_tables[[comparison_name]]
      
      data.frame(
        Comparison = comparison_name,
        RankMethod = unique(
          current_table$RankMethod
        )[1],
        RankedEnsemblGenes = length(
          go_rank_lists[[comparison_name]]
        ),
        RankedKEGGEnsemblGenes = length(
          kegg_rank_lists[[comparison_name]]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write_tsv(
  rank_summary,
  file.path(
    output_dir,
    "00_ranked_gene_list_summary.tsv"
  )
)

print(
  rank_summary,
  row.names = FALSE
)


# =============================================================================
# 7. Build a more complete species-aware KEGG TERM2GENE from MSigDB
#
# The prior OrgDb PATH mapping omitted several selected KEGG pathways.
# Here the human MSigDB KEGG_LEGACY and KEGG_MEDICUS collections are mapped
# to the detected target species by msigdbr. GO analysis remains unchanged.
#
# This section only supplies the KEGG gene-set database. The final figure still
# combines GO and KEGG rows in one Sheet5-selected plot.
# =============================================================================

normalize_pathway_name <- function(x) {
  x <- toupper(
    clean_text(x)
  )
  
  x <- gsub(
    "KEGG|MEDICUS|LEGACY|REFERENCE",
    "",
    x
  )
  
  x <- gsub(
    "SIGNALLING|SIGNALING|PATHWAY",
    "",
    x
  )
  
  x <- gsub(
    "KAPPA",
    "K",
    x
  )
  
  x <- gsub(
    "[^A-Z0-9]",
    "",
    x
  )
  
  x
}


extract_prefixed_kegg_code <- function(x) {
  x <- toupper(
    ifelse(
      is.na(x),
      "",
      as.character(x)
    )
  )
  
  match_position <- regexpr(
    "(HSA|MMU|KO)[0-9]{5}",
    x,
    perl = TRUE
  )
  
  match_length <- attr(
    match_position,
    "match.length"
  )
  
  answer <- rep(
    NA_character_,
    length(x)
  )
  
  matched <- (
    match_position > 0 &
      match_length > 0
  )
  
  if (
    any(
      matched
    )
  ) {
    matched_tokens <- substr(
      x[
        matched
      ],
      match_position[
        matched
      ],
      match_position[
        matched
      ] +
        match_length[
          matched
        ] -
        1
    )
    
    answer[
      matched
    ] <- sub(
      "^(HSA|MMU|KO)",
      "",
      matched_tokens
    )
  }
  
  answer
}


msig_target_species <- if (
  analysis_species ==
  "mouse"
) {
  "Mus musculus"
} else {
  "Homo sapiens"
}

msig_collection_table <- as.data.frame(
  msigdbr::msigdbr_collections(
    db_species = "HS"
  ),
  stringsAsFactors = FALSE
)

kegg_subcollections <- unique(
  msig_collection_table$gs_subcollection[
    msig_collection_table$gs_collection ==
      "C2" &
      grepl(
        "KEGG",
        msig_collection_table$gs_subcollection,
        ignore.case = TRUE
      )
  ]
)

if (
  length(
    kegg_subcollections
  ) == 0
) {
  stop(
    "No KEGG subcollection was found in the installed msigdbr database."
  )
}

msig_kegg_list <- lapply(
  kegg_subcollections,
  function(subcollection_name) {
    message(
      "Loading ",
      subcollection_name,
      " from human MSigDB and mapping to ",
      msig_target_species,
      "..."
    )
    
    current_result <- tryCatch(
      msigdbr::msigdbr(
        db_species = "HS",
        species = msig_target_species,
        collection = "C2",
        subcollection = subcollection_name
      ),
      error = function(e) {
        warning(
          subcollection_name,
          " could not be loaded: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (
      is.null(
        current_result
      )
    ) {
      return(
        NULL
      )
    }
    
    current_result <- as.data.frame(
      current_result,
      stringsAsFactors = FALSE
    )
    
    current_result$RequestedSubcollection <- subcollection_name
    
    current_result
  }
)

msig_kegg_list <- Filter(
  Negate(
    is.null
  ),
  msig_kegg_list
)

if (
  length(
    msig_kegg_list
  ) == 0
) {
  stop(
    "No MSigDB KEGG collection could be loaded."
  )
}

msig_kegg <- do.call(
  rbind,
  msig_kegg_list
)

required_msig_columns <- c(
  "gs_name",
  "ensembl_gene"
)

missing_msig_columns <- setdiff(
  required_msig_columns,
  colnames(
    msig_kegg
  )
)

if (
  length(
    missing_msig_columns
  ) > 0
) {
  stop(
    "msigdbr output lacks columns: ",
    paste(
      missing_msig_columns,
      collapse = ", "
    )
  )
}

kegg_source_text <- paste(
  if (
    "gs_exact_source" %in%
    colnames(
      msig_kegg
    )
  ) {
    msig_kegg$gs_exact_source
  } else {
    ""
  },
  if (
    "gs_url" %in%
    colnames(
      msig_kegg
    )
  ) {
    msig_kegg$gs_url
  } else {
    ""
  },
  msig_kegg$gs_name
)

msig_kegg$KeggCode <- extract_prefixed_kegg_code(
  kegg_source_text
)

# When both KEGG_MEDICUS and KEGG_LEGACY contain the same KEGG code,
# prefer the MEDICUS definition to avoid merging two versions of one pathway.
msig_kegg$SourcePriority <- ifelse(
  grepl(
    "MEDICUS",
    msig_kegg$RequestedSubcollection,
    ignore.case = TRUE
  ),
  1L,
  ifelse(
    grepl(
      "LEGACY",
      msig_kegg$RequestedSubcollection,
      ignore.case = TRUE
    ),
    2L,
    3L
  )
)

coded_positions <- which(
  !is.na(
    msig_kegg$KeggCode
  )
)

if (
  length(
    coded_positions
  ) > 0
) {
  minimum_priority_by_code <- tapply(
    msig_kegg$SourcePriority[
      coded_positions
    ],
    msig_kegg$KeggCode[
      coded_positions
    ],
    min
  )
  
  keep_coded_row <- (
    is.na(
      msig_kegg$KeggCode
    ) |
      msig_kegg$SourcePriority ==
      unname(
        minimum_priority_by_code[
          msig_kegg$KeggCode
        ]
      )
  )
  
  keep_coded_row[
    is.na(
      keep_coded_row
    )
  ] <- TRUE
  
  msig_kegg <- msig_kegg[
    keep_coded_row,
    ,
    drop = FALSE
  ]
}

msig_kegg$NormalizedGsName <- normalize_pathway_name(
  msig_kegg$gs_name
)

msig_kegg$NormalizedGsDescription <- normalize_pathway_name(
  if (
    "gs_description" %in%
    colnames(
      msig_kegg
    )
  ) {
    msig_kegg$gs_description
  } else {
    msig_kegg$gs_name
  }
)

msig_kegg$TermID <- ifelse(
  !is.na(
    msig_kegg$KeggCode
  ),
  paste0(
    kegg_prefix,
    msig_kegg$KeggCode
  ),
  paste0(
    "MSIGDB:",
    msig_kegg$gs_name
  )
)

msig_kegg$GeneID <- sub(
  "\\.[0-9]+$",
  "",
  clean_text(
    msig_kegg$ensembl_gene
  )
)

msig_kegg <- msig_kegg[
  !is.na(
    msig_kegg$TermID
  ) &
    !is.na(
      msig_kegg$GeneID
    ) &
    grepl(
      valid_ensembl_pattern,
      msig_kegg$GeneID
    ),
  ,
  drop = FALSE
]

if (
  nrow(
    msig_kegg
  ) == 0
) {
  stop(
    "The MSigDB KEGG mapping produced no valid target-species Ensembl genes."
  )
}

kegg_term2gene <- unique(
  data.frame(
    Term = msig_kegg$TermID,
    Gene = msig_kegg$GeneID,
    stringsAsFactors = FALSE
  )
)

kegg_term_metadata <- unique(
  data.frame(
    Term = msig_kegg$TermID,
    KeggCode = msig_kegg$KeggCode,
    GsName = clean_text(
      msig_kegg$gs_name
    ),
    GsDescription = clean_text(
      if (
        "gs_description" %in%
        colnames(
          msig_kegg
        )
      ) {
        msig_kegg$gs_description
      } else {
        msig_kegg$gs_name
      }
    ),
    NormalizedGsName = msig_kegg$NormalizedGsName,
    NormalizedGsDescription = msig_kegg$NormalizedGsDescription,
    Subcollection = clean_text(
      msig_kegg$RequestedSubcollection
    ),
    stringsAsFactors = FALSE
  )
)

kegg_term_name_rows <- lapply(
  unique(
    kegg_term_metadata$Term
  ),
  function(term_id) {
    current_rows <- kegg_term_metadata[
      kegg_term_metadata$Term ==
        term_id,
      ,
      drop = FALSE
    ]
    
    current_names <- c(
      current_rows$GsDescription,
      current_rows$GsName
    )
    
    current_names <- unique(
      current_names[
        !is.na(
          current_names
        ) &
          current_names !=
          ""
      ]
    )
    
    data.frame(
      Term = term_id,
      Name = if (
        length(
          current_names
        ) > 0
      ) {
        current_names[1]
      } else {
        term_id
      },
      stringsAsFactors = FALSE
    )
  }
)

kegg_term2name <- do.call(
  rbind,
  kegg_term_name_rows
)

write_tsv(
  kegg_term2gene,
  file.path(
    output_dir,
    "00_MSigDB_KEGG_species_mapped_TERM2GENE.tsv"
  )
)

write_tsv(
  kegg_term_metadata,
  file.path(
    output_dir,
    "00_MSigDB_KEGG_term_metadata.tsv"
  )
)

write_tsv(
  kegg_term2name,
  file.path(
    output_dir,
    "00_MSigDB_KEGG_TERM2NAME.tsv"
  )
)


# =============================================================================
# 8. Run complete GO/KEGG GSEA
# =============================================================================

standardize_gsea <- function(
    object,
    comparison,
    database
) {
  empty <- data.frame(
    Comparison = character(0),
    Database = character(0),
    ID = character(0),
    Description = character(0),
    setSize = numeric(0),
    enrichmentScore = numeric(0),
    NES = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    qvalues = numeric(0),
    rank = numeric(0),
    leading_edge = character(0),
    core_enrichment = character(0),
    stringsAsFactors = FALSE
  )
  
  if (is.null(object)) {
    return(empty)
  }
  
  x <- as.data.frame(object)
  
  if (nrow(x) == 0) {
    return(empty)
  }
  
  data.frame(
    Comparison = comparison,
    Database = database,
    ID = clean_text(x$ID),
    Description = clean_text(x$Description),
    setSize = suppressWarnings(
      as.numeric(x$setSize)
    ),
    enrichmentScore = suppressWarnings(
      as.numeric(x$enrichmentScore)
    ),
    NES = suppressWarnings(
      as.numeric(x$NES)
    ),
    pvalue = suppressWarnings(
      as.numeric(x$pvalue)
    ),
    p.adjust = suppressWarnings(
      as.numeric(x$p.adjust)
    ),
    qvalues = if ("qvalues" %in% colnames(x)) {
      suppressWarnings(
        as.numeric(x$qvalues)
      )
    } else {
      NA_real_
    },
    rank = if ("rank" %in% colnames(x)) {
      suppressWarnings(
        as.numeric(x$rank)
      )
    } else {
      NA_real_
    },
    leading_edge = if ("leading_edge" %in% colnames(x)) {
      as.character(x$leading_edge)
    } else {
      NA_character_
    },
    core_enrichment = if ("core_enrichment" %in% colnames(x)) {
      as.character(x$core_enrichment)
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}


set.seed(
  gsea_seed
)

all_gsea_tables <- list()

for (comparison_name in comparison_order) {
  go_gene_list <- go_rank_lists[[comparison_name]]
  kegg_gene_list <- kegg_rank_lists[[comparison_name]]
  
  for (ontology in c("BP", "MF", "CC")) {
    database_name <- paste0(
      "GO_",
      ontology
    )
    
    gsea_object <- tryCatch(
      clusterProfiler::gseGO(
        geneList = go_gene_list,
        OrgDb = OrgDb,
        keyType = "ENSEMBL",
        ont = ontology,
        exponent = 1,
        minGSSize = minimum_gene_set_size,
        maxGSSize = maximum_gene_set_size,
        eps = gsea_eps,
        pvalueCutoff = gsea_pvalue_cutoff,
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = TRUE,
        by = "fgsea"
      ),
      error = function(e) {
        warning(
          comparison_name,
          " / ",
          database_name,
          " GSEA failed: ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    gsea_table <- standardize_gsea(
      gsea_object,
      comparison_name,
      database_name
    )
    
    write_tsv(
      gsea_table,
      file.path(
        complete_gsea_dir,
        paste0(
          gsub(
            " ",
            "_",
            comparison_name
          ),
          "_",
          database_name,
          "_complete_GSEA.tsv"
        )
      )
    )
    
    all_gsea_tables[[paste(
      comparison_name,
      database_name,
      sep = " | "
    )]] <- gsea_table
  }
  
  kegg_object <- tryCatch(
    clusterProfiler::GSEA(
      geneList = kegg_gene_list,
      TERM2GENE = kegg_term2gene,
      TERM2NAME = kegg_term2name,
      exponent = 1,
      minGSSize = minimum_gene_set_size,
      maxGSSize = maximum_gene_set_size,
      eps = gsea_eps,
      pvalueCutoff = gsea_pvalue_cutoff,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    ),
    error = function(e) {
      warning(
        comparison_name,
        " / KEGG GSEA failed: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  kegg_table <- standardize_gsea(
    kegg_object,
    comparison_name,
    "KEGG"
  )
  
  write_tsv(
    kegg_table,
    file.path(
      complete_gsea_dir,
      paste0(
        gsub(
          " ",
          "_",
          comparison_name
        ),
        "_KEGG_complete_GSEA.tsv"
      )
    )
  )
  
  all_gsea_tables[[paste(
    comparison_name,
    "KEGG",
    sep = " | "
  )]] <- kegg_table
}

all_complete_gsea <- do.call(
  rbind,
  all_gsea_tables
)

rownames(
  all_complete_gsea
) <- NULL

if (nrow(all_complete_gsea) == 0) {
  stop(
    "No complete GSEA results were generated."
  )
}

write_tsv(
  all_complete_gsea,
  file.path(
    complete_gsea_dir,
    "all_complete_GSEA_results.tsv"
  )
)


# =============================================================================
# 9. Read selected pathway IDs from the assigned Excel sheets
# =============================================================================

if (!file.exists(selected_pathway_excel)) {
  stop(
    "Cannot find pathway selector Excel:\n",
    selected_pathway_excel
  )
}

available_sheets <- readxl::excel_sheets(
  selected_pathway_excel
)

required_sheets <- names(
  selected_sheet_direction
)

missing_sheets <- setdiff(
  required_sheets,
  available_sheets
)

if (length(missing_sheets) > 0) {
  stop(
    "Excel is missing sheets: ",
    paste(
      missing_sheets,
      collapse = ", "
    )
  )
}

selected_list <- vector(
  "list",
  length(required_sheets)
)

ontology_check_list <- vector(
  "list",
  length(required_sheets)
)

for (sheet_index in seq_along(required_sheets)) {
  sheet_name <- required_sheets[sheet_index]
  
  selector_direction <- unname(
    selected_sheet_direction[
      sheet_name
    ]
  )
  
  x <- readxl::read_xlsx(
    selected_pathway_excel,
    sheet = sheet_name,
    col_names = TRUE,
    .name_repair = "minimal"
  )
  
  x <- as.data.frame(
    x,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cleaned_names <- clean_text(
    colnames(x)
  )
  
  cleaned_names[
    is.na(cleaned_names)
  ] <- "BlankColumn"
  
  colnames(x) <- make.unique(
    cleaned_names,
    sep = "_"
  )
  
  x <- remove_blank_columns(x)
  
  required_columns <- c(
    "ID",
    "Description"
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(x)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      sheet_name,
      " is missing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }
  
  id <- clean_text(
    x[
      ,
      "ID",
      drop = TRUE
    ]
  )
  
  database_from_excel <- normalize_database(
    database = if ("Database" %in% colnames(x)) {
      x$Database
    } else {
      NA
    },
    source_file = if ("SourceFile" %in% colnames(x)) {
      x$SourceFile
    } else {
      NA
    },
    panel = if ("Panel" %in% colnames(x)) {
      x$Panel
    } else {
      NA
    },
    id = id
  )
  
  database_from_annotation <- infer_go_database_from_id(
    id
  )
  
  database <- database_from_excel
  
  database[
    !is.na(
      database_from_annotation
    )
  ] <- database_from_annotation[
    !is.na(
      database_from_annotation
    )
  ]
  
  kegg_positions <- which(
    !is.na(database) &
      database == "KEGG" &
      !is.na(id)
  )
  
  if (length(kegg_positions) > 0) {
    id[
      kegg_positions
    ] <- normalize_kegg_pathway_id(
      id[
        kegg_positions
      ]
    )
  }
  
  if (any(is.na(database))) {
    stop(
      sheet_name,
      " contains unrecognized pathway IDs: ",
      paste(
        unique(
          id[
            is.na(database)
          ]
        ),
        collapse = ", "
      )
    )
  }
  
  description <- clean_text(
    x[
      ,
      "Description",
      drop = TRUE
    ]
  )
  
  default_label <- paste0(
    id,
    ", ",
    description
  )
  
  plot_label <- rep(
    NA_character_,
    nrow(x)
  )
  
  if (
    use_excel_display_label &
    "DisplayLabel" %in%
    colnames(x)
  ) {
    plot_label <- clean_text(
      x[
        ,
        "DisplayLabel",
        drop = TRUE
      ]
    )
  }
  
  plot_label[
    is.na(plot_label)
  ] <- default_label[
    is.na(plot_label)
  ]
  
  plot_order <- if (
    "PlotOrder" %in%
    colnames(x)
  ) {
    suppressWarnings(
      as.numeric(
        x[
          ,
          "PlotOrder",
          drop = TRUE
        ]
      )
    )
  } else {
    rep(
      NA_real_,
      nrow(x)
    )
  }
  
  source_row <- seq_len(
    nrow(x)
  ) + 1L
  
  plot_order[
    is.na(plot_order)
  ] <- source_row[
    is.na(plot_order)
  ]
  
  ontology_check_list[[sheet_index]] <- data.frame(
    SourceSheet = sheet_name,
    SourceExcelRow = source_row,
    ID = id,
    DatabaseFromExcelContext = database_from_excel,
    DatabaseFromCurrentGOAnnotation = database_from_annotation,
    DatabaseUsed = database,
    stringsAsFactors = FALSE
  )
  
  selected_list[[sheet_index]] <- data.frame(
    SourceSheet = sheet_name,
    SelectorDirection = selector_direction,
    Database = database,
    ID = id,
    SelectedDescription = description,
    PlotLabel = plot_label,
    PlotOrder = plot_order,
    SourceExcelRow = source_row,
    stringsAsFactors = FALSE
  )
}

selected_targets <- do.call(
  rbind,
  selected_list
)

ontology_check_table <- do.call(
  rbind,
  ontology_check_list
)

selected_targets <- selected_targets[
  !is.na(
    selected_targets$ID
  ) &
    !is.na(
      selected_targets$PlotLabel
    ),
  ,
  drop = FALSE
]

selected_key <- paste(
  selected_targets$SelectorDirection,
  selected_targets$Database,
  selected_targets$ID,
  sep = " | "
)

selected_targets <- selected_targets[
  !duplicated(
    selected_key
  ),
  ,
  drop = FALSE
]


# Map each selected KEGG row to the actual MSigDB KEGG term used for GSEA.
# GO rows retain their original GO IDs.
selected_targets$GSEAID <- selected_targets$ID
selected_targets$KEGGMatchMethod <- NA_character_
selected_targets$KEGGMatchedGeneSet <- NA_character_
selected_targets$KEGGMappedGenes <- NA_integer_


map_selected_kegg_term <- function(
    selected_id,
    selected_description
) {
  selected_id <- clean_text(
    selected_id
  )
  
  selected_description_normalized <- normalize_pathway_name(
    selected_description
  )
  
  # First choice: exact term ID, for example mmu04657.
  exact_id_hits <- unique(
    kegg_term_metadata$Term[
      kegg_term_metadata$Term ==
        selected_id
    ]
  )
  
  exact_id_hits <- exact_id_hits[
    !is.na(
      exact_id_hits
    )
  ]
  
  if (
    length(
      exact_id_hits
    ) > 0
  ) {
    return(
      list(
        Term = exact_id_hits[1],
        Method = "exact species KEGG ID"
      )
    )
  }
  
  selected_code <- sub(
    "^(hsa|mmu|ko)",
    "",
    tolower(
      selected_id
    )
  )
  
  if (
    grepl(
      "^[0-9]{5}$",
      selected_code
    )
  ) {
    code_hits <- unique(
      kegg_term_metadata$Term[
        !is.na(
          kegg_term_metadata$KeggCode
        ) &
          kegg_term_metadata$KeggCode ==
          selected_code
      ]
    )
    
    code_hits <- code_hits[
      !is.na(
        code_hits
      )
    ]
    
    if (
      length(
        code_hits
      ) > 0
    ) {
      return(
        list(
          Term = code_hits[1],
          Method = "exact KEGG numeric code"
        )
      )
    }
  }
  
  exact_name_hits <- unique(
    kegg_term_metadata$Term[
      kegg_term_metadata$NormalizedGsName ==
        selected_description_normalized |
        kegg_term_metadata$NormalizedGsDescription ==
        selected_description_normalized
    ]
  )
  
  exact_name_hits <- exact_name_hits[
    !is.na(
      exact_name_hits
    )
  ]
  
  if (
    length(
      exact_name_hits
    ) > 0
  ) {
    return(
      list(
        Term = exact_name_hits[1],
        Method = "exact normalized pathway name"
      )
    )
  }
  
  candidate_terms <- unique(
    kegg_term_metadata$Term
  )
  
  candidate_names <- vapply(
    candidate_terms,
    function(current_term) {
      current_rows <- kegg_term_metadata[
        kegg_term_metadata$Term ==
          current_term,
        ,
        drop = FALSE
      ]
      
      candidate_values <- c(
        current_rows$NormalizedGsName,
        current_rows$NormalizedGsDescription
      )
      
      candidate_values <- unique(
        candidate_values[
          !is.na(
            candidate_values
          ) &
            candidate_values !=
            ""
        ]
      )
      
      if (
        length(
          candidate_values
        ) == 0
      ) {
        return(
          ""
        )
      }
      
      candidate_values[1]
    },
    FUN.VALUE = character(1)
  )
  
  if (
    is.na(
      selected_description_normalized
    ) ||
    selected_description_normalized ==
    ""
  ) {
    return(
      list(
        Term = NA_character_,
        Method = "unmatched"
      )
    )
  }
  
  distance <- as.numeric(
    adist(
      selected_description_normalized,
      candidate_names
    )
  )
  
  denominator <- pmax(
    nchar(
      selected_description_normalized
    ),
    nchar(
      candidate_names
    ),
    1
  )
  
  normalized_distance <- distance /
    denominator
  
  best_index <- which.min(
    normalized_distance
  )
  
  if (
    length(
      best_index
    ) == 1 &&
    is.finite(
      normalized_distance[
        best_index
      ]
    ) &&
    normalized_distance[
      best_index
    ] <=
    0.30
  ) {
    return(
      list(
        Term = candidate_terms[
          best_index
        ],
        Method = paste0(
          "fuzzy pathway-name match; distance=",
          round(
            normalized_distance[
              best_index
            ],
            3
          )
        )
      )
    )
  }
  
  list(
    Term = NA_character_,
    Method = "unmatched"
  )
}


kegg_selected_positions <- which(
  selected_targets$Database ==
    "KEGG"
)

if (
  length(
    kegg_selected_positions
  ) > 0
) {
  for (
    position_index in kegg_selected_positions
  ) {
    mapping_result <- map_selected_kegg_term(
      selected_id = selected_targets$ID[
        position_index
      ],
      selected_description = selected_targets$SelectedDescription[
        position_index
      ]
    )
    
    selected_targets$GSEAID[
      position_index
    ] <- mapping_result$Term
    
    selected_targets$KEGGMatchMethod[
      position_index
    ] <- mapping_result$Method
    
    if (
      !is.na(
        mapping_result$Term
      )
    ) {
      metadata_rows <- kegg_term_metadata[
        kegg_term_metadata$Term ==
          mapping_result$Term,
        ,
        drop = FALSE
      ]
      
      matched_names <- unique(
        metadata_rows$GsName[
          !is.na(
            metadata_rows$GsName
          )
        ]
      )
      
      if (
        length(
          matched_names
        ) > 0
      ) {
        selected_targets$KEGGMatchedGeneSet[
          position_index
        ] <- matched_names[1]
      }
      
      selected_targets$KEGGMappedGenes[
        position_index
      ] <- length(
        unique(
          kegg_term2gene$Gene[
            kegg_term2gene$Term ==
              mapping_result$Term
          ]
        )
      )
    }
  }
}


write_tsv(
  selected_targets[
    selected_targets$Database ==
      "KEGG",
    ,
    drop = FALSE
  ],
  file.path(
    output_dir,
    "00_Sheet5_KEGG_mapping_diagnostic.tsv"
  )
)


write_tsv(
  selected_targets,
  file.path(
    output_dir,
    "00_selected_pathway_targets.tsv"
  )
)

write_tsv(
  ontology_check_table,
  file.path(
    output_dir,
    "00_GO_ontology_database_check.tsv"
  )
)


# =============================================================================
# 10. Build pathway gene maps and Sheet5 UP-DEG overlap counts
#
# GO gene sets come from the species-matched OrgDb.
# KEGG gene sets come from the MSigDB KEGG collections mapped to the detected
# target species. Both are represented with Ensembl IDs, so DEG Count is
# computed consistently for GO and KEGG.
# =============================================================================

selected_go_ids <- unique(
  selected_targets$GSEAID[
    grepl(
      "^GO_",
      selected_targets$Database
    ) &
      !is.na(
        selected_targets$GSEAID
      )
  ]
)

if (
  length(
    selected_go_ids
  ) > 0
) {
  valid_selected_go_ids <- intersect(
    selected_go_ids,
    AnnotationDbi::keys(
      OrgDb,
      keytype = "GOALL"
    )
  )
  
  go_pathway_gene_map <- suppressMessages(
    AnnotationDbi::select(
      OrgDb,
      keys = valid_selected_go_ids,
      keytype = "GOALL",
      columns = c(
        "ENSEMBL",
        "ONTOLOGYALL"
      )
    )
  )
  
  go_pathway_gene_map <- go_pathway_gene_map[
    !is.na(
      go_pathway_gene_map$GOALL
    ) &
      !is.na(
        go_pathway_gene_map$ENSEMBL
      ) &
      !is.na(
        go_pathway_gene_map$ONTOLOGYALL
      ),
    c(
      "GOALL",
      "ENSEMBL",
      "ONTOLOGYALL"
    ),
    drop = FALSE
  ]
  
  go_pathway_gene_map$ENSEMBL <- sub(
    "\\.[0-9]+$",
    "",
    clean_text(
      go_pathway_gene_map$ENSEMBL
    )
  )
  
  go_pathway_gene_map <- go_pathway_gene_map[
    grepl(
      valid_ensembl_pattern,
      go_pathway_gene_map$ENSEMBL
    ),
    ,
    drop = FALSE
  ]
  
  go_pathway_gene_map <- data.frame(
    Database = paste0(
      "GO_",
      go_pathway_gene_map$ONTOLOGYALL
    ),
    GSEAID = go_pathway_gene_map$GOALL,
    GeneID = go_pathway_gene_map$ENSEMBL,
    stringsAsFactors = FALSE
  )
  
  go_pathway_gene_map <- unique(
    go_pathway_gene_map
  )
  
} else {
  go_pathway_gene_map <- data.frame(
    Database = character(0),
    GSEAID = character(0),
    GeneID = character(0),
    stringsAsFactors = FALSE
  )
}


selected_kegg_gsea_ids <- unique(
  selected_targets$GSEAID[
    selected_targets$Database ==
      "KEGG" &
      !is.na(
        selected_targets$GSEAID
      )
  ]
)

kegg_pathway_gene_map <- kegg_term2gene[
  kegg_term2gene$Term %in%
    selected_kegg_gsea_ids,
  c(
    "Term",
    "Gene"
  ),
  drop = FALSE
]

if (
  nrow(
    kegg_pathway_gene_map
  ) > 0
) {
  kegg_pathway_gene_map <- data.frame(
    Database = "KEGG",
    GSEAID = kegg_pathway_gene_map$Term,
    GeneID = clean_text(
      kegg_pathway_gene_map$Gene
    ),
    stringsAsFactors = FALSE
  )
  
  kegg_pathway_gene_map <- unique(
    kegg_pathway_gene_map
  )
  
} else {
  kegg_pathway_gene_map <- data.frame(
    Database = character(0),
    GSEAID = character(0),
    GeneID = character(0),
    stringsAsFactors = FALSE
  )
}


selected_pathway_gene_map <- rbind(
  go_pathway_gene_map,
  kegg_pathway_gene_map
)

write_tsv(
  selected_pathway_gene_map,
  file.path(
    output_dir,
    "00_selected_pathway_gene_map_for_DEGCount.tsv"
  )
)


up_deg_ensembl_by_comparison <- list()

for (
  comparison_name in comparison_order
) {
  current_table <- all_rank_tables[[comparison_name]]
  
  up_rows <- (
    current_table$DirectionForCount ==
      "Up"
  )
  
  up_rows[
    is.na(
      up_rows
    )
  ] <- FALSE
  
  up_deg_ensembl_by_comparison[[comparison_name]] <- unique(
    current_table$EnsemblIDForGSEA[
      up_rows &
        !is.na(
          current_table$EnsemblIDForGSEA
        )
    ]
  )
}


calculate_up_deg_count <- function(
    comparison_name,
    database_name,
    gsea_id
) {
  if (
    is.na(
      gsea_id
    ) ||
    gsea_id ==
    ""
  ) {
    return(
      0L
    )
  }
  
  pathway_genes <- selected_pathway_gene_map$GeneID[
    selected_pathway_gene_map$Database ==
      database_name &
      selected_pathway_gene_map$GSEAID ==
      gsea_id
  ]
  
  pathway_genes <- unique(
    pathway_genes[
      !is.na(
        pathway_genes
      )
    ]
  )
  
  if (
    length(
      pathway_genes
    ) == 0
  ) {
    return(
      0L
    )
  }
  
  current_degs <- up_deg_ensembl_by_comparison[[comparison_name]]
  
  as.integer(
    length(
      intersect(
        pathway_genes,
        current_degs
      )
    )
  )
}


deg_count_summary <- do.call(
  rbind,
  lapply(
    comparison_order,
    function(comparison_name) {
      data.frame(
        Comparison = comparison_name,
        UpDEGEnsembl = length(
          up_deg_ensembl_by_comparison[[comparison_name]]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

write_tsv(
  deg_count_summary,
  file.path(
    output_dir,
    "00_UP_DEG_count_summary.tsv"
  )
)


# =============================================================================
# 11. Expand selected GO and KEGG pathways to all three comparisons
# =============================================================================

filled_rows <- vector(
  "list",
  nrow(
    selected_targets
  ) *
    length(
      comparison_order
    )
)

filled_index <- 0L

for (
  target_index in seq_len(
    nrow(
      selected_targets
    )
  )
) {
  target <- selected_targets[
    target_index,
    ,
    drop = FALSE
  ]
  
  target_gsea_id <- clean_text(
    target$GSEAID
  )
  
  for (
    comparison_name in comparison_order
  ) {
    filled_index <- filled_index + 1L
    
    if (
      is.na(
        target_gsea_id
      )
    ) {
      matched <- all_complete_gsea[
        FALSE,
        ,
        drop = FALSE
      ]
    } else {
      matched <- all_complete_gsea[
        all_complete_gsea$Comparison ==
          comparison_name &
          all_complete_gsea$Database ==
          target$Database &
          all_complete_gsea$ID ==
          target_gsea_id,
        ,
        drop = FALSE
      ]
    }
    
    if (
      nrow(
        matched
      ) > 1
    ) {
      matched <- matched[
        order(
          matched$p.adjust,
          -abs(
            matched$NES
          )
        ),
        ,
        drop = FALSE
      ]
      
      matched <- matched[
        1,
        ,
        drop = FALSE
      ]
    }
    
    deg_count <- calculate_up_deg_count(
      comparison_name = comparison_name,
      database_name = target$Database,
      gsea_id = target_gsea_id
    )
    
    if (
      nrow(
        matched
      ) == 1
    ) {
      filled_rows[[filled_index]] <- data.frame(
        Comparison = comparison_name,
        SelectorDirection = target$SelectorDirection,
        SourceSheet = target$SourceSheet,
        Database = target$Database,
        ID = target$ID,
        GSEAID = target_gsea_id,
        PlotLabel = target$PlotLabel,
        PlotOrder = target$PlotOrder,
        FoundInCompleteGSEA = TRUE,
        MissingReason = NA_character_,
        setSize = matched$setSize,
        enrichmentScore = matched$enrichmentScore,
        NES = matched$NES,
        pvalue = matched$pvalue,
        p.adjust = matched$p.adjust,
        qvalues = matched$qvalues,
        rank = matched$rank,
        leading_edge = matched$leading_edge,
        core_enrichment = matched$core_enrichment,
        DEGCount = deg_count,
        stringsAsFactors = FALSE
      )
      
    } else {
      missing_reason <- if (
        target$Database ==
        "KEGG" &&
        is.na(
          target_gsea_id
        )
      ) {
        "Selected KEGG pathway could not be matched to the MSigDB KEGG collection"
      } else {
        paste0(
          "Pathway not returned by complete GSEA; ",
          "check ranked-gene overlap and minGSSize/maxGSSize"
        )
      }
      
      filled_rows[[filled_index]] <- data.frame(
        Comparison = comparison_name,
        SelectorDirection = target$SelectorDirection,
        SourceSheet = target$SourceSheet,
        Database = target$Database,
        ID = target$ID,
        GSEAID = target_gsea_id,
        PlotLabel = target$PlotLabel,
        PlotOrder = target$PlotOrder,
        FoundInCompleteGSEA = FALSE,
        MissingReason = missing_reason,
        setSize = NA_real_,
        enrichmentScore = NA_real_,
        NES = NA_real_,
        pvalue = NA_real_,
        p.adjust = NA_real_,
        qvalues = NA_real_,
        rank = NA_real_,
        leading_edge = NA_character_,
        core_enrichment = NA_character_,
        DEGCount = deg_count,
        stringsAsFactors = FALSE
      )
    }
  }
}

plot_table <- do.call(
  rbind,
  filled_rows
)

plot_table$Comparison <- factor(
  plot_table$Comparison,
  levels = comparison_order
)

plot_table$NegLog10FDR <- ifelse(
  is.na(
    plot_table$p.adjust
  ),
  NA_real_,
  -log10(
    pmax(
      plot_table$p.adjust,
      .Machine$double.xmin
    )
  )
)

write_tsv(
  plot_table,
  file.path(
    plot_dir,
    "selected_pathway_GSEA_plot_data.tsv"
  )
)

write_tsv(
  plot_table[
    !plot_table$FoundInCompleteGSEA,
    ,
    drop = FALSE
  ],
  file.path(
    plot_dir,
    "missing_GSEA_diagnostics.tsv"
  )
)


# =============================================================================
# 12. Comparative Sheet5-selected GSEA plot
#
# Color = GSEA BH-adjusted P value
#   dark red = smaller adjusted P value
#   light red = larger adjusted P value
#
# Size = DEG Count
#   number of threshold-defined UP DEGs overlapping the pathway
#
# NES remains in selected_pathway_GSEA_plot_data.tsv and must be used to
# interpret the enrichment direction. The requested color scale does not encode
# whether NES is positive or negative.
# =============================================================================

positive_adjusted_p <- plot_table$p.adjust[
  is.finite(
    plot_table$p.adjust
  ) &
    plot_table$p.adjust >
    0
]

smallest_positive_adjusted_p <- if (
  length(
    positive_adjusted_p
  ) > 0
) {
  min(
    positive_adjusted_p
  )
} else {
  1e-300
}

plot_table$PAdjustForColour <- ifelse(
  !is.na(
    plot_table$p.adjust
  ) &
    plot_table$p.adjust <=
    0,
  smallest_positive_adjusted_p /
    10,
  plot_table$p.adjust
)

global_adjusted_p_values <- plot_table$PAdjustForColour[
  is.finite(
    plot_table$PAdjustForColour
  ) &
    plot_table$PAdjustForColour >
    0
]

if (
  length(
    global_adjusted_p_values
  ) == 0
) {
  global_colour_limits <- c(
    1e-300,
    1
  )
} else {
  global_colour_limits <- range(
    global_adjusted_p_values
  )
  
  if (
    global_colour_limits[1] ==
    global_colour_limits[2]
  ) {
    global_colour_limits <- c(
      global_colour_limits[1] /
        10,
      min(
        1,
        global_colour_limits[2] *
          10
      )
    )
  }
}

maximum_deg_count <- max(
  plot_table$DEGCount,
  na.rm = TRUE
)

if (
  !is.finite(
    maximum_deg_count
  ) ||
  maximum_deg_count <=
  0
) {
  maximum_deg_count <- 1
}

# Rewrite the plot table after all requested display variables are added.
write_tsv(
  plot_table,
  file.path(
    plot_dir,
    "selected_pathway_GSEA_plot_data.tsv"
  )
)

write_tsv(
  plot_table[
    !plot_table$FoundInCompleteGSEA,
    ,
    drop = FALSE
  ],
  file.path(
    plot_dir,
    "missing_GSEA_diagnostics.tsv"
  )
)


make_sheet4_gsea_plot <- function() {
  x <- plot_table
  
  if (
    nrow(
      x
    ) == 0
  ) {
    stop(
      "Sheet5 produced no selected pathways."
    )
  }
  
  order_table <- aggregate(
    PlotOrder ~ PlotLabel,
    data = x,
    FUN = min
  )
  
  order_table <- order_table[
    order(
      order_table$PlotOrder,
      order_table$PlotLabel
    ),
    ,
    drop = FALSE
  ]
  
  x$PlotLabel <- factor(
    x$PlotLabel,
    levels = rev(
      order_table$PlotLabel
    )
  )
  
  p <- ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = Comparison,
      y = PlotLabel
    )
  ) +
    ggplot2::geom_point(
      data = x[
        x$FoundInCompleteGSEA,
        ,
        drop = FALSE
      ],
      ggplot2::aes(
        size = DEGCount,
        colour = PAdjustForColour
      ),
      alpha = 1
    ) +
    ggplot2::geom_point(
      data = x[
        !x$FoundInCompleteGSEA,
        ,
        drop = FALSE
      ],
      shape = 1,
      size = 3,
      stroke = 0.8,
      colour = "grey55"
    ) +
    ggplot2::scale_colour_gradientn(
      colours = c(
        "#A50F15",
        "#CB181D",
        "#EF3B2C",
        "#FB6A4A",
        "#FC9272",
        "#FEE5D9"
      ),
      trans = "log10",
      limits = global_colour_limits,
      name = "Adjusted P value"
    ) +
    ggplot2::scale_size_continuous(
      range = c(
        3,
        10
      ),
      limits = c(
        0,
        maximum_deg_count
      ),
      breaks = scales::breaks_pretty(
        n = 5
      ),
      name = "DEG Count"
    ) +
    ggplot2::scale_x_discrete(
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Sheet5-selected GO/KEGG GSEA",
      subtitle = paste0(
        "Complete ranked-gene GSEA; GO = species OrgDb; KEGG = MSigDB mapped to mouse; ",
        "color = BH-adjusted P value; size = pathway-overlapping UP DEG count"
      ),
      caption = paste0(
        "NES is retained in the exported table for direction interpretation. ",
        "A count of 0 is drawn at the minimum size so a valid GSEA result remains visible. ",
        "Grey hollow marker = pathway still not returned after GO/MSigDB-KEGG mapping and size filtering."
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_bw(
      base_size = 11
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        colour = "grey30"
      ),
      plot.caption = ggplot2::element_text(
        hjust = 0,
        colour = "grey35"
      ),
      axis.text = ggplot2::element_text(
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
  
  height <- max(
    5.5,
    0.38 *
      length(
        unique(
          x$PlotLabel
        )
      ) +
      2.5
  )
  
  save_plot(
    p,
    "01_Sheet5_selected_GO_KEGG_GSEA_AdjustedP_DEGCount_dotplot",
    width = 11.5,
    height = height
  )
  
  p
}


sheet4_gsea_plot <- make_sheet4_gsea_plot()


# =============================================================================
# 13. Reproducibility records
# =============================================================================

settings <- data.frame(
  Setting = c(
    "Dataset",
    "DetectedSpecies",
    "OrgDbPackage",
    "IMPAllGenesTSV",
    "Y180AllGenesTSV",
    "Y180VsIMP321AllGenesTSV",
    "SelectedPathwayExcel",
    "SelectedSheet",
    "ColourMetric",
    "SizeMetric",
    "DEGpadjCutoff",
    "DEGabsShrunkenLog2FCCutoff",
    "RankMetric",
    "GOKeyType",
    "KEGGKeyType",
    "KEGGSource",
    "MSigDBVersion",
    "MSigDBKEGGSubcollections",
    "GSEApvalueCutoff",
    "GSEAeps",
    "MinimumGeneSetSize",
    "MaximumGeneSetSize"
  ),
  Value = c(
    dataset_name,
    analysis_species,
    orgdb_package,
    all_gene_files[["IMP vs Ctrl"]],
    all_gene_files[["Y180 vs Ctrl"]],
    all_gene_files[["Y180 vs IMP321"]],
    selected_pathway_excel,
    "Sheet5",
    "GSEA BH-adjusted P value",
    "Number of threshold-defined UP DEGs overlapping each pathway",
    deg_padj_cutoff,
    deg_abs_shrunken_log2fc_cutoff,
    "DESeq2 Wald statistic; signed -log10(P) only if stat is absent",
    "ENSEMBL",
    "ENSEMBL",
    "Human MSigDB KEGG collections mapped to target species by msigdbr",
    paste(
      unique(
        msig_kegg$db_version
      ),
      collapse = "; "
    ),
    paste(
      kegg_subcollections,
      collapse = "; "
    ),
    gsea_pvalue_cutoff,
    gsea_eps,
    minimum_gene_set_size,
    maximum_gene_set_size
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  settings,
  file.path(
    output_dir,
    "analysis_settings.tsv"
  )
)

writeLines(
  c(
    "GSEA uses every ranked tested gene; no DEG threshold is applied to calculate enrichment statistics.",
    "The primary ranking metric is the signed DESeq2 Wald statistic.",
    "GO GSEA uses the detected species OrgDb.",
    "KEGG GSEA uses the human MSigDB KEGG_LEGACY/KEGG_MEDICUS collections mapped to the detected target species by msigdbr.",
    "GO and KEGG GSEA are completed before Sheet5 pathway IDs are retained.",
    "Sheet5 is a display selector only; its GeneName column is not used as the ranked list.",
    "Bubble color is the GSEA BH-adjusted P value.",
    "Bubble size is the number of threshold-defined UP DEGs overlapping each pathway.",
    "The DEG-count threshold is separate from the GSEA significance calculation.",
    "NES is retained in the output table: positive means numerator-enriched and negative means denominator-enriched.",
    "For Y180 vs IMP321: positive NES means Y180-enriched; negative NES means IMP321-enriched.",
    "The KEGG mapping diagnostic is stored in 00_Sheet5_KEGG_mapping_diagnostic.tsv."
  ),
  file.path(
    output_dir,
    "README_GSEA_logic.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "R_sessionInfo.txt"
  )
)

message(
  "GSEA completed: ",
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = FALSE
  )
)

message(
  "Main plot: 02_selected_pathway_plots/01_Sheet5_selected_GO_KEGG_GSEA_AdjustedP_DEGCount_dotplot.pdf"
)


# =============================================================================
# PBMC: repair KEGG in the EXISTING Sheet5 GO/KEGG GSEA figure
#
# This script does NOT recalculate the GO part.
#
# It:
#   1. Reads the existing Sheet5 integrated GO/KEGG GSEA plot-data table.
#   2. Downloads/caches the complete Enrichr KEGG_2021_Human library.
#   3. Runs complete human KEGG GSEA for:
#        IMP vs Ctrl
#        Y180 vs Ctrl
#        Y180 vs IMP321
#   4. Replaces only KEGG rows in the existing combined table.
#   5. Puts KEGG pathway(s) at the TOP and GO pathways below.
#
# Plot definition preserved:
#   color = GSEA BH-adjusted P value
#   size  = pathway-overlapping UP DEG count
#
# For direction use NES in the exported table:
#   NES > 0 = numerator group enriched
#   NES < 0 = denominator group enriched
# =============================================================================


# =============================================================================
# 1. User settings
# =============================================================================

all_gene_files <- c(
  "IMP vs Ctrl" = paste0(
    "D:/个人记录/科研/Y180/PBMC Paired_GO_KEGG_results/",
    "IMP_vs_Ctrl/IMP_vs_Ctrl_paired_all_genes.tsv"
  ),
  "Y180 vs Ctrl" = paste0(
    "D:/个人记录/科研/Y180/PBMC Paired_GO_KEGG_results/",
    "Y180_vs_Ctrl/Y180_vs_Ctrl_paired_all_genes.tsv"
  ),
  "Y180 vs IMP321" = paste0(
    "D:/个人记录/科研/Y180/PBMC Paired_GO_KEGG_results/",
    "Y180_vs_IMP321/Y180_vs_IMP321_paired_all_genes.tsv"
  )
)

selected_excel <- paste0(
  "D:/个人记录/科研/Y180/",
  "PBMC Paired_GO_KEGG_results/final 10 selected.xlsx"
)

selected_sheet <- "Sheet5"

comparison_order <- c(
  "IMP vs Ctrl",
  "Y180 vs Ctrl",
  "Y180 vs IMP321"
)

# Directory produced by the existing PBMC integrated GSEA.
existing_output_dir <- paste0(
  "D:/个人记录/科研/Y180/",
  "PBMC Paired_GO_KEGG_results/",
  "selected_pathways_GSEA_Sheet5_AdjustedP_DEGCount_KEGG_integrated"
)

# Prefer the original Stage-1 plot table because it definitely contains GO rows.
existing_plot_data_file <- file.path(
  existing_output_dir,
  "02_selected_pathway_plots",
  "selected_pathway_GSEA_plot_data.tsv"
)

patched_plot_dir <- file.path(
  existing_output_dir,
  "03_KEGG2021_Human_patched_combined_plot"
)

dir.create(
  patched_plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Same DEG-count definition as before.
deg_padj_cutoff <- 0.05
deg_abs_shrunken_log2fc_cutoff <- 1

# GSEA settings.
minimum_gene_set_size <- 5L
maximum_gene_set_size <- 1000L
gsea_pvalue_cutoff <- 1
gsea_seed <- 20260827L

# Complete human KEGG library from Enrichr.
kegg_library_name <- "KEGG_2021_Human"

kegg_library_url <- paste0(
  "https://maayanlab.cloud/Enrichr/",
  "geneSetLibrary?mode=text&libraryName=",
  kegg_library_name
)

cached_kegg_file <- file.path(
  existing_output_dir,
  paste0(
    "00_",
    kegg_library_name,
    ".gmt"
  )
)

redownload_kegg_library <- FALSE

padjust_red_palette <- c(
  "#A50F15",
  "#CB181D",
  "#EF3B2C",
  "#FB6A4A",
  "#FC9272",
  "#FEE5D9"
)


# =============================================================================
# 2. Packages
# =============================================================================

cran_packages <- c(
  "readxl",
  "ggplot2",
  "scales"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0) {
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

bioc_packages <- c(
  "clusterProfiler",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_bioc <- bioc_packages[
  !vapply(
    bioc_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc) > 0) {
  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(scales)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})


# =============================================================================
# 3. Helpers
# =============================================================================

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("　", " ", x, fixed = TRUE)
  x <- trimws(x)
  x[x %in% c("", "NA", "N/A", "NULL", "-")] <- NA_character_
  x
}


normalize_name <- function(x) {
  x <- toupper(
    clean_text(x)
  )
  
  x <- gsub(
    "SIGNALLING",
    "SIGNALING",
    x
  )
  
  x <- gsub(
    "[^A-Z0-9]",
    "",
    x
  )
  
  x
}


write_tsv <- function(x, path) {
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
}


set_column_if_present <- function(
    data,
    row_index,
    column_name,
    value
) {
  if (column_name %in% colnames(data)) {
    data[row_index, column_name] <- value
  }
  
  data
}


# =============================================================================
# 4. Read existing PBMC GO/KEGG plot data
# =============================================================================

if (!file.exists(existing_plot_data_file)) {
  stop(
    paste0(
      "Cannot find existing PBMC plot data:\n",
      existing_plot_data_file,
      "\nPlease run PBMC_GSEA_SHEET5_GO_KEGG_INTEGRATED_COMPLETE.R first."
    )
  )
}

plot_table_old <- read.delim(
  existing_plot_data_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_plot_columns <- c(
  "Comparison",
  "Database",
  "ID",
  "PlotLabel",
  "PlotOrder",
  "FoundInCompleteGSEA",
  "p.adjust",
  "DEGCount"
)

missing_plot_columns <- setdiff(
  required_plot_columns,
  colnames(plot_table_old)
)

if (length(missing_plot_columns) > 0) {
  stop(
    paste0(
      "Existing plot table is missing columns:\n",
      paste(
        missing_plot_columns,
        collapse = ", "
      )
    )
  )
}

go_rows <- plot_table_old[
  plot_table_old$Database != "KEGG",
  ,
  drop = FALSE
]

old_kegg_rows <- plot_table_old[
  plot_table_old$Database == "KEGG",
  ,
  drop = FALSE
]

if (nrow(old_kegg_rows) == 0) {
  stop(
    "The existing Sheet5 plot table contains no KEGG rows."
  )
}


# =============================================================================
# 5. Download/cache KEGG_2021_Human
# =============================================================================

if (
  redownload_kegg_library ||
  !file.exists(cached_kegg_file)
) {
  options(
    timeout = max(
      600,
      getOption("timeout")
    )
  )
  
  options(
    download.file.method = "libcurl"
  )
  
  message(
    "Downloading ",
    kegg_library_name,
    " ..."
  )
  
  download_ok <- tryCatch(
    {
      download.file(
        url = kegg_library_url,
        destfile = cached_kegg_file,
        method = "libcurl",
        mode = "wb",
        quiet = FALSE
      )
      
      TRUE
    },
    error = function(e) {
      warning(
        "Download failed: ",
        conditionMessage(e)
      )
      
      FALSE
    }
  )
  
  if (
    !download_ok ||
    !file.exists(cached_kegg_file)
  ) {
    stop(
      paste0(
        "Could not download ",
        kegg_library_name,
        ".\n",
        "If necessary, manually download the library and save it as:\n",
        cached_kegg_file
      )
    )
  }
}

kegg_lines <- readLines(
  cached_kegg_file,
  warn = FALSE,
  encoding = "UTF-8"
)

kegg_lines <- kegg_lines[
  nzchar(
    trimws(
      kegg_lines
    )
  )
]

if (length(kegg_lines) == 0) {
  stop(
    "The KEGG_2021_Human library file is empty."
  )
}


# =============================================================================
# 6. Parse KEGG library into TERM2GENE
# =============================================================================

term2gene_rows <- vector(
  "list",
  length(kegg_lines)
)

term2name_rows <- vector(
  "list",
  length(kegg_lines)
)

for (line_index in seq_along(kegg_lines)) {
  fields <- strsplit(
    kegg_lines[line_index],
    split = "\t",
    fixed = TRUE
  )[[1]]
  
  fields <- clean_text(fields)
  
  pathway_name <- fields[1]
  
  pathway_genes <- unique(
    toupper(
      fields[
        !is.na(fields) &
          seq_along(fields) > 1
      ]
    )
  )
  
  pathway_genes <- pathway_genes[
    !is.na(pathway_genes) &
      pathway_genes != ""
  ]
  
  if (
    is.na(pathway_name) ||
    length(pathway_genes) == 0
  ) {
    next
  }
  
  term2gene_rows[[line_index]] <- data.frame(
    Term = pathway_name,
    Gene = pathway_genes,
    stringsAsFactors = FALSE
  )
  
  term2name_rows[[line_index]] <- data.frame(
    Term = pathway_name,
    Name = pathway_name,
    stringsAsFactors = FALSE
  )
}

term2gene_rows <- Filter(
  Negate(is.null),
  term2gene_rows
)

term2name_rows <- Filter(
  Negate(is.null),
  term2name_rows
)

if (length(term2gene_rows) == 0) {
  stop(
    "No valid KEGG_2021_Human pathway-to-gene mapping was parsed."
  )
}

kegg_term2gene <- unique(
  do.call(
    rbind,
    term2gene_rows
  )
)

kegg_term2name <- unique(
  do.call(
    rbind,
    term2name_rows
  )
)

write_tsv(
  kegg_term2gene,
  file.path(
    patched_plot_dir,
    "00_KEGG_2021_Human_TERM2GENE.tsv"
  )
)


# =============================================================================
# 7. Read Sheet5 KEGG pathway(s) and match by pathway name
# =============================================================================

if (!file.exists(selected_excel)) {
  stop(
    "Cannot find:\n",
    selected_excel
  )
}

selected_sheet_data <- as.data.frame(
  read_excel(
    selected_excel,
    sheet = selected_sheet
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(selected_sheet_data) <- make.unique(
  clean_text(
    colnames(selected_sheet_data)
  ),
  sep = "_"
)

if (
  !"ID" %in% colnames(selected_sheet_data) ||
  !"Description" %in% colnames(selected_sheet_data)
) {
  stop(
    "Sheet5 must contain ID and Description columns."
  )
}

selected_id <- clean_text(
  selected_sheet_data$ID
)

is_kegg_row <- grepl(
  "^(HSA|KO)?[0-9]{5}$",
  toupper(
    selected_id
  )
)

if ("Database" %in% colnames(selected_sheet_data)) {
  is_kegg_row <- (
    is_kegg_row |
      grepl(
        "KEGG",
        toupper(
          clean_text(
            selected_sheet_data$Database
          )
        )
      )
  )
}

selected_kegg <- selected_sheet_data[
  is_kegg_row &
    !is.na(selected_id),
  ,
  drop = FALSE
]

if (nrow(selected_kegg) == 0) {
  stop(
    "No KEGG pathway was found in Sheet5."
  )
}

selected_kegg$SelectedID <- clean_text(
  selected_kegg$ID
)

selected_kegg$SelectedDescription <- clean_text(
  selected_kegg$Description
)

selected_kegg$NormalizedDescription <- normalize_name(
  selected_kegg$SelectedDescription
)

library_name_map <- data.frame(
  Term = kegg_term2name$Term,
  NormalizedDescription = normalize_name(
    kegg_term2name$Name
  ),
  stringsAsFactors = FALSE
)

# ============================================================
# Robust Sheet5 KEGG matching
#
# Human PBMC:
#   hsa04151 = PI3K-Akt signaling pathway
#
# Matching order:
#   1. exact normalized name
#   2. explicit hsa04151 alias
#   3. PI3K/AKT keyword match
#   4. fuzzy normalized-name match
# ============================================================


# ------------------------------------------------------------
# 1. Normalize names
# ------------------------------------------------------------

normalize_kegg_name <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  # standardise spelling variants
  x <- gsub(
    "SIGNALLING",
    "SIGNALING",
    x
  )
  
  # remove species suffixes if present
  x <- gsub(
    "HOMO SAPIENS",
    "",
    x
  )
  
  x <- gsub(
    "HUMAN",
    "",
    x
  )
  
  # remove punctuation
  x <- gsub(
    "[^A-Z0-9]",
    "",
    x
  )
  
  x
}


# ------------------------------------------------------------
# 2. Prepare library names
# ------------------------------------------------------------

library_name_map <- data.frame(
  
  Term =
    as.character(
      kegg_term2name$Term
    ),
  
  NormalizedTerm =
    normalize_kegg_name(
      kegg_term2name$Term
    ),
  
  NormalizedName =
    normalize_kegg_name(
      kegg_term2name$Name
    ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 3. Prepare selected Sheet5 KEGG names
# ------------------------------------------------------------

selected_kegg$NormalizedDescription <-
  normalize_kegg_name(
    selected_kegg$SelectedDescription
  )


selected_kegg$MatchedTerm <- NA_character_

selected_kegg$MatchMethod <- NA_character_


# ------------------------------------------------------------
# 4. Match one selected KEGG pathway at a time
# ------------------------------------------------------------

for (
  row_index in seq_len(
    nrow(
      selected_kegg
    )
  )
) {
  
  current_id <-
    as.character(
      selected_kegg$SelectedID[
        row_index
      ]
    )
  
  
  current_description <-
    as.character(
      selected_kegg$SelectedDescription[
        row_index
      ]
    )
  
  
  current_normalized <-
    selected_kegg$NormalizedDescription[
      row_index
    ]
  
  
  matched_term <- NA_character_
  
  match_method <- NA_character_
  
  
  # ==========================================================
  # A. Exact normalized name
  # ==========================================================
  
  exact_hits <- library_name_map$Term[
    
    library_name_map$NormalizedName ==
      current_normalized |
      
      library_name_map$NormalizedTerm ==
      current_normalized
  ]
  
  
  exact_hits <- unique(
    exact_hits[
      !is.na(
        exact_hits
      )
    ]
  )
  
  
  if (
    length(
      exact_hits
    ) > 0
  ) {
    
    matched_term <- exact_hits[1]
    
    match_method <- "exact normalized pathway name"
  }
  
  
  # ==========================================================
  # B. Explicit hsa04151 safeguard
  # ==========================================================
  
  if (
    is.na(
      matched_term
    ) &&
    current_id ==
    "hsa04151"
  ) {
    
    pi3k_hits <- library_name_map$Term[
      
      grepl(
        "PI3K",
        library_name_map$Term,
        ignore.case = TRUE
      ) &
        
        grepl(
          "AKT",
          library_name_map$Term,
          ignore.case = TRUE
        )
    ]
    
    
    pi3k_hits <- unique(
      pi3k_hits[
        !is.na(
          pi3k_hits
        )
      ]
    )
    
    
    if (
      length(
        pi3k_hits
      ) > 0
    ) {
      
      matched_term <- pi3k_hits[1]
      
      match_method <-
        "explicit hsa04151 PI3K-AKT keyword match"
    }
  }
  
  
  # ==========================================================
  # C. Search normalized library names for PI3KAKT
  # ==========================================================
  
  if (
    is.na(
      matched_term
    ) &&
    grepl(
      "PI3K",
      current_normalized
    ) &&
    grepl(
      "AKT",
      current_normalized
    )
  ) {
    
    keyword_hits <- library_name_map$Term[
      
      grepl(
        "PI3K",
        library_name_map$NormalizedName
      ) &
        
        grepl(
          "AKT",
          library_name_map$NormalizedName
        )
    ]
    
    
    keyword_hits <- unique(
      keyword_hits[
        !is.na(
          keyword_hits
        )
      ]
    )
    
    
    if (
      length(
        keyword_hits
      ) > 0
    ) {
      
      matched_term <- keyword_hits[1]
      
      match_method <-
        "PI3K-AKT normalized keyword match"
    }
  }
  
  
  # ==========================================================
  # D. Fuzzy match
  # ==========================================================
  
  if (
    is.na(
      matched_term
    )
  ) {
    
    candidate_names <-
      library_name_map$NormalizedName
    
    
    distance <- adist(
      
      current_normalized,
      
      candidate_names
      
    )
    
    
    denominator <- pmax(
      
      nchar(
        current_normalized
      ),
      
      nchar(
        candidate_names
      ),
      
      1
    )
    
    
    normalized_distance <-
      as.numeric(
        distance
      ) /
      denominator
    
    
    best_index <-
      which.min(
        normalized_distance
      )
    
    
    if (
      length(
        best_index
      ) == 1 &&
      is.finite(
        normalized_distance[
          best_index
        ]
      ) &&
      normalized_distance[
        best_index
      ] <=
      0.35
    ) {
      
      matched_term <-
        library_name_map$Term[
          best_index
        ]
      
      
      match_method <-
        paste0(
          "fuzzy match; distance=",
          round(
            normalized_distance[
              best_index
            ],
            3
          )
        )
    }
  }
  
  
  selected_kegg$MatchedTerm[
    row_index
  ] <- matched_term
  
  
  selected_kegg$MatchMethod[
    row_index
  ] <- match_method
}


# ============================================================
# 5. Count pathway genes
# ============================================================

selected_kegg$MappedGeneCount <- vapply(
  
  selected_kegg$MatchedTerm,
  
  function(term_name) {
    
    if (
      is.na(
        term_name
      )
    ) {
      
      return(
        0L
      )
    }
    
    
    as.integer(
      
      length(
        
        unique(
          
          kegg_term2gene$Gene[
            
            kegg_term2gene$Term ==
              term_name
          ]
        )
      )
    )
  },
  
  FUN.VALUE = integer(1)
)


# ============================================================
# 6. Diagnostic table
# ============================================================

selected_mapping <- selected_kegg[
  ,
  c(
    "SelectedID",
    "SelectedDescription",
    "MatchedTerm",
    "MatchMethod",
    "MappedGeneCount"
  ),
  drop = FALSE
]


write.table(
  
  selected_mapping,
  
  file = file.path(
    patched_plot_dir,
    "00_Sheet5_KEGG_2021_Human_mapping.tsv"
  ),
  
  sep = "\t",
  
  quote = FALSE,
  
  row.names = FALSE
  
)


cat(
  "\n================ SHEET5 KEGG MAPPING ================\n"
)


print(
  selected_mapping,
  row.names = FALSE
)


# ============================================================
# 7. Do NOT stop immediately
# ============================================================

if (
  any(
    is.na(
      selected_mapping$MatchedTerm
    )
  )
) {
  
  warning(
    paste0(
      "Some Sheet5 KEGG pathways remain unmatched. ",
      "The script will continue. See ",
      "00_Sheet5_KEGG_2021_Human_mapping.tsv"
    )
  )
}


# ============================================================
# 8. Specific check for PI3K-Akt
# ============================================================

pi3k_check <- selected_mapping[
  selected_mapping$SelectedID ==
    "hsa04151",
  ,
  drop = FALSE
]


cat(
  "\n================ hsa04151 CHECK ================\n"
)


print(
  pi3k_check,
  row.names = FALSE
)


if (
  nrow(
    pi3k_check
  ) == 0
) {
  
  stop(
    "hsa04151 was not found in Sheet5."
  )
}


if (
  is.na(
    pi3k_check$MatchedTerm[1]
  )
) {
  
  cat(
    "\nAvailable KEGG library terms containing PI3K or AKT:\n"
  )
  
  
  candidate_terms <- library_name_map$Term[
    
    grepl(
      "PI3K|AKT",
      library_name_map$Term,
      ignore.case = TRUE
    )
  ]
  
  
  print(
    unique(
      candidate_terms
    )
  )
  
  
  stop(
    paste0(
      "hsa04151 still could not be matched. ",
      "The PI3K/AKT candidate terms were printed above."
    )
  )
}


# =============================================================================
# 8. Read three complete PBMC DE tables and build SYMBOL rankings
# =============================================================================

read_complete_de_table <- function(
    file_path,
    comparison_name
) {
  if (!file.exists(file_path)) {
    stop(
      "Cannot find:\n",
      file_path
    )
  }
  
  x <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  
  colnames(x) <- make.unique(
    clean_text(
      colnames(x)
    ),
    sep = "_"
  )
  
  ensembl_column <- if (
    "EnsemblID_clean" %in% colnames(x)
  ) {
    "EnsemblID_clean"
  } else if (
    "Geneid" %in% colnames(x)
  ) {
    "Geneid"
  } else {
    stop(
      comparison_name,
      " has no EnsemblID_clean/Geneid column."
    )
  }
  
  x$EnsemblIDForMap <- sub(
    "\\.[0-9]+$",
    "",
    clean_text(
      x[
        ,
        ensembl_column,
        drop = TRUE
      ]
    )
  )
  
  x <- x[
    !is.na(
      x$EnsemblIDForMap
    ) &
      grepl(
        "^ENSG[0-9]+$",
        x$EnsemblIDForMap
      ),
    ,
    drop = FALSE
  ]
  
  if ("GeneName" %in% colnames(x)) {
    symbol_value <- clean_text(
      x$GeneName
    )
  } else {
    symbol_value <- rep(
      NA_character_,
      nrow(x)
    )
  }
  
  missing_symbol <- (
    is.na(symbol_value) |
      symbol_value ==
      ""
  )
  
  if (any(missing_symbol)) {
    symbol_map <- suppressMessages(
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = unique(
          x$EnsemblIDForMap[
            missing_symbol
          ]
        ),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
    
    symbol_value[
      missing_symbol
    ] <- unname(
      symbol_map[
        x$EnsemblIDForMap[
          missing_symbol
        ]
      ]
    )
  }
  
  x$SymbolForGSEA <- toupper(
    clean_text(
      symbol_value
    )
  )
  
  if ("stat" %in% colnames(x)) {
    x$RankScore <- suppressWarnings(
      as.numeric(
        x$stat
      )
    )
    
    rank_method <- "DESeq2 Wald statistic"
    
  } else if (
    all(
      c(
        "log2FoldChange",
        "pvalue"
      ) %in% colnames(x)
    )
  ) {
    lfc <- suppressWarnings(
      as.numeric(
        x$log2FoldChange
      )
    )
    
    p_value <- suppressWarnings(
      as.numeric(
        x$pvalue
      )
    )
    
    x$RankScore <- sign(lfc) *
      -log10(
        pmax(
          p_value,
          .Machine$double.xmin
        )
      )
    
    rank_method <- "signed -log10(P value) fallback"
    
  } else {
    stop(
      comparison_name,
      " has neither stat nor log2FoldChange+pvalue."
    )
  }
  
  x <- x[
    !is.na(
      x$SymbolForGSEA
    ) &
      x$SymbolForGSEA !=
      "" &
      is.finite(
        x$RankScore
      ),
    ,
    drop = FALSE
  ]
  
  x <- x[
    order(
      -abs(
        x$RankScore
      )
    ),
    ,
    drop = FALSE
  ]
  
  x <- x[
    !duplicated(
      x$SymbolForGSEA
    ),
    ,
    drop = FALSE
  ]
  
  if ("Direction" %in% colnames(x)) {
    direction_upper <- toupper(
      clean_text(
        x$Direction
      )
    )
    
    x$DirectionForCount <- "Not_significant"
    x$DirectionForCount[
      direction_upper ==
        "UP"
    ] <- "Up"
    x$DirectionForCount[
      direction_upper ==
        "DOWN"
    ] <- "Down"
    
  } else {
    required_direction_columns <- c(
      "padj",
      "log2FoldChange_shrunken"
    )
    
    missing_direction_columns <- setdiff(
      required_direction_columns,
      colnames(x)
    )
    
    if (length(missing_direction_columns) > 0) {
      stop(
        comparison_name,
        " cannot reconstruct Direction because columns are missing: ",
        paste(
          missing_direction_columns,
          collapse = ", "
        )
      )
    }
    
    padj_value <- suppressWarnings(
      as.numeric(
        x$padj
      )
    )
    
    shrunken_lfc <- suppressWarnings(
      as.numeric(
        x$log2FoldChange_shrunken
      )
    )
    
    x$DirectionForCount <- "Not_significant"
    
    x$DirectionForCount[
      !is.na(padj_value) &
        padj_value <
        deg_padj_cutoff &
        !is.na(shrunken_lfc) &
        shrunken_lfc >=
        deg_abs_shrunken_log2fc_cutoff
    ] <- "Up"
    
    x$DirectionForCount[
      !is.na(padj_value) &
        padj_value <
        deg_padj_cutoff &
        !is.na(shrunken_lfc) &
        shrunken_lfc <=
        -deg_abs_shrunken_log2fc_cutoff
    ] <- "Down"
  }
  
  x$Comparison <- comparison_name
  x$RankMethod <- rank_method
  
  x
}


de_tables <- lapply(
  names(all_gene_files),
  function(comparison_name) {
    read_complete_de_table(
      all_gene_files[[comparison_name]],
      comparison_name
    )
  }
)

names(de_tables) <- names(all_gene_files)


rank_lists <- lapply(
  de_tables,
  function(x) {
    rank_vector <- x$RankScore
    names(rank_vector) <- x$SymbolForGSEA
    
    sort(
      rank_vector,
      decreasing = TRUE
    )
  }
)


up_deg_lists <- lapply(
  de_tables,
  function(x) {
    unique(
      x$SymbolForGSEA[
        x$DirectionForCount ==
          "Up"
      ]
    )
  }
)


rank_summary <- do.call(
  rbind,
  lapply(
    comparison_order,
    function(comparison_name) {
      data.frame(
        Comparison = comparison_name,
        RankedSymbols = length(
          rank_lists[[comparison_name]]
        ),
        UpDEGSymbols = length(
          up_deg_lists[[comparison_name]]
        ),
        RankMethod = unique(
          de_tables[[comparison_name]]$RankMethod
        )[1],
        stringsAsFactors = FALSE
      )
    }
  )
)

write_tsv(
  rank_summary,
  file.path(
    patched_plot_dir,
    "00_PBMC_KEGG_rank_and_UP_DEG_summary.tsv"
  )
)


# =============================================================================
# 9. Run complete KEGG_2021_Human GSEA
# =============================================================================

run_kegg_gsea <- function(
    ranked_gene_list,
    comparison_name
) {
  gsea_arguments <- list(
    geneList = ranked_gene_list,
    exponent = 1,
    minGSSize = minimum_gene_set_size,
    maxGSSize = maximum_gene_set_size,
    pvalueCutoff = gsea_pvalue_cutoff,
    pAdjustMethod = "BH",
    TERM2GENE = kegg_term2gene,
    TERM2NAME = kegg_term2name,
    verbose = FALSE
  )
  
  available_arguments <- names(
    formals(
      clusterProfiler::GSEA
    )
  )
  
  if ("eps" %in% available_arguments) {
    gsea_arguments$eps <- 0
  }
  
  if ("seed" %in% available_arguments) {
    gsea_arguments$seed <- TRUE
  }
  
  if ("by" %in% available_arguments) {
    gsea_arguments$by <- "fgsea"
  }
  
  set.seed(
    gsea_seed
  )
  
  result_object <- tryCatch(
    do.call(
      clusterProfiler::GSEA,
      gsea_arguments
    ),
    error = function(e) {
      stop(
        comparison_name,
        " KEGG_2021_Human GSEA failed: ",
        conditionMessage(e)
      )
    }
  )
  
  result_table <- as.data.frame(
    result_object
  )
  
  if (nrow(result_table) > 0) {
    result_table$Comparison <- comparison_name
  }
  
  result_table
}


complete_kegg_results <- list()

for (comparison_name in comparison_order) {
  message(
    "Running ",
    kegg_library_name,
    " GSEA: ",
    comparison_name
  )
  
  complete_kegg_results[[comparison_name]] <- run_kegg_gsea(
    rank_lists[[comparison_name]],
    comparison_name
  )
}

all_complete_kegg_gsea <- do.call(
  rbind,
  complete_kegg_results
)

rownames(
  all_complete_kegg_gsea
) <- NULL

if (nrow(all_complete_kegg_gsea) == 0) {
  stop(
    "KEGG_2021_Human GSEA returned no pathways."
  )
}

write_tsv(
  all_complete_kegg_gsea,
  file.path(
    patched_plot_dir,
    "all_complete_KEGG_2021_Human_GSEA.tsv"
  )
)


# =============================================================================
# 10. Patch KEGG row(s) in the existing plot table
# =============================================================================

pathway_gene_lists <- split(
  kegg_term2gene$Gene,
  kegg_term2gene$Term
)

pathway_gene_lists <- lapply(
  pathway_gene_lists,
  unique
)


calculate_deg_count <- function(
    comparison_name,
    pathway_term
) {
  pathway_genes <- pathway_gene_lists[[pathway_term]]
  
  if (is.null(pathway_genes)) {
    return(0L)
  }
  
  as.integer(
    length(
      intersect(
        pathway_genes,
        up_deg_lists[[comparison_name]]
      )
    )
  )
}


patched_kegg_rows <- old_kegg_rows

if (!"KEGGSourceTerm" %in% colnames(patched_kegg_rows)) {
  patched_kegg_rows$KEGGSourceTerm <- NA_character_
}

for (row_index in seq_len(nrow(patched_kegg_rows))) {
  comparison_name <- as.character(
    patched_kegg_rows$Comparison[row_index]
  )
  
  selected_pathway_id <- as.character(
    patched_kegg_rows$ID[row_index]
  )
  
  pathway_term <- selected_kegg$MatchedTerm[
    match(
      selected_pathway_id,
      selected_kegg$SelectedID
    )
  ]
  
  if (
    length(pathway_term) == 0 ||
    is.na(pathway_term)
  ) {
    patched_kegg_rows$FoundInCompleteGSEA[row_index] <- FALSE
    patched_kegg_rows$MissingReason[row_index] <- paste0(
      selected_pathway_id,
      " was not found in ",
      kegg_library_name
    )
    
    next
  }
  
  current_result <- all_complete_kegg_gsea[
    all_complete_kegg_gsea$Comparison ==
      comparison_name &
      all_complete_kegg_gsea$ID ==
      pathway_term,
    ,
    drop = FALSE
  ]
  
  if (nrow(current_result) > 1) {
    current_result <- current_result[
      order(
        current_result$p.adjust,
        -abs(
          current_result$NES
        )
      ),
      ,
      drop = FALSE
    ]
    
    current_result <- current_result[
      1,
      ,
      drop = FALSE
    ]
  }
  
  deg_count <- calculate_deg_count(
    comparison_name,
    pathway_term
  )
  
  patched_kegg_rows$DEGCount[row_index] <- deg_count
  patched_kegg_rows$KEGGSourceTerm[row_index] <- pathway_term
  
  if (nrow(current_result) == 1) {
    patched_kegg_rows$FoundInCompleteGSEA[row_index] <- TRUE
    patched_kegg_rows$MissingReason[row_index] <- NA_character_
    
    patched_kegg_rows <- set_column_if_present(
      patched_kegg_rows,
      row_index,
      "GSEAID",
      selected_pathway_id
    )
    
    patched_kegg_rows <- set_column_if_present(
      patched_kegg_rows,
      row_index,
      "setSize",
      suppressWarnings(
        as.numeric(
          current_result$setSize[1]
        )
      )
    )
    
    patched_kegg_rows <- set_column_if_present(
      patched_kegg_rows,
      row_index,
      "enrichmentScore",
      suppressWarnings(
        as.numeric(
          current_result$enrichmentScore[1]
        )
      )
    )
    
    patched_kegg_rows <- set_column_if_present(
      patched_kegg_rows,
      row_index,
      "NES",
      suppressWarnings(
        as.numeric(
          current_result$NES[1]
        )
      )
    )
    
    patched_kegg_rows <- set_column_if_present(
      patched_kegg_rows,
      row_index,
      "pvalue",
      suppressWarnings(
        as.numeric(
          current_result$pvalue[1]
        )
      )
    )
    
    patched_kegg_rows <- set_column_if_present(
      patched_kegg_rows,
      row_index,
      "p.adjust",
      suppressWarnings(
        as.numeric(
          current_result$p.adjust[1]
        )
      )
    )
    
    if ("qvalues" %in% colnames(current_result)) {
      patched_kegg_rows <- set_column_if_present(
        patched_kegg_rows,
        row_index,
        "qvalues",
        suppressWarnings(
          as.numeric(
            current_result$qvalues[1]
          )
        )
      )
    }
    
    if ("rank" %in% colnames(current_result)) {
      patched_kegg_rows <- set_column_if_present(
        patched_kegg_rows,
        row_index,
        "rank",
        suppressWarnings(
          as.numeric(
            current_result$rank[1]
          )
        )
      )
    }
    
    if ("leading_edge" %in% colnames(current_result)) {
      patched_kegg_rows <- set_column_if_present(
        patched_kegg_rows,
        row_index,
        "leading_edge",
        as.character(
          current_result$leading_edge[1]
        )
      )
    }
    
    if ("core_enrichment" %in% colnames(current_result)) {
      patched_kegg_rows <- set_column_if_present(
        patched_kegg_rows,
        row_index,
        "core_enrichment",
        as.character(
          current_result$core_enrichment[1]
        )
      )
    }
    
  } else {
    patched_kegg_rows$FoundInCompleteGSEA[row_index] <- FALSE
    patched_kegg_rows$MissingReason[row_index] <- paste0(
      pathway_term,
      " was not returned after ranked-gene overlap/size filtering"
    )
  }
}


# =============================================================================
# 11. Merge unchanged GO + repaired KEGG
# =============================================================================

if (!"KEGGSourceTerm" %in% colnames(go_rows)) {
  go_rows$KEGGSourceTerm <- NA_character_
}

missing_in_go <- setdiff(
  colnames(patched_kegg_rows),
  colnames(go_rows)
)

for (column_name in missing_in_go) {
  go_rows[, column_name] <- NA
}

missing_in_kegg <- setdiff(
  colnames(go_rows),
  colnames(patched_kegg_rows)
)

for (column_name in missing_in_kegg) {
  patched_kegg_rows[, column_name] <- NA
}

go_rows <- go_rows[
  ,
  colnames(patched_kegg_rows),
  drop = FALSE
]

patched_plot_table <- rbind(
  go_rows,
  patched_kegg_rows
)

patched_plot_table$Comparison <- factor(
  patched_plot_table$Comparison,
  levels = comparison_order
)


# =============================================================================
# 12. Put KEGG at TOP and GO below
# =============================================================================

term_info <- unique(
  patched_plot_table[
    ,
    c(
      "Database",
      "PlotLabel",
      "PlotOrder"
    ),
    drop = FALSE
  ]
)

term_info$PathwayClass <- ifelse(
  term_info$Database == "KEGG",
  "KEGG",
  "GO"
)

go_terms <- term_info[
  term_info$PathwayClass == "GO",
  ,
  drop = FALSE
]

kegg_terms <- term_info[
  term_info$PathwayClass == "KEGG",
  ,
  drop = FALSE
]

go_terms <- go_terms[
  order(
    go_terms$PlotOrder,
    go_terms$PlotLabel
  ),
  ,
  drop = FALSE
]

kegg_terms <- kegg_terms[
  order(
    kegg_terms$PlotOrder,
    kegg_terms$PlotLabel
  ),
  ,
  drop = FALSE
]

go_levels <- rev(
  go_terms$PlotLabel
)

kegg_levels <- rev(
  kegg_terms$PlotLabel
)

y_levels <- c(
  go_levels,
  kegg_levels
)

patched_plot_table$PlotLabel <- factor(
  patched_plot_table$PlotLabel,
  levels = y_levels
)


# =============================================================================
# 13. Color/size scales
# =============================================================================

positive_adjusted_p <- patched_plot_table$p.adjust[
  is.finite(
    patched_plot_table$p.adjust
  ) &
    patched_plot_table$p.adjust > 0
]

smallest_positive <- if (
  length(
    positive_adjusted_p
  ) > 0
) {
  min(
    positive_adjusted_p
  )
} else {
  1e-300
}

patched_plot_table$PAdjustForColour <- ifelse(
  !is.na(
    patched_plot_table$p.adjust
  ) &
    patched_plot_table$p.adjust <= 0,
  smallest_positive / 10,
  patched_plot_table$p.adjust
)

valid_colour_values <- patched_plot_table$PAdjustForColour[
  is.finite(
    patched_plot_table$PAdjustForColour
  ) &
    patched_plot_table$PAdjustForColour > 0
]

colour_limits <- if (
  length(
    valid_colour_values
  ) == 0
) {
  c(
    1e-300,
    1
  )
} else {
  range(
    valid_colour_values
  )
}

if (
  colour_limits[1] ==
  colour_limits[2]
) {
  colour_limits <- c(
    colour_limits[1] / 10,
    min(
      1,
      colour_limits[2] * 10
    )
  )
}

maximum_deg_count <- max(
  patched_plot_table$DEGCount,
  na.rm = TRUE
)

if (
  !is.finite(
    maximum_deg_count
  ) ||
  maximum_deg_count <= 0
) {
  maximum_deg_count <- 1
}


# =============================================================================
# 14. Draw repaired combined figure
# =============================================================================

p_final <- ggplot(
  patched_plot_table,
  aes(
    x = Comparison,
    y = PlotLabel
  )
) +
  geom_point(
    data = patched_plot_table[
      patched_plot_table$FoundInCompleteGSEA,
      ,
      drop = FALSE
    ],
    aes(
      size = DEGCount,
      colour = PAdjustForColour
    ),
    alpha = 1
  ) +
  geom_point(
    data = patched_plot_table[
      !patched_plot_table$FoundInCompleteGSEA,
      ,
      drop = FALSE
    ],
    shape = 1,
    size = 3,
    stroke = 0.8,
    colour = "grey55"
  ) +
  scale_colour_gradientn(
    colours = padjust_red_palette,
    trans = "log10",
    limits = colour_limits,
    name = "Adjusted P value"
  ) +
  scale_size_continuous(
    range = c(
      3,
      10
    ),
    limits = c(
      0,
      maximum_deg_count
    ),
    breaks = scales::breaks_pretty(
      n = 5
    ),
    name = "DEG Count"
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  labs(
    title = "Sheet5-selected GO/KEGG GSEA",
    subtitle = paste0(
      "GO = human OrgDb; KEGG = Enrichr KEGG_2021_Human; ",
      "color = BH-adjusted P value; size = pathway-overlapping UP DEG count"
    ),
    caption = paste0(
      "GSEA uses all ranked genes. NES remains in the exported table for direction. ",
      "A DEG count of 0 is drawn at minimum size."
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      colour = "grey30"
    ),
    plot.caption = element_text(
      hjust = 0,
      colour = "grey35"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

plot_height <- max(
  5.5,
  0.38 *
    length(
      unique(
        patched_plot_table$PlotLabel
      )
    ) +
    2.5
)

print(
  p_final
)


# =============================================================================
# 15. Save
# =============================================================================

write_tsv(
  patched_plot_table,
  file.path(
    patched_plot_dir,
    "selected_pathway_GSEA_plot_data_KEGG2021_Human_patched.tsv"
  )
)

write_tsv(
  patched_plot_table[
    !patched_plot_table$FoundInCompleteGSEA,
    ,
    drop = FALSE
  ],
  file.path(
    patched_plot_dir,
    "missing_GSEA_diagnostics_KEGG2021_Human_patched.tsv"
  )
)

ggsave(
  filename = file.path(
    patched_plot_dir,
    "01_Sheet5_GO_KEGG_GSEA_KEGG2021_Human_FIXED.png"
  ),
  plot = p_final,
  width = 12,
  height = plot_height,
  units = "in",
  dpi = 320,
  bg = "white"
)

ggsave(
  filename = file.path(
    patched_plot_dir,
    "01_Sheet5_GO_KEGG_GSEA_KEGG2021_Human_FIXED.pdf"
  ),
  plot = p_final,
  width = 12,
  height = plot_height,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)


# =============================================================================
# 16. Final KEGG diagnostic
# =============================================================================

kegg_check_columns <- intersect(
  c(
    "Comparison",
    "ID",
    "KEGGSourceTerm",
    "FoundInCompleteGSEA",
    "setSize",
    "NES",
    "pvalue",
    "p.adjust",
    "DEGCount",
    "MissingReason"
  ),
  colnames(patched_kegg_rows)
)

cat(
  "\n================ FINAL PBMC KEGG RESULTS ================\n"
)

print(
  patched_kegg_rows[
    ,
    kegg_check_columns,
    drop = FALSE
  ],
  row.names = FALSE
)

cat(
  "\nFinal figure:\n",
  file.path(
    patched_plot_dir,
    "01_Sheet5_GO_KEGG_GSEA_KEGG2021_Human_FIXED.pdf"
  ),
  "\n",
  sep = ""
)
