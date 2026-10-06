# =============================================================================
# BMDC COMPLETE GSEA PIPELINE
# Starts from complete DESeq2 all-gene tables and generates final Sheet4
# GO/KEGG GSEA bubble plot.
#
# GO   = gseGO with org.Mm.eg.db
# KEGG = Enrichr KEGG_2019_Mouse + clusterProfiler::GSEA
# color = BH-adjusted P value
# size  = pathway-overlapping UP DEG Count
# KEGG pathways at top; GO below
# =============================================================================

# -----------------------------
# 1. User settings
# -----------------------------

all_gene_files <- c(
  "IMP vs Ctrl" =
    "D:/个人记录/科研/Y180/BMDC Paired_GO_KEGG_results/IMP_vs_Ctrl/IMP_vs_Ctrl_paired_all_genes.tsv",
  "Y180 vs Ctrl" =
    "D:/个人记录/科研/Y180/BMDC Paired_GO_KEGG_results/Y180_vs_Ctrl/Y180_vs_Ctrl_paired_all_genes.tsv",
  "Y180 vs IMP321" =
    "D:/个人记录/科研/Y180/BMDC Paired_GO_KEGG_results/Y180_vs_IMP321/Y180_vs_IMP321_paired_all_genes.tsv"
)

selected_pathway_excel <-
  "D:/个人记录/科研/Y180/BMDC Paired_GO_KEGG_results/final 10 selected BMDC.xlsx"

selected_sheet <- "Sheet4"

comparison_order <- c(
  "IMP vs Ctrl",
  "Y180 vs Ctrl",
  "Y180 vs IMP321"
)

deg_padj_cutoff <- 0.05
deg_abs_shrunken_log2fc_cutoff <- 1

minimum_gene_set_size <- 5L
maximum_gene_set_size <- 1000L
gsea_pvalue_cutoff <- 1
gsea_seed <- 20260827L

kegg_library_name <- "KEGG_2019_Mouse"

kegg_library_url <- paste0(
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=",
  kegg_library_name
)

output_dir <-
  "D:/个人记录/科研/Y180/BMDC Paired_GO_KEGG_results/BMDC_GSEA_Sheet4_complete"

complete_gsea_dir <- file.path(output_dir, "01_complete_GSEA")
plot_dir <- file.path(output_dir, "02_selected_pathway_plots")

dir.create(complete_gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cached_kegg_file <- file.path(
  output_dir,
  paste0("00_", kegg_library_name, ".gmt")
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

# -----------------------------
# 2. Packages
# -----------------------------

cran_packages <- c("readxl", "ggplot2", "scales")

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "clusterProfiler",
  "AnnotationDbi",
  "org.Mm.eg.db"
)

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(scales)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

OrgDb <- org.Mm.eg.db::org.Mm.eg.db

# -----------------------------
# 3. Helpers
# -----------------------------

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

normalize_database <- function(database, source_file, panel, id) {
  database <- clean_text(database)
  source_file <- clean_text(source_file)
  panel <- clean_text(panel)
  id <- clean_text(id)
  
  context <- toupper(paste(database, source_file, panel, sep = " | "))
  
  answer <- rep(NA_character_, length(id))
  
  answer[grepl("^(MMU|KO)[0-9]+", toupper(id))] <- "KEGG"
  answer[is.na(answer) & grepl("KEGG", context)] <- "KEGG"
  answer[is.na(answer) & grepl("GO_?BP|BIOLOGICAL_PROCESS", context)] <- "GO_BP"
  answer[is.na(answer) & grepl("GO_?MF|MOLECULAR_FUNCTION", context)] <- "GO_MF"
  answer[is.na(answer) & grepl("GO_?CC|CELLULAR_COMPONENT", context)] <- "GO_CC"
  
  answer
}

infer_go_database_from_id <- function(go_ids) {
  go_ids <- clean_text(go_ids)
  
  answer <- rep(NA_character_, length(go_ids))
  
  go_positions <- which(
    !is.na(go_ids) &
      grepl("^GO:[0-9]{7}$", go_ids)
  )
  
  if (length(go_positions) == 0) {
    return(answer)
  }
  
  candidate_ids <- unique(go_ids[go_positions])
  
  valid_ids <- intersect(
    candidate_ids,
    AnnotationDbi::keys(OrgDb, keytype = "GOALL")
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
    !is.na(ontology_table$GOALL) &
      !is.na(ontology_table$ONTOLOGYALL),
    c("GOALL", "ONTOLOGYALL"),
    drop = FALSE
  ]
  
  ontology_table <- unique(ontology_table)
  
  ontology_by_id <- tapply(
    ontology_table$ONTOLOGYALL,
    ontology_table$GOALL,
    function(x) unique(x)[1]
  )
  
  matched_ontology <- unname(
    ontology_by_id[
      go_ids[go_positions]
    ]
  )
  
  answer[go_positions] <- ifelse(
    is.na(matched_ontology),
    NA_character_,
    paste0("GO_", matched_ontology)
  )
  
  answer
}

normalize_kegg_name <- function(x) {
  x <- toupper(clean_text(x))
  x <- gsub("SIGNALLING", "SIGNALING", x)
  x <- gsub("MUS MUSCULUS|MOUSE", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

standardize_gsea <- function(object, comparison, database) {
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
    setSize = suppressWarnings(as.numeric(x$setSize)),
    enrichmentScore = suppressWarnings(as.numeric(x$enrichmentScore)),
    NES = suppressWarnings(as.numeric(x$NES)),
    pvalue = suppressWarnings(as.numeric(x$pvalue)),
    p.adjust = suppressWarnings(as.numeric(x$p.adjust)),
    qvalues = if ("qvalues" %in% colnames(x)) {
      suppressWarnings(as.numeric(x$qvalues))
    } else {
      NA_real_
    },
    rank = if ("rank" %in% colnames(x)) {
      suppressWarnings(as.numeric(x$rank))
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

# -----------------------------
# 4. Read complete DE tables and build rankings
# -----------------------------

read_complete_de_table <- function(file_path, comparison_name) {
  if (!file.exists(file_path)) {
    stop("Cannot find:\n", file_path)
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
  
  colnames(x) <- make.unique(clean_text(colnames(x)), sep = "_")
  
  ensembl_column <- if ("EnsemblID_clean" %in% colnames(x)) {
    "EnsemblID_clean"
  } else if ("Geneid" %in% colnames(x)) {
    "Geneid"
  } else {
    stop(comparison_name, " has no EnsemblID_clean/Geneid column.")
  }
  
  x$EnsemblID <- sub(
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
    !is.na(x$EnsemblID) &
      grepl("^ENSMUSG[0-9]+$", x$EnsemblID),
    ,
    drop = FALSE
  ]
  
  if ("GeneName" %in% colnames(x)) {
    symbol_value <- clean_text(x$GeneName)
  } else {
    symbol_value <- rep(NA_character_, nrow(x))
  }
  
  missing_symbol <- is.na(symbol_value) | symbol_value == ""
  
  if (any(missing_symbol)) {
    symbol_map <- suppressMessages(
      AnnotationDbi::mapIds(
        OrgDb,
        keys = unique(x$EnsemblID[missing_symbol]),
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
    
    symbol_value[missing_symbol] <- unname(
      symbol_map[
        x$EnsemblID[missing_symbol]
      ]
    )
  }
  
  x$SymbolForGSEA <- toupper(clean_text(symbol_value))
  
  if ("stat" %in% colnames(x)) {
    x$RankScore <- suppressWarnings(as.numeric(x$stat))
    rank_method <- "DESeq2 Wald statistic"
  } else if (
    all(
      c("log2FoldChange", "pvalue") %in%
      colnames(x)
    )
  ) {
    lfc <- suppressWarnings(as.numeric(x$log2FoldChange))
    p_value <- suppressWarnings(as.numeric(x$pvalue))
    
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
    !is.na(x$SymbolForGSEA) &
      x$SymbolForGSEA != "" &
      is.finite(x$RankScore),
    ,
    drop = FALSE
  ]
  
  x <- x[
    order(-abs(x$RankScore)),
    ,
    drop = FALSE
  ]
  
  x <- x[
    !duplicated(x$SymbolForGSEA),
    ,
    drop = FALSE
  ]
  
  if ("Direction" %in% colnames(x)) {
    direction_upper <- toupper(clean_text(x$Direction))
    
    x$DirectionForCount <- "Not_significant"
    x$DirectionForCount[direction_upper == "UP"] <- "Up"
    x$DirectionForCount[direction_upper == "DOWN"] <- "Down"
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
        paste(missing_direction_columns, collapse = ", ")
      )
    }
    
    padj_value <- suppressWarnings(as.numeric(x$padj))
    shrunken_lfc <- suppressWarnings(
      as.numeric(x$log2FoldChange_shrunken)
    )
    
    x$DirectionForCount <- "Not_significant"
    
    x$DirectionForCount[
      !is.na(padj_value) &
        padj_value < deg_padj_cutoff &
        !is.na(shrunken_lfc) &
        shrunken_lfc >= deg_abs_shrunken_log2fc_cutoff
    ] <- "Up"
    
    x$DirectionForCount[
      !is.na(padj_value) &
        padj_value < deg_padj_cutoff &
        !is.na(shrunken_lfc) &
        shrunken_lfc <= -deg_abs_shrunken_log2fc_cutoff
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
    
    sort(rank_vector, decreasing = TRUE)
  }
)

up_deg_lists <- lapply(
  de_tables,
  function(x) {
    unique(
      x$SymbolForGSEA[
        x$DirectionForCount == "Up"
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
  file.path(output_dir, "00_ranked_gene_list_summary.tsv")
)

# -----------------------------
# 5. Read Sheet4 selected pathways
# -----------------------------

if (!file.exists(selected_pathway_excel)) {
  stop("Cannot find:\n", selected_pathway_excel)
}

selected_raw <- as.data.frame(
  read_excel(
    selected_pathway_excel,
    sheet = selected_sheet
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(selected_raw) <- make.unique(
  clean_text(colnames(selected_raw)),
  sep = "_"
)

if (
  !"ID" %in% colnames(selected_raw) ||
  !"Description" %in% colnames(selected_raw)
) {
  stop("Sheet4 must contain ID and Description columns.")
}

if ("Selected" %in% colnames(selected_raw)) {
  selected_value <- toupper(clean_text(selected_raw$Selected))
  
  recognized_yes <- selected_value %in%
    c(
      "YES",
      "Y",
      "TRUE",
      "1",
      "是",
      "保留"
    )
  
  if (any(recognized_yes, na.rm = TRUE)) {
    selected_raw <- selected_raw[
      recognized_yes,
      ,
      drop = FALSE
    ]
  }
}

selected_id <- clean_text(selected_raw$ID)

database_from_excel <- normalize_database(
  database = if ("Database" %in% colnames(selected_raw)) {
    selected_raw$Database
  } else {
    NA
  },
  source_file = if ("SourceFile" %in% colnames(selected_raw)) {
    selected_raw$SourceFile
  } else {
    NA
  },
  panel = if ("Panel" %in% colnames(selected_raw)) {
    selected_raw$Panel
  } else {
    NA
  },
  id = selected_id
)

database_from_annotation <- infer_go_database_from_id(
  selected_id
)

database <- database_from_excel

database[
  !is.na(database_from_annotation)
] <- database_from_annotation[
  !is.na(database_from_annotation)
]

if (any(is.na(database))) {
  stop(
    "Some Sheet4 pathway IDs could not be classified as GO/KEGG:\n",
    paste(
      unique(selected_id[is.na(database)]),
      collapse = ", "
    )
  )
}

selected_targets <- data.frame(
  Database = database,
  ID = selected_id,
  Description = clean_text(selected_raw$Description),
  PlotLabel = paste0(
    selected_id,
    ", ",
    clean_text(selected_raw$Description)
  ),
  PlotOrder = if (
    "PlotOrder" %in% colnames(selected_raw)
  ) {
    suppressWarnings(as.numeric(selected_raw$PlotOrder))
  } else {
    seq_len(nrow(selected_raw))
  },
  SourceExcelRow = seq_len(nrow(selected_raw)) + 1L,
  stringsAsFactors = FALSE
)

selected_targets$PlotOrder[
  !is.finite(selected_targets$PlotOrder)
] <- selected_targets$SourceExcelRow[
  !is.finite(selected_targets$PlotOrder)
]

selected_targets <- selected_targets[
  !is.na(selected_targets$ID) &
    !is.na(selected_targets$Description),
  ,
  drop = FALSE
]

selected_targets <- selected_targets[
  !duplicated(
    paste(
      selected_targets$Database,
      selected_targets$ID,
      sep = " | "
    )
  ),
  ,
  drop = FALSE
]

write_tsv(
  selected_targets,
  file.path(output_dir, "00_selected_pathway_targets.tsv")
)

# -----------------------------
# 6. Complete GO GSEA
# -----------------------------

complete_gsea_tables <- list()

for (comparison_name in comparison_order) {
  
  for (ontology in c("BP", "MF", "CC")) {
    
    database_name <- paste0("GO_", ontology)
    
    set.seed(gsea_seed)
    
    go_object <- tryCatch(
      clusterProfiler::gseGO(
        geneList = rank_lists[[comparison_name]],
        OrgDb = OrgDb,
        keyType = "SYMBOL",
        ont = ontology,
        exponent = 1,
        minGSSize = minimum_gene_set_size,
        maxGSSize = maximum_gene_set_size,
        eps = 0,
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
    
    go_table <- standardize_gsea(
      go_object,
      comparison_name,
      database_name
    )
    
    complete_gsea_tables[[
      paste(
        comparison_name,
        database_name,
        sep = " | "
      )
    ]] <- go_table
    
    write_tsv(
      go_table,
      file.path(
        complete_gsea_dir,
        paste0(
          gsub(" ", "_", comparison_name),
          "_",
          database_name,
          "_complete_GSEA.tsv"
        )
      )
    )
  }
}

# -----------------------------
# 7. Download and parse mouse KEGG
# -----------------------------

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
        "KEGG library download failed: ",
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
      "Could not download ",
      kegg_library_name,
      "."
    )
  }
}

kegg_lines <- readLines(
  cached_kegg_file,
  warn = FALSE,
  encoding = "UTF-8"
)

kegg_lines <- kegg_lines[
  nzchar(trimws(kegg_lines))
]

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
    "No valid KEGG_2019_Mouse pathway-to-gene mapping was parsed."
  )
}

kegg_term2gene <- unique(
  do.call(rbind, term2gene_rows)
)

kegg_term2name <- unique(
  do.call(rbind, term2name_rows)
)

write_tsv(
  kegg_term2gene,
  file.path(
    output_dir,
    "00_KEGG_2019_Mouse_TERM2GENE.tsv"
  )
)

# -----------------------------
# 8. Match Sheet4 KEGG terms
# -----------------------------

selected_kegg <- selected_targets[
  selected_targets$Database == "KEGG",
  ,
  drop = FALSE
]

selected_kegg$NormalizedDescription <- normalize_kegg_name(
  selected_kegg$Description
)

library_name_map <- data.frame(
  Term = kegg_term2name$Term,
  NormalizedTerm = normalize_kegg_name(
    kegg_term2name$Term
  ),
  NormalizedName = normalize_kegg_name(
    kegg_term2name$Name
  ),
  stringsAsFactors = FALSE
)

selected_kegg$GSEAID <- NA_character_
selected_kegg$KEGGMatchMethod <- NA_character_

for (
  row_index in seq_len(nrow(selected_kegg))
) {
  
  current_id <- as.character(
    selected_kegg$ID[row_index]
  )
  
  current_name <- selected_kegg$NormalizedDescription[row_index]
  
  matched_term <- NA_character_
  match_method <- NA_character_
  
  exact_hits <- library_name_map$Term[
    library_name_map$NormalizedName == current_name |
      library_name_map$NormalizedTerm == current_name
  ]
  
  exact_hits <- unique(
    exact_hits[
      !is.na(exact_hits)
    ]
  )
  
  if (length(exact_hits) > 0) {
    matched_term <- exact_hits[1]
    match_method <- "exact normalized pathway name"
  }
  
  if (is.na(matched_term)) {
    
    explicit_name <- switch(
      current_id,
      "mmu04657" = "IL-17 signaling pathway",
      "mmu04668" = "TNF signaling pathway",
      "mmu04064" = "NF-kappa B signaling pathway",
      "mmu04620" = "Toll-like receptor signaling pathway",
      NA_character_
    )
    
    if (!is.na(explicit_name)) {
      
      explicit_hits <- library_name_map$Term[
        library_name_map$NormalizedName ==
          normalize_kegg_name(explicit_name) |
          library_name_map$NormalizedTerm ==
          normalize_kegg_name(explicit_name)
      ]
      
      explicit_hits <- unique(
        explicit_hits[
          !is.na(explicit_hits)
        ]
      )
      
      if (length(explicit_hits) > 0) {
        matched_term <- explicit_hits[1]
        match_method <- paste0(
          current_id,
          " explicit pathway-name match"
        )
      }
    }
  }
  
  if (
    is.na(matched_term) &&
    !is.na(current_name) &&
    current_name != ""
  ) {
    
    candidate_names <- library_name_map$NormalizedName
    
    distance <- as.numeric(
      adist(
        current_name,
        candidate_names
      )
    )
    
    denominator <- pmax(
      nchar(current_name),
      nchar(candidate_names),
      1
    )
    
    normalized_distance <- distance / denominator
    
    best_index <- which.min(
      normalized_distance
    )
    
    if (
      length(best_index) == 1 &&
      is.finite(normalized_distance[best_index]) &&
      normalized_distance[best_index] <= 0.35
    ) {
      
      matched_term <- library_name_map$Term[
        best_index
      ]
      
      match_method <- paste0(
        "fuzzy normalized-name match; distance=",
        round(
          normalized_distance[best_index],
          3
        )
      )
    }
  }
  
  selected_kegg$GSEAID[row_index] <- matched_term
  selected_kegg$KEGGMatchMethod[row_index] <- match_method
}

selected_kegg$KEGGMappedGenes <- vapply(
  selected_kegg$GSEAID,
  function(term_name) {
    
    if (is.na(term_name)) {
      return(0L)
    }
    
    as.integer(
      length(
        unique(
          kegg_term2gene$Gene[
            kegg_term2gene$Term == term_name
          ]
        )
      )
    )
  },
  FUN.VALUE = integer(1)
)

write_tsv(
  selected_kegg[
    ,
    c(
      "ID",
      "Description",
      "GSEAID",
      "KEGGMatchMethod",
      "KEGGMappedGenes"
    ),
    drop = FALSE
  ],
  file.path(
    output_dir,
    "00_Sheet4_KEGG_mapping_diagnostic.tsv"
  )
)

# -----------------------------
# 9. Complete KEGG GSEA
# -----------------------------

for (comparison_name in comparison_order) {
  
  gsea_arguments <- list(
    geneList = rank_lists[[comparison_name]],
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
  
  set.seed(gsea_seed)
  
  kegg_object <- tryCatch(
    do.call(
      clusterProfiler::GSEA,
      gsea_arguments
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
  
  complete_gsea_tables[[
    paste(
      comparison_name,
      "KEGG",
      sep = " | "
    )
  ]] <- kegg_table
  
  write_tsv(
    kegg_table,
    file.path(
      complete_gsea_dir,
      paste0(
        gsub(" ", "_", comparison_name),
        "_KEGG_complete_GSEA.tsv"
      )
    )
  )
}

# -----------------------------
# 10. Combine complete GSEA
# -----------------------------

all_complete_gsea <- do.call(
  rbind,
  complete_gsea_tables
)

rownames(all_complete_gsea) <- NULL

if (nrow(all_complete_gsea) == 0) {
  stop("No complete GSEA results were generated.")
}

write_tsv(
  all_complete_gsea,
  file.path(
    complete_gsea_dir,
    "all_complete_GSEA_results.tsv"
  )
)

# -----------------------------
# 11. Pathway gene maps for DEG Count
# -----------------------------

selected_go_ids <- unique(
  selected_targets$ID[
    grepl(
      "^GO_",
      selected_targets$Database
    )
  ]
)

if (length(selected_go_ids) > 0) {
  
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
        "SYMBOL",
        "ONTOLOGYALL"
      )
    )
  )
  
  go_pathway_gene_map <- go_pathway_gene_map[
    !is.na(go_pathway_gene_map$GOALL) &
      !is.na(go_pathway_gene_map$SYMBOL) &
      !is.na(go_pathway_gene_map$ONTOLOGYALL),
    c(
      "GOALL",
      "SYMBOL",
      "ONTOLOGYALL"
    ),
    drop = FALSE
  ]
  
  go_pathway_gene_map$SYMBOL <- toupper(
    clean_text(go_pathway_gene_map$SYMBOL)
  )
  
  go_pathway_gene_map <- data.frame(
    Database = paste0(
      "GO_",
      go_pathway_gene_map$ONTOLOGYALL
    ),
    ID = go_pathway_gene_map$GOALL,
    Gene = go_pathway_gene_map$SYMBOL,
    stringsAsFactors = FALSE
  )
  
  go_pathway_gene_map <- unique(
    go_pathway_gene_map
  )
  
} else {
  
  go_pathway_gene_map <- data.frame(
    Database = character(0),
    ID = character(0),
    Gene = character(0),
    stringsAsFactors = FALSE
  )
}

kegg_pathway_gene_map <- data.frame(
  Database = "KEGG",
  ID = selected_kegg$ID,
  GSEAID = selected_kegg$GSEAID,
  stringsAsFactors = FALSE
)

kegg_gene_rows <- list()
kegg_gene_index <- 0L

for (row_index in seq_len(nrow(kegg_pathway_gene_map))) {
  
  current_gsea_id <- kegg_pathway_gene_map$GSEAID[row_index]
  
  if (is.na(current_gsea_id)) {
    next
  }
  
  current_genes <- unique(
    kegg_term2gene$Gene[
      kegg_term2gene$Term == current_gsea_id
    ]
  )
  
  if (length(current_genes) == 0) {
    next
  }
  
  kegg_gene_index <- kegg_gene_index + 1L
  
  kegg_gene_rows[[kegg_gene_index]] <- data.frame(
    Database = "KEGG",
    ID = kegg_pathway_gene_map$ID[row_index],
    Gene = current_genes,
    stringsAsFactors = FALSE
  )
}

if (length(kegg_gene_rows) > 0) {
  kegg_pathway_gene_map_expanded <- unique(
    do.call(
      rbind,
      kegg_gene_rows
    )
  )
} else {
  kegg_pathway_gene_map_expanded <- data.frame(
    Database = character(0),
    ID = character(0),
    Gene = character(0),
    stringsAsFactors = FALSE
  )
}

selected_pathway_gene_map <- rbind(
  go_pathway_gene_map,
  kegg_pathway_gene_map_expanded
)

write_tsv(
  selected_pathway_gene_map,
  file.path(
    output_dir,
    "00_selected_pathway_gene_map_for_DEGCount.tsv"
  )
)

calculate_deg_count <- function(
    comparison_name,
    database_name,
    pathway_id
) {
  
  pathway_genes <- unique(
    selected_pathway_gene_map$Gene[
      selected_pathway_gene_map$Database == database_name &
        selected_pathway_gene_map$ID == pathway_id
    ]
  )
  
  pathway_genes <- pathway_genes[
    !is.na(pathway_genes)
  ]
  
  if (length(pathway_genes) == 0) {
    return(0L)
  }
  
  current_up_degs <- up_deg_lists[[comparison_name]]
  
  as.integer(
    length(
      intersect(
        pathway_genes,
        current_up_degs
      )
    )
  )
}

# -----------------------------
# 12. Build selected-pathway plot table
# -----------------------------

selected_for_plot <- selected_targets

selected_for_plot$GSEAID <- selected_for_plot$ID

selected_for_plot$GSEAID[
  selected_for_plot$Database == "KEGG"
] <- selected_kegg$GSEAID[
  match(
    selected_for_plot$ID[
      selected_for_plot$Database == "KEGG"
    ],
    selected_kegg$ID
  )
]

plot_rows <- list()
plot_index <- 0L

for (target_index in seq_len(nrow(selected_for_plot))) {
  
  target <- selected_for_plot[
    target_index,
    ,
    drop = FALSE
  ]
  
  for (comparison_name in comparison_order) {
    
    plot_index <- plot_index + 1L
    
    match_rows <- all_complete_gsea[
      all_complete_gsea$Comparison == comparison_name &
        all_complete_gsea$Database == target$Database &
        all_complete_gsea$ID == target$GSEAID,
      ,
      drop = FALSE
    ]
    
    if (nrow(match_rows) > 1) {
      
      match_rows <- match_rows[
        order(
          match_rows$p.adjust,
          -abs(match_rows$NES)
        ),
        ,
        drop = FALSE
      ]
      
      match_rows <- match_rows[
        1,
        ,
        drop = FALSE
      ]
    }
    
    deg_count <- calculate_deg_count(
      comparison_name,
      target$Database,
      target$ID
    )
    
    if (nrow(match_rows) == 1) {
      
      plot_rows[[plot_index]] <- data.frame(
        Comparison = comparison_name,
        Database = target$Database,
        ID = target$ID,
        GSEAID = target$GSEAID,
        PlotLabel = target$PlotLabel,
        PlotOrder = target$PlotOrder,
        FoundInCompleteGSEA = TRUE,
        MissingReason = NA_character_,
        setSize = match_rows$setSize,
        enrichmentScore = match_rows$enrichmentScore,
        NES = match_rows$NES,
        pvalue = match_rows$pvalue,
        p.adjust = match_rows$p.adjust,
        qvalues = match_rows$qvalues,
        rank = match_rows$rank,
        leading_edge = match_rows$leading_edge,
        core_enrichment = match_rows$core_enrichment,
        DEGCount = deg_count,
        stringsAsFactors = FALSE
      )
      
    } else {
      
      plot_rows[[plot_index]] <- data.frame(
        Comparison = comparison_name,
        Database = target$Database,
        ID = target$ID,
        GSEAID = target$GSEAID,
        PlotLabel = target$PlotLabel,
        PlotOrder = target$PlotOrder,
        FoundInCompleteGSEA = FALSE,
        MissingReason = if (
          target$Database == "KEGG" &&
          is.na(target$GSEAID)
        ) {
          "Selected KEGG pathway could not be matched to KEGG_2019_Mouse"
        } else {
          "Pathway was not returned after annotation/size filtering"
        },
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
  plot_rows
)

plot_table$Comparison <- factor(
  plot_table$Comparison,
  levels = comparison_order
)

# -----------------------------
# 13. Put KEGG at TOP and GO below
# -----------------------------

term_info <- unique(
  plot_table[
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

go_levels <- rev(go_terms$PlotLabel)
kegg_levels <- rev(kegg_terms$PlotLabel)

y_levels <- c(
  go_levels,
  kegg_levels
)

plot_table$PlotLabel <- factor(
  plot_table$PlotLabel,
  levels = y_levels
)

# -----------------------------
# 14. Color and size scales
# -----------------------------

positive_adjusted_p <- plot_table$p.adjust[
  is.finite(plot_table$p.adjust) &
    plot_table$p.adjust > 0
]

smallest_positive <- if (length(positive_adjusted_p) > 0) {
  min(positive_adjusted_p)
} else {
  1e-300
}

plot_table$PAdjustForColour <- ifelse(
  !is.na(plot_table$p.adjust) &
    plot_table$p.adjust <= 0,
  smallest_positive / 10,
  plot_table$p.adjust
)

valid_colour_values <- plot_table$PAdjustForColour[
  is.finite(plot_table$PAdjustForColour) &
    plot_table$PAdjustForColour > 0
]

colour_limits <- if (length(valid_colour_values) == 0) {
  c(1e-300, 1)
} else {
  range(valid_colour_values)
}

if (colour_limits[1] == colour_limits[2]) {
  colour_limits <- c(
    colour_limits[1] / 10,
    min(
      1,
      colour_limits[2] * 10
    )
  )
}

maximum_deg_count <- max(
  plot_table$DEGCount,
  na.rm = TRUE
)

if (
  !is.finite(maximum_deg_count) ||
  maximum_deg_count <= 0
) {
  maximum_deg_count <- 1
}

# -----------------------------
# 15. Final plot
# -----------------------------

p_final <- ggplot(
  plot_table,
  aes(
    x = Comparison,
    y = PlotLabel
  )
) +
  geom_point(
    data = plot_table[
      plot_table$FoundInCompleteGSEA,
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
    data = plot_table[
      !plot_table$FoundInCompleteGSEA,
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
    range = c(3, 10),
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
    title = "Sheet4-selected GO/KEGG GSEA",
    subtitle = paste0(
      "GO = mouse OrgDb; KEGG = Enrichr KEGG_2019_Mouse; ",
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
      unique(plot_table$PlotLabel)
    ) +
    2.5
)

print(p_final)

# -----------------------------
# 16. Save
# -----------------------------

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

ggsave(
  filename = file.path(
    plot_dir,
    "01_Sheet4_GO_KEGG_GSEA_KEGG_TOP_GO_BOTTOM.png"
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
    plot_dir,
    "01_Sheet4_GO_KEGG_GSEA_KEGG_TOP_GO_BOTTOM.pdf"
  ),
  plot = p_final,
  width = 12,
  height = plot_height,
  units = "in",
  device = grDevices::cairo_pdf,
  bg = "white"
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "R_sessionInfo.txt"
  )
)

cat(
  "\nBMDC GSEA finished.\n"
)

cat(
  "Output directory:\n",
  output_dir,
  "\n",
  sep = ""
)
