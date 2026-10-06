# ============================================================================
# Human PBMC paired RNA-seq: UP/DOWN selected GO/KEGG module dot plots
#
# X axis  : Ctrl, IMP, Y180
# Y axis  : manually selected GO/KEGG terms (UP and DOWN plotted separately)
# Color   : mean gene-wise VST z-score (NO subtraction of Ctrl)
# Size    : percent of selected module genes detected in the group
#
# IMPORTANT
# 1) This is a bulk-RNA-seq adaptation of the reference figure.
# 2) Bulk RNA-seq cannot calculate "% of cells expressing a module".
#    Circle size here is "% of selected module genes detected".
# 3) The GeneName column from an ORA table contains overlapping selected DEGs,
#    not necessarily every gene annotated to the complete GO/KEGG pathway.
# 4) final 10 selected.xlsx contains the complete 22-column selected-term table:
#      Sheet2 = UP pathways; Sheet1 = DOWN pathways.
# 5) Differential-expression/pathway selection remains paired. The plot itself
#    is descriptive: it shows VST expression patterns for Ctrl, IMP and Y180.
# ============================================================================

# ----------------------------------------------------------------------------
# 1. Modify these paths and parameters
# ----------------------------------------------------------------------------

count_file <- "D:/个人记录/科研/Y180/BMDC.reads.count.txt"
module_excel <- paste0(
  "D:/个人记录/科研/Y180/",
  "BMDC Paired_GO_KEGG_results/final 10 selected BMDC.xlsx"
)

# Excel worksheet assignment requested by the user.
sheet_specs <- data.frame(
  Sheet = c("Sheet4"),
  DirectionPanel = c("UP"),
  stringsAsFactors = FALSE
)

output_dir <- file.path(
  dirname(module_excel),
  "Selected_module_dotplot_UP_DOWN"
)

# Column names in the raw-count table
gene_id_column <- "Geneid"
gene_name_column <- "GeneName"

# Nine paired samples
sample_order <- c(
  "Ctrl_1", "Ctrl_2", "Ctrl_3",
  "IMP_1", "IMP_2", "IMP_3",
  "Y180_1", "Y180_2", "Y180_3"
)

group_order <- c("Ctrl", "IMP", "Y180")

# Low-count filtering before VST
minimum_count <- 10
minimum_samples <- 3

# A module gene is considered detected in a group when its normalized count is
# at least detection_count in at least detection_min_donors of the 3 donors.
detection_count <- 50
detection_min_donors <- 2

# Minimum matched genes required to draw a pathway/module
minimum_module_genes <- 2

# Show GO/KEGG ID before the pathway label.
show_id_in_label <- TRUE

# Fixed symmetric color limit. Use NULL for automatic limits.
# Example: fixed_color_limit <- 1.5
fixed_color_limit <- NULL

# ----------------------------------------------------------------------------
# 2. Check packages
# ----------------------------------------------------------------------------

required_packages <- c("DESeq2", "readxl", "ggplot2", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R packages: ", paste(missing_packages, collapse = ", "), "\n",
      "Install CRAN packages with install.packages(), and DESeq2 with ",
      "BiocManager::install(\"DESeq2\")."
    )
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(count_file)) {
  stop(
    paste0(
      "Raw-count file was not found:\n",
      count_file,
      "\nUse forward slashes in the Windows path."
    )
  )
}

if (!file.exists(module_excel)) {
  stop(
    paste0(
      "Selected-pathway Excel file was not found:\n",
      module_excel,
      "\nUse forward slashes in the Windows path."
    )
  )
}

# ----------------------------------------------------------------------------
# 3. Helper functions
# ----------------------------------------------------------------------------

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u3000", "", x, fixed = TRUE)
  trimws(x)
}

read_table_by_extension <- function(file_path) {
  extension <- tolower(tools::file_ext(file_path))
  
  if (extension %in% c("xlsx", "xls")) {
    return(
      as.data.frame(
        readxl::read_excel(
          path = file_path,
          col_names = TRUE,
          .name_repair = "minimal"
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }
  
  as.data.frame(
    read.delim(
      file = file_path,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

split_gene_symbols <- function(x) {
  x <- clean_text(x)
  
  if (is.na(x) || x == "") {
    return(character(0))
  }
  
  genes <- unlist(
    strsplit(x, split = "[/;,|]", perl = TRUE),
    use.names = FALSE
  )
  
  genes <- toupper(clean_text(genes))
  unique(
    genes[
      !is.na(genes) &
        genes != "" &
        !genes %in% c("NA", "N/A", "NULL")
    ]
  )
}

safe_filename <- function(x, index) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  
  if (is.na(x) || x == "") {
    x <- paste0("panel_", index)
  }
  
  paste0(sprintf("%02d", index), "_", x)
}

# ----------------------------------------------------------------------------
# 4. Read and validate raw counts
# ----------------------------------------------------------------------------

expression_data <- read_table_by_extension(count_file)

required_count_columns <- c(
  gene_id_column,
  gene_name_column,
  sample_order
)

missing_count_columns <- setdiff(
  required_count_columns,
  colnames(expression_data)
)

if (length(missing_count_columns) > 0) {
  stop(
    paste0(
      "Raw-count table is missing columns:\n",
      paste(missing_count_columns, collapse = ", ")
    )
  )
}

count_frame <- expression_data[, sample_order, drop = FALSE]

for (sample_name in sample_order) {
  count_frame[[sample_name]] <- suppressWarnings(
    as.numeric(count_frame[[sample_name]])
  )
}

if (anyNA(count_frame)) {
  bad_position <- which(is.na(as.matrix(count_frame)), arr.ind = TRUE)[1, ]
  
  stop(
    paste0(
      "Raw-count matrix contains a non-numeric or missing value at row ",
      bad_position[1], ", sample ",
      colnames(count_frame)[bad_position[2]], "."
    )
  )
}

count_matrix_raw <- as.matrix(count_frame)
storage.mode(count_matrix_raw) <- "numeric"

if (any(count_matrix_raw < 0)) {
  stop("Raw-count matrix contains negative values.")
}

non_integer <- abs(count_matrix_raw - round(count_matrix_raw)) > 1e-8

if (any(non_integer)) {
  bad_position <- which(non_integer, arr.ind = TRUE)[1, ]
  
  stop(
    paste0(
      "Counts must be integers. Example: row ", bad_position[1],
      ", sample ", colnames(count_matrix_raw)[bad_position[2]],
      " = ", count_matrix_raw[bad_position[1], bad_position[2]], "."
    )
  )
}

gene_id_clean <- sub(
  "\\..*$",
  "",
  clean_text(expression_data[[gene_id_column]])
)

valid_gene_id <- !is.na(gene_id_clean) & gene_id_clean != ""
gene_id_clean <- gene_id_clean[valid_gene_id]
count_matrix_raw <- count_matrix_raw[valid_gene_id, , drop = FALSE]

# Sum duplicated Geneid rows while preserving Geneid for upstream processing.
count_matrix <- rowsum(
  count_matrix_raw,
  group = gene_id_clean,
  reorder = FALSE
)
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

# One GeneName annotation per cleaned Geneid. Prefer a non-empty symbol.
annotation_source <- data.frame(
  Geneid = gene_id_clean,
  GeneName = clean_text(expression_data[[gene_name_column]][valid_gene_id]),
  stringsAsFactors = FALSE
)

annotation_source$NameMissing <- (
  is.na(annotation_source$GeneName) |
    annotation_source$GeneName == ""
)

annotation_source <- annotation_source[
  order(annotation_source$Geneid, annotation_source$NameMissing),
  ,
  drop = FALSE
]

annotation <- annotation_source[
  !duplicated(annotation_source$Geneid),
  c("Geneid", "GeneName"),
  drop = FALSE
]

# Low-expression filter
keep_gene <- rowSums(count_matrix >= minimum_count) >= minimum_samples
count_matrix_filtered <- count_matrix[keep_gene, , drop = FALSE]

if (nrow(count_matrix_filtered) < 2) {
  stop("Fewer than two genes remain after low-expression filtering.")
}

message(
  "Genes retained for VST: ", nrow(count_matrix_filtered),
  " / ", nrow(count_matrix)
)

# ----------------------------------------------------------------------------
# 5. Paired sample information and VST
# ----------------------------------------------------------------------------

sample_info <- data.frame(
  Group = factor(
    c(rep("Ctrl", 3), rep("IMP", 3), rep("Y180", 3)),
    levels = group_order
  ),
  Donor = factor(
    c(
      "Donor1", "Donor2", "Donor3",
      "Donor1", "Donor2", "Donor3",
      "Donor1", "Donor2", "Donor3"
    ),
    levels = c("Donor1", "Donor2", "Donor3")
  ),
  row.names = sample_order,
  stringsAsFactors = FALSE
)

if (!identical(colnames(count_matrix_filtered), rownames(sample_info))) {
  stop("Count-matrix columns and sample-information rows are not in the same order.")
}

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = sample_info,
  design = ~ Donor + Group
)

dds <- DESeq2::DESeq(dds, quiet = TRUE)
vsd <- DESeq2::vst(dds, blind = FALSE)

vst_matrix <- SummarizedExperiment::assay(vsd)
normalized_counts <- DESeq2::counts(dds, normalized = TRUE)

# ----------------------------------------------------------------------------
# 6. Convert Geneid to GeneName only after VST
#    If several Geneids map to one GeneName, keep the most highly expressed ID.
# ----------------------------------------------------------------------------

symbol_map <- annotation$GeneName[
  match(rownames(vst_matrix), annotation$Geneid)
]

symbol_key <- toupper(clean_text(symbol_map))
mean_normalized_count <- rowMeans(normalized_counts)

mapping_table <- data.frame(
  Geneid = rownames(vst_matrix),
  GeneName = symbol_map,
  GeneKey = symbol_key,
  MeanNormalizedCount = mean_normalized_count,
  stringsAsFactors = FALSE
)

valid_symbol <- !is.na(mapping_table$GeneKey) & mapping_table$GeneKey != ""
mapping_table <- mapping_table[valid_symbol, , drop = FALSE]

mapping_table <- mapping_table[
  order(
    mapping_table$GeneKey,
    -mapping_table$MeanNormalizedCount
  ),
  ,
  drop = FALSE
]

representative_gene <- !duplicated(mapping_table$GeneKey)
mapping_table <- mapping_table[representative_gene, , drop = FALSE]

representative_position <- match(
  mapping_table$Geneid,
  rownames(vst_matrix)
)

vst_by_symbol <- vst_matrix[representative_position, , drop = FALSE]
normalized_by_symbol <- normalized_counts[
  representative_position,
  ,
  drop = FALSE
]

rownames(vst_by_symbol) <- mapping_table$GeneKey
rownames(normalized_by_symbol) <- mapping_table$GeneKey

# Gene-wise z-scores across all nine samples
gene_mean <- rowMeans(vst_by_symbol)
gene_sd <- apply(vst_by_symbol, 1, stats::sd)
variable_gene <- is.finite(gene_sd) & gene_sd > 0

vst_by_symbol <- vst_by_symbol[variable_gene, , drop = FALSE]
normalized_by_symbol <- normalized_by_symbol[variable_gene, , drop = FALSE]
gene_mean <- gene_mean[variable_gene]
gene_sd <- gene_sd[variable_gene]

z_matrix <- sweep(vst_by_symbol, 1, gene_mean, FUN = "-")
z_matrix <- sweep(z_matrix, 1, gene_sd, FUN = "/")

# IMPORTANT: do not subtract matched Ctrl here.
# The paired design (~ Donor + Group) is retained for DESeq2 normalization and
# dispersion modeling, while this descriptive plot uses each sample's own
# gene-wise VST z-score. Therefore Ctrl, IMP and Y180 all have visible values.

# ----------------------------------------------------------------------------
# 7. Read and validate the selected GO/KEGG Excel table
# ----------------------------------------------------------------------------

expected_module_columns <- c(
  "Panel", "Comparison", "Direction", "Database", "SourceFile",
  "PlotOrder", "Selected", "DisplayLabel", "Theme", "ID",
  "Description", "GeneRatio", "BgRatio", "RichFactor",
  "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
  "geneID", "Count", "GeneName"
)

available_sheets <- readxl::excel_sheets(module_excel)
missing_sheets <- setdiff(sheet_specs$Sheet, available_sheets)

if (length(missing_sheets) > 0) {
  stop(
    paste0(
      "The following worksheet(s) were not found in final 10 selected.xlsx:\n",
      paste(missing_sheets, collapse = ", "),
      "\nAvailable worksheets: ",
      paste(available_sheets, collapse = ", ")
    )
  )
}

read_selected_module_sheet <- function(sheet_name, direction_panel) {
  current_table <- as.data.frame(
    readxl::read_xlsx(
      path = module_excel,
      sheet = sheet_name,
      col_names = TRUE,
      .name_repair = "minimal"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(current_table) <- clean_text(colnames(current_table))
  
  missing_columns <- setdiff(
    expected_module_columns,
    colnames(current_table)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Worksheet ", sheet_name,
        " has incorrect or merged headers. Missing:\n",
        paste(missing_columns, collapse = ", "),
        "\n\nRequired headers:\n",
        paste(expected_module_columns, collapse = " | ")
      )
    )
  }
  
  # Ignore completely blank extra Excel columns after GeneName.
  current_table <- current_table[
    ,
    expected_module_columns,
    drop = FALSE
  ]
  
  current_table$ExcelRow <- seq_len(nrow(current_table)) + 1L
  description_text <- clean_text(current_table$Description)
  current_table <- current_table[
    !is.na(description_text) & description_text != "",
    ,
    drop = FALSE
  ]
  
  if (nrow(current_table) == 0) {
    stop("No non-empty pathway rows were found in worksheet ", sheet_name, ".")
  }
  
  # If the Selected column contains recognizable Yes values, keep only those.
  # If it is entirely blank/unrecognized, treat all non-empty rows as selected.
  selected_text <- tolower(clean_text(current_table$Selected))
  selected_yes <- selected_text %in% c(
    "yes", "y", "true", "1", "是", "保留"
  )
  
  if (any(selected_yes)) {
    current_table <- current_table[selected_yes, , drop = FALSE]
  } else {
    message(
      "Worksheet ", sheet_name,
      ": no recognizable Selected=Yes value; all non-empty rows are used."
    )
  }
  
  if (nrow(current_table) == 0) {
    stop("No selected pathway rows remain in worksheet ", sheet_name, ".")
  }
  
  gene_name_text <- clean_text(current_table$GeneName)
  
  if (all(is.na(gene_name_text) | gene_name_text == "")) {
    stop(
      "Every GeneName cell is empty in worksheet ", sheet_name,
      ". Module scoring requires gene symbols."
    )
  }
  
  empty_gene_rows <- which(is.na(gene_name_text) | gene_name_text == "")
  
  if (length(empty_gene_rows) > 0) {
    warning(
      "Worksheet ", sheet_name,
      " has empty GeneName at Excel row(s): ",
      paste(current_table$ExcelRow[empty_gene_rows], collapse = ", ")
    )
  }
  
  plot_order <- suppressWarnings(as.numeric(current_table$PlotOrder))
  missing_plot_order <- !is.finite(plot_order)
  
  if (any(missing_plot_order)) {
    plot_order[missing_plot_order] <- seq_len(nrow(current_table))[
      missing_plot_order
    ]
    warning(
      "Worksheet ", sheet_name,
      ": missing/non-numeric PlotOrder values were replaced by Excel row order."
    )
  }
  
  display_label <- clean_text(current_table$DisplayLabel)
  missing_display <- is.na(display_label) | display_label == ""
  display_label[missing_display] <- clean_text(
    current_table$Description[missing_display]
  )
  
  if (show_id_in_label) {
    id_text <- clean_text(current_table$ID)
    has_id <- !is.na(id_text) & id_text != ""
    display_label[has_id] <- paste0(
      id_text[has_id], ", ", display_label[has_id]
    )
  }
  
  current_table$ExcelSheet <- sheet_name
  current_table$ExcelDirection <- clean_text(current_table$Direction)
  current_table$Direction <- direction_panel
  current_table$PanelClean <- direction_panel
  current_table$PlotOrderNumeric <- plot_order
  current_table$PlotLabel <- display_label
  
  current_table
}

module_table_list <- lapply(
  seq_len(nrow(sheet_specs)),
  function(sheet_index) {
    read_selected_module_sheet(
      sheet_name = sheet_specs$Sheet[sheet_index],
      direction_panel = sheet_specs$DirectionPanel[sheet_index]
    )
  }
)

module_table <- do.call(rbind, module_table_list)
rownames(module_table) <- NULL

# ----------------------------------------------------------------------------
# 8. Calculate module scores and detection percentages
# ----------------------------------------------------------------------------

sample_level_rows <- list()
group_level_rows <- list()
coverage_rows <- list()
sample_row_index <- 1L
group_row_index <- 1L
coverage_row_index <- 1L

for (module_index in seq_len(nrow(module_table))) {
  requested_genes <- split_gene_symbols(
    module_table$GeneName[module_index]
  )
  
  matched_genes <- intersect(
    requested_genes,
    rownames(z_matrix)
  )
  
  missing_genes <- setdiff(
    requested_genes,
    rownames(z_matrix)
  )
  
  coverage_rows[[coverage_row_index]] <- data.frame(
    ModuleRow = module_index,
    Panel = module_table$PanelClean[module_index],
    ExcelSheet = module_table$ExcelSheet[module_index],
    ExcelRow = module_table$ExcelRow[module_index],
    PlotOrder = module_table$PlotOrderNumeric[module_index],
    ID = module_table$ID[module_index],
    Description = module_table$Description[module_index],
    PlotLabel = module_table$PlotLabel[module_index],
    Database = module_table$Database[module_index],
    GeneRatio = module_table$GeneRatio[module_index],
    BgRatio = module_table$BgRatio[module_index],
    RichFactor = module_table$RichFactor[module_index],
    FoldEnrichment = module_table$FoldEnrichment[module_index],
    EnrichmentZScore = module_table$zScore[module_index],
    EnrichmentPValue = module_table$pvalue[module_index],
    EnrichmentAdjustedP = module_table$p.adjust[module_index],
    EnrichmentQValue = module_table$qvalue[module_index],
    EnrichmentCount = module_table$Count[module_index],
    RequestedGeneNumber = length(requested_genes),
    MatchedGeneNumber = length(matched_genes),
    CoveragePercent = if (length(requested_genes) > 0) {
      100 * length(matched_genes) / length(requested_genes)
    } else {
      0
    },
    MatchedGenes = paste(matched_genes, collapse = "/"),
    MissingGenes = paste(missing_genes, collapse = "/"),
    stringsAsFactors = FALSE
  )
  coverage_row_index <- coverage_row_index + 1L
  
  if (length(matched_genes) < minimum_module_genes) {
    warning(
      "Skipped module '", module_table$PlotLabel[module_index],
      "': only ", length(matched_genes), " genes matched."
    )
    next
  }
  
  # One descriptive module-expression score per sample.
  # This uses VST gene-wise z-scores directly; Ctrl is NOT subtracted.
  module_sample_score <- colMeans(
    z_matrix[matched_genes, , drop = FALSE],
    na.rm = TRUE
  )
  
  # Also export the uncentered mean VST value for auditing. It is not used for
  # the red-white-blue plot because genes have different baseline abundances.
  module_sample_mean_vst <- colMeans(
    vst_by_symbol[matched_genes, , drop = FALSE],
    na.rm = TRUE
  )
  
  for (sample_name in sample_order) {
    sample_level_rows[[sample_row_index]] <- data.frame(
      ModuleRow = module_index,
      Panel = module_table$PanelClean[module_index],
      ExcelSheet = module_table$ExcelSheet[module_index],
      ExcelRow = module_table$ExcelRow[module_index],
      PlotOrder = module_table$PlotOrderNumeric[module_index],
      ID = module_table$ID[module_index],
      Description = module_table$Description[module_index],
      PlotLabel = module_table$PlotLabel[module_index],
      Sample = sample_name,
      Donor = as.character(sample_info[sample_name, "Donor"]),
      Group = as.character(sample_info[sample_name, "Group"]),
      ModuleExpressionZ = module_sample_score[sample_name],
      MeanVSTExpression = module_sample_mean_vst[sample_name],
      MatchedGeneNumber = length(matched_genes),
      stringsAsFactors = FALSE
    )
    sample_row_index <- sample_row_index + 1L
  }
  
  # Group means and percent of selected genes detected
  for (current_group in group_order) {
    group_samples <- rownames(sample_info)[
      sample_info$Group == current_group
    ]
    
    group_scores <- module_sample_score[group_samples]
    group_mean_vst <- module_sample_mean_vst[group_samples]
    
    detected_gene <- rowSums(
      normalized_by_symbol[
        matched_genes,
        group_samples,
        drop = FALSE
      ] >= detection_count
    ) >= detection_min_donors
    
    group_level_rows[[group_row_index]] <- data.frame(
      ModuleRow = module_index,
      Panel = module_table$PanelClean[module_index],
      ExcelSheet = module_table$ExcelSheet[module_index],
      ExcelRow = module_table$ExcelRow[module_index],
      PlotOrder = module_table$PlotOrderNumeric[module_index],
      ID = module_table$ID[module_index],
      Description = module_table$Description[module_index],
      PlotLabel = module_table$PlotLabel[module_index],
      Comparison = module_table$Comparison[module_index],
      Direction = module_table$Direction[module_index],
      Database = module_table$Database[module_index],
      Theme = module_table$Theme[module_index],
      GeneRatio = module_table$GeneRatio[module_index],
      BgRatio = module_table$BgRatio[module_index],
      RichFactor = module_table$RichFactor[module_index],
      FoldEnrichment = module_table$FoldEnrichment[module_index],
      EnrichmentZScore = module_table$zScore[module_index],
      EnrichmentPValue = module_table$pvalue[module_index],
      EnrichmentAdjustedP = module_table$p.adjust[module_index],
      EnrichmentQValue = module_table$qvalue[module_index],
      EnrichmentCount = module_table$Count[module_index],
      Group = current_group,
      ModuleExpressionZ = mean(group_scores, na.rm = TRUE),
      ModuleExpressionZSE = if (length(group_scores) > 1) {
        stats::sd(group_scores, na.rm = TRUE) / sqrt(length(group_scores))
      } else {
        NA_real_
      },
      MeanVSTExpression = mean(group_mean_vst, na.rm = TRUE),
      MeanVSTExpressionSE = if (length(group_mean_vst) > 1) {
        stats::sd(group_mean_vst, na.rm = TRUE) /
          sqrt(length(group_mean_vst))
      } else {
        NA_real_
      },
      PercentModuleGenesDetected = 100 * mean(detected_gene),
      MatchedGeneNumber = length(matched_genes),
      stringsAsFactors = FALSE
    )
    group_row_index <- group_row_index + 1L
  }
}

coverage_table <- do.call(rbind, coverage_rows)

write.table(
  coverage_table,
  file = file.path(output_dir, "module_gene_matching_report.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (length(group_level_rows) == 0) {
  stop(
    "No module had enough matched genes. Check GeneName symbols and the ",
    "module_gene_matching_report.tsv file."
  )
}

sample_score_table <- do.call(rbind, sample_level_rows)
group_score_table <- do.call(rbind, group_level_rows)

write.table(
  sample_score_table,
  file = file.path(output_dir, "module_scores_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  group_score_table,
  file = file.path(output_dir, "module_dotplot_group_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 9. Draw one reference-style dot plot per Panel
# ----------------------------------------------------------------------------

panel_levels <- unique(group_score_table$Panel)

# Use one common symmetric color scale for both UP and DOWN plots so their
# colors are directly comparable.
if (is.null(fixed_color_limit)) {
  global_color_limit <- max(
    abs(group_score_table$ModuleExpressionZ),
    na.rm = TRUE
  )
  
  if (!is.finite(global_color_limit) || global_color_limit == 0) {
    global_color_limit <- 1
  }
  
  global_color_limit <- ceiling(global_color_limit * 2) / 2
} else {
  global_color_limit <- fixed_color_limit
}

for (panel_index in seq_along(panel_levels)) {
  current_panel <- panel_levels[panel_index]
  plot_data <- group_score_table[
    group_score_table$Panel == current_panel,
    ,
    drop = FALSE
  ]
  
  pathway_order_table <- unique(
    plot_data[, c("ModuleRow", "PlotOrder", "PlotLabel"), drop = FALSE]
  )
  pathway_order_table <- pathway_order_table[
    order(pathway_order_table$PlotOrder, pathway_order_table$ModuleRow),
    ,
    drop = FALSE
  ]
  
  # ModuleRow is an internal unique key. This prevents two selected rows with
  # the same visible label from being accidentally merged on the y axis.
  pathway_levels <- paste0("Module_", pathway_order_table$ModuleRow)
  wrapped_pathway_labels <- vapply(
    pathway_order_table$PlotLabel,
    function(current_label) {
      paste(
        strwrap(current_label, width = 48),
        collapse = "\n"
      )
    },
    character(1)
  )
  label_dictionary <- setNames(
    wrapped_pathway_labels,
    pathway_levels
  )
  
  plot_data$PathwayKey <- paste0("Module_", plot_data$ModuleRow)
  plot_data$PathwayKey <- factor(
    plot_data$PathwayKey,
    levels = rev(pathway_levels)
  )
  plot_data$Group <- factor(plot_data$Group, levels = group_order)
  
  module_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = Group,
      y = PathwayKey,
      size = PercentModuleGenesDetected,
      fill = ModuleExpressionZ
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "grey40",
      stroke = 0.45,
      alpha = 0.75
    ) +
    ggplot2::scale_y_discrete(labels = label_dictionary) +
    ggplot2::scale_fill_gradient2(
      name = "Mean gene-wise\nVST z-score",
      low = "#1749E6",
      mid = "#FFFFFF",
      high = "#FF2A0A",
      midpoint = 0,
      limits = c(-global_color_limit, global_color_limit),
      oob = scales::squish
    ) +
    ggplot2::scale_size_continuous(
      name = "% selected module\ngenes detected",
      limits = c(0, 100),
      range = c(1.5, 10),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    ggplot2::labs(
      title = paste(current_panel, "selected GO/KEGG pathways"),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5,
        size = 17
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(
        colour = "black",
        size = 11
      ),
      axis.line = ggplot2::element_line(colour = "black"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(10, 25, 10, 10)
    )
  
  plot_height <- max(5, 0.65 * length(pathway_levels) + 2.5)
  plot_width <- max(
    10,
    min(
      15,
      9 + max(nchar(pathway_order_table$PlotLabel), na.rm = TRUE) / 30
    )
  )
  output_stem <- safe_filename(current_panel, panel_index)
  
  # Export the exact numerical data used for this panel. These files make it
  # possible to audit every circle back to the three donor-level samples.
  panel_sample_data <- sample_score_table[
    sample_score_table$Panel == current_panel,
    ,
    drop = FALSE
  ]
  
  write.table(
    panel_sample_data,
    file = file.path(
      output_dir,
      paste0(output_stem, "_module_scores_by_sample.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    plot_data,
    file = file.path(
      output_dir,
      paste0(output_stem, "_module_dotplot_group_summary.tsv")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(output_stem, "_module_dotplot.png")
    ),
    plot = module_plot,
    width = plot_width,
    height = plot_height,
    dpi = 600,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(output_stem, "_module_dotplot.pdf")
    ),
    plot = module_plot,
    width = plot_width,
    height = plot_height,
    device = "pdf",
    bg = "white"
  )
}

message(
  "Finished. Generated ", length(panel_levels), " panel plot(s).\n",
  "Output directory: ",
  normalizePath(output_dir, winslash = "/", mustWork = FALSE)
)