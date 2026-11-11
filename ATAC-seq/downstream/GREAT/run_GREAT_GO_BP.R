#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 180)

# GO Biological Process enrichment for differential accessibility regions.
# Usage: Rscript run_GREAT_GO_BP.R <DAR_results.csv> <output_dir> [FDR] [abs_log2FC]
# Input: CSV with region IDs in column 1 and chrom, start, end,
#        log2FoldChange, padj columns. Coordinates are 0-based, half-open.
# Positive log2FoldChange = higher accessibility in condition B vs A.

required_packages <- c(
    "rGREAT", "GenomicRanges", "IRanges", "S4Vectors", "GO.db",
    "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene"
)
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
    stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
    library(rGREAT)
    library(GenomicRanges)
    library(IRanges)
    library(S4Vectors)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 4L) {
    stop("Usage: Rscript run_GREAT_GO_BP.R <DAR_results.csv> <output_dir> [FDR] [abs_log2FC]")
}
input_file <- args[[1]]
output_dir <- args[[2]]
fdr_cutoff <- if (length(args) >= 3L) as.numeric(args[[3]]) else 0.05
log2fc_cutoff <- if (length(args) >= 4L) as.numeric(args[[4]]) else 0.5
analysis_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))
if (!is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff > 1 ||
    !is.finite(log2fc_cutoff) || log2fc_cutoff < 0 ||
    is.na(analysis_cores) || analysis_cores < 1) {
    stop("Invalid cutoff or core count.")
}

tss_source <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
standard_chromosomes <- paste0("chr", c(1:19, "X", "Y"))

log_message <- function(...) {
    message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
            paste0(..., collapse = ""))
}
write_tsv <- function(x, path) {
    write.table(x, file = path, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE, na = "NA")
}
to_granges <- function(x) {
    regions <- GRanges(
        seqnames = x$chrom,
        ranges = IRanges(start = as.integer(x$start) + 1L,
                         end = as.integer(x$end)),
        strand = "*"
    )
    names(regions) <- x$region_id
    regions
}
flatten_column <- function(x) {
    if (is.list(x) || inherits(x, "List")) {
        return(vapply(as.list(x), function(value) {
            paste(as.character(value), collapse = ";")
        }, character(1)))
    }
    as.character(x)
}
associations_to_dataframe <- function(associations) {
    output <- data.frame(
        region_id = if (is.null(names(associations))) {
            paste0("region_", seq_along(associations))
        } else {
            names(associations)
        },
        chromosome = as.character(seqnames(associations)),
        start_1based = start(associations),
        end = end(associations)
    )
    annotation <- mcols(associations)
    for (column in colnames(annotation)) {
        output[[column]] <- flatten_column(annotation[[column]])
    }
    output
}

if (!file.exists(input_file)) stop("Missing input file: ", input_file)
dat <- read.csv(input_file, row.names = 1, check.names = FALSE)
required_columns <- c("chrom", "start", "end", "log2FoldChange", "padj")
if (!all(required_columns %in% colnames(dat))) {
    stop("Missing columns: ", paste(setdiff(required_columns, colnames(dat)), collapse = ", "))
}
if (anyNA(dat[, c("chrom", "start", "end")]) ||
    any(!is.finite(dat$start)) || any(!is.finite(dat$end)) ||
    any(dat$start != floor(dat$start)) || any(dat$end != floor(dat$end)) ||
    any(dat$start < 0) || any(dat$end <= dat$start)) {
    stop("Invalid genomic coordinates.")
}
dat$region_id <- rownames(dat)
keep <- dat$chrom %in% standard_chromosomes
log_message("Removing ", sum(!keep), " regions on nonstandard chromosomes.")
dat <- dat[keep, , drop = FALSE]
if (!nrow(dat)) stop("No regions remain after chromosome filtering.")
if (anyDuplicated(dat[, c("chrom", "start", "end")])) {
    stop("Duplicate genomic intervals found in input.")
}

background <- to_granges(dat)
is_dar <- is.finite(dat$padj) & is.finite(dat$log2FoldChange) &
    dat$padj < fdr_cutoff & abs(dat$log2FoldChange) >= log2fc_cutoff
foreground_tables <- list(
    Lower_in_B = dat[is_dar & dat$log2FoldChange < 0, , drop = FALSE],
    Higher_in_B = dat[is_dar & dat$log2FoldChange > 0, , drop = FALSE]
)
log_message("Background: ", length(background), " tested regions.")
for (label in names(foreground_tables)) {
    log_message(label, ": ", nrow(foreground_tables[[label]]), " DARs.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_one_analysis <- function(foreground, label) {
    regions <- to_granges(foreground)
    log_message("Starting ", label, " — GO:BP.")
    started <- Sys.time()
    great_result <- rGREAT::great(
        gr = regions,
        gene_sets = "GO:BP",
        tss_source = tss_source,
        min_gene_set_size = 5,
        mode = "basalPlusExt",
        extend_from = "TSS",
        basal_upstream = 5000,
        basal_downstream = 1000,
        extension = 1000000,
        background = background,
        exclude = NULL,
        cores = analysis_cores,
        verbose = TRUE
    )
    enrichment <- rGREAT::getEnrichmentTable(great_result, min_region_hits = 1)
    write_tsv(enrichment, file.path(output_dir, paste0(label, "_GO_BP_all_terms.tsv")))
    associations <- rGREAT::getRegionGeneAssociations(great_result, use_symbols = TRUE)
    association_table <- associations_to_dataframe(associations)
    write_tsv(association_table,
              file.path(output_dir, paste0(label, "_region_gene_associations.tsv")))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
    log_message("Completed ", label, " in ", round(elapsed, 1), " minutes.")
    data.frame(
        direction = label,
        foreground_regions = length(regions),
        background_regions = length(background),
        genome = "mm10",
        gene_sets = "GO:BP",
        tss_source = tss_source,
        association_rule = "basalPlusExt",
        basal_upstream_bp = 5000,
        basal_downstream_bp = 1000,
        maximum_extension_bp = 1000000,
        adjusted_p_cutoff = fdr_cutoff,
        absolute_log2FC_cutoff = log2fc_cutoff,
        returned_terms = nrow(enrichment),
        elapsed_minutes = round(elapsed, 2),
        status = "completed",
        error_message = ""
    )
}

run_records <- list()
for (label in names(foreground_tables)) {
    foreground <- foreground_tables[[label]]
    if (!nrow(foreground)) {
        log_message("Skipping ", label, ": no DARs.")
        next
    }
    run_records[[label]] <- tryCatch(
        run_one_analysis(foreground, label),
        error = function(error) {
            log_message("ERROR: ", label, ": ", conditionMessage(error))
            data.frame(
                direction = label,
                foreground_regions = nrow(foreground),
                background_regions = length(background),
                genome = "mm10",
                gene_sets = "GO:BP",
                tss_source = tss_source,
                association_rule = "basalPlusExt",
                basal_upstream_bp = 5000,
                basal_downstream_bp = 1000,
                maximum_extension_bp = 1000000,
                adjusted_p_cutoff = fdr_cutoff,
                absolute_log2FC_cutoff = log2fc_cutoff,
                returned_terms = NA_integer_,
                elapsed_minutes = NA_real_,
                status = "failed",
                error_message = conditionMessage(error)
            )
        }
    )
    invisible(gc())
}
if (!length(run_records)) stop("No DAR sets were available for analysis.")
run_summary <- do.call(rbind, run_records)
rownames(run_summary) <- NULL
write_tsv(run_summary, file.path(output_dir, "GO_BP_run_summary.tsv"))
cat("\nRUN SUMMARY\n")
print(run_summary)
cat("\nSESSION INFORMATION\n")
print(sessionInfo())
if (any(run_summary$status == "failed")) {
    stop("One or more analyses failed. See GO_BP_run_summary.tsv.")
}
log_message("GO:BP analyses completed. Tables: ", output_dir)
