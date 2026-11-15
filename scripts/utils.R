# =============================================================================
# utils.R -- shared helpers for the Seurat scRNA-seq workflow
#
# Sourced by every numbered script. Run scripts from the repository root:
#   Rscript scripts/01_qc.R --config config.yaml
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
})

get_config_path <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- grep("^--config", args)
  if (length(hit) > 0) {
    a <- args[hit[1]]
    if (grepl("=", a)) return(sub("^--config=", "", a))
    if (length(args) > hit[1]) return(args[hit[1] + 1])
  }
  env <- Sys.getenv("PIPELINE_CONFIG", unset = "")
  if (nzchar(env)) return(env)
  "config.yaml"
}

load_config <- function(path = get_config_path()) {
  if (!file.exists(path)) {
    stop("Configuration file not found: ", path,
         "\nRun the scripts from the repository root, or pass --config <file>.")
  }
  cfg <- yaml::read_yaml(path)
  cfg$config_file <- normalizePath(path)
  res <- cfg$paths$results
  cfg$dirs <- list(results = res,
                   figures = file.path(res, "figures"),
                   tables = file.path(res, "tables"),
                   rds = file.path(res, "rds"),
                   logs = file.path(res, "logs"))
  for (d in cfg$dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (is.null(cfg$seed)) cfg$seed <- 42
  set.seed(cfg$seed)
  cfg
}

p_fig <- function(cfg, f) file.path(cfg$dirs$figures, f)
p_table <- function(cfg, f) file.path(cfg$dirs$tables, f)
p_rds <- function(cfg, f) file.path(cfg$dirs$rds, f)

require_files <- function(...) {
  f <- c(...)
  missing <- f[!file.exists(f)]
  if (length(missing) > 0) {
    stop("Required input(s) not found:\n  ", paste(missing, collapse = "\n  "),
         "\nRun the preceding pipeline step(s) first.")
  }
  invisible(TRUE)
}

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

log_step <- function(title) {
  cat("\n", strrep("=", 70), "\n", title, "\n", strrep("=", 70), "\n", sep = "")
}

write_session_info <- function(cfg, script_name) {
  out <- file.path(cfg$dirs$logs, paste0(script_name, "_sessionInfo.txt"))
  writeLines(c(paste("Script:", script_name),
               paste("Config:", cfg$config_file),
               paste("Seed:", cfg$seed),
               paste("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
               "", utils::capture.output(utils::sessionInfo())), out)
  log_msg("sessionInfo written to ", out)
}

#' Save a ggplot / patchwork object as PNG.
save_plot <- function(cfg, plot, name, width, height, dpi = 300) {
  force(plot)
  f <- p_fig(cfg, paste0(name, ".png"))
  ggplot2::ggsave(f, plot, width = width, height = height, dpi = dpi, bg = "white")
  log_msg("Saved ", basename(f))
  invisible(f)
}

#' Save a figure drawn with base / grid graphics as PNG.
save_base_plot <- function(cfg, draw, name, width, height, dpi = 300) {
  f <- p_fig(cfg, paste0(name, ".png"))
  grDevices::png(f, width = width, height = height, units = "in", res = dpi)
  tryCatch(draw(), finally = grDevices::dev.off())
  log_msg("Saved ", basename(f))
  invisible(f)
}

#' Load a SingleR reference by name: a celldex dataset or "custom"
#' (annotation$custom_reference_rds).
load_reference <- function(cfg, name) {
  if (identical(name, "custom")) {
    f <- cfg$annotation$custom_reference_rds
    if (is.null(f)) stop("annotation$custom_reference_rds is not set")
    return(readRDS(f))
  }
  fun <- switch(name,
                HumanPrimaryCellAtlas = "HumanPrimaryCellAtlasData",
                BlueprintEncode = "BlueprintEncodeData",
                MonacoImmune = "MonacoImmuneData",
                DatabaseImmuneCellExpression = "DatabaseImmuneCellExpressionData",
                NovershternHematopoietic = "NovershternHematopoieticData",
                stop("Unknown celldex reference: ", name))
  if (!requireNamespace("celldex", quietly = TRUE)) stop("Package 'celldex' is required.")
  getExportedValue("celldex", fun)()
}

#' Names of the references configured for annotation (celldex names, plus
#' "custom" when a custom reference file is given).
reference_names <- function(cfg) {
  refs <- unlist(cfg$annotation$references)
  if (!is.null(cfg$annotation$custom_reference_rds)) refs <- c(refs, "custom")
  if (length(refs) == 0) stop("No annotation references configured.")
  refs
}
