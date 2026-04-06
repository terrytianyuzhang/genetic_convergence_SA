work_directory_yao_2023 <- file.path(work_directory, "yao_2023")
go_dir <- file.path(work_directory_yao_2023, "data", "final_data", "module_GO")
report_dir <- file.path(work_directory_yao_2023, "report")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

selected_terms <- list(
  `6` = c("GO:0006954", "GO:0001816", "GO:0002250"),
  `12` = c("GO:0042110", "GO:0046649", "GO:0050851", "GO:0050852"),
  `14` = c("GO:0001568", "GO:0001501", "GO:0042692"),
  `31` = c("GO:0034340", "GO:0051607", "GO:0140888", "GO:0001816"),
  `11` = c("GO:0019752", "GO:0009117", "GO:0044282"),
  `20` = c("GO:0050900", "GO:0006935", "GO:2000406"),
  `5` = c("GO:0015031", "GO:0006886", "GO:0072657"),
  `13` = c("GO:0007155", "GO:0072678", "GO:0007229"),
  `17` = c("GO:0098771", "GO:0055074", "GO:0030001")
)

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

build_module_row <- function(module_id, go_ids, p_adjust_cutoff = 0.05) {
  file_path <- file.path(go_dir, paste0("module_", module_id, ".csv"))
  if (!file.exists(file_path)) {
    stop(glue("GO result file not found for module {module_id}: {file_path}"))
  }

  go_dt <- fread(file_path)
  go_dt <- go_dt[ONTOLOGY == "BP"]
  selected_dt <- go_dt[ID %in% go_ids]

  missing_ids <- setdiff(go_ids, selected_dt$ID)
  if (length(missing_ids) > 0) {
    stop(glue("Module {module_id} is missing selected GO IDs: {paste(missing_ids, collapse = ', ')}"))
  }

  nonsignificant_ids <- selected_dt[p.adjust > p_adjust_cutoff, ID]
  if (length(nonsignificant_ids) > 0) {
    stop(glue("Module {module_id} has selected GO IDs with p.adjust > {p_adjust_cutoff}: {paste(nonsignificant_ids, collapse = ', ')}"))
  }

  selected_dt[, order_id := match(ID, go_ids)]
  setorder(selected_dt, order_id)

  selected_terms_text <- paste0(selected_dt$Description, " (", sub("^GO:", "", selected_dt$ID), ")")
  data.table(
    module_index = module_id,
    representative_terms = paste(selected_terms_text, collapse = "; ")
  )
}

table_dt <- rbindlist(
  lapply(names(selected_terms), function(module_id) {
    build_module_row(module_id, selected_terms[[module_id]])
  })
)
table_dt[, module_index := as.integer(module_index)]
setorder(table_dt, module_index)

csv_path <- file.path(report_dir, "33_selected_module_GO_terms.csv")
fwrite(table_dt, csv_path)

latex_lines <- c(
  "\\begin{table}[!tbp]",
  "\\caption{Representative GO terms associated with selected modules in Yao et al. (2023). For each module, representative significant biological process terms are reported with GO IDs.}",
  "\\vspace{3mm}",
  "\\centering",
  "\\setlength{\\tabcolsep}{12pt}",
  "\\renewcommand{\\arraystretch}{1.3}",
  "\\label{tab:GO_of_modules_Yao_2023_selected}",
  "\\begin{tabular}{>{\\raggedright\\arraybackslash}p{0.15\\linewidth} >{\\raggedright\\arraybackslash}p{0.75\\linewidth}}",
  "\\toprule",
  "\\textbf{Module Index} & \\textbf{Representative Gene Ontology Terms} \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(table_dt))) {
  latex_lines <- c(
    latex_lines,
    glue("{table_dt$module_index[i]} & {latex_escape(table_dt$representative_terms[i])} \\\\")
  )
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

tex_path <- file.path(report_dir, "33_selected_module_GO_terms.tex")
writeLines(latex_lines, tex_path)

message(glue("Saved CSV table to: {csv_path}"))
message(glue("Saved LaTeX table to: {tex_path}"))
