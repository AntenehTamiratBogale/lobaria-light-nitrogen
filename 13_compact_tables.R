# Keep Supplementary Tables S1–S6 in six worksheets.
# Sections A, B and C stay within their respective table worksheet.
cat("\nCombine supplementary table sections into six worksheets\n")
table_path <- file.path(OUTDIR, "Supplementary_Tables_S1-S6.xlsx")
table_names <- getSheetNames(table_path)
table_workbook <- createWorkbook()
for (table_number in 1:6) {
  name <- paste0("S", table_number)
  sections <- table_names[grepl(paste0("^", name, "[A-Z_]"), table_names)]
  addWorksheet(table_workbook, name)
  writeData(table_workbook, name, paste("Supplementary Table", name), startRow = 1)
  addStyle(table_workbook, name, createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  next_row <- 3
  for (section in sections) {
    values <- read.xlsx(table_path, sheet = section, startRow = 3, check.names = FALSE)
    names(values) <- sub("Functional.subset", "Protein.group", names(values), fixed = TRUE)
    for (column in names(values)) {
      if (is.character(values[[column]])) {
        values[[column]] <- gsub("Photobiont", "Phycobiont", values[[column]], fixed = TRUE)
        if (section == "S3A_Keywords" && !grepl("keyword|search", column, ignore.case = TRUE))
          values[[column]] <- gsub("final subset", "final functional group", values[[column]], fixed = TRUE)
      }
    }
    writeData(table_workbook, name, section, startRow = next_row)
    writeData(table_workbook, name, values, startRow = next_row + 1)
    addStyle(table_workbook, name, createStyle(textDecoration = "bold"), rows = next_row, cols = 1)
    addStyle(table_workbook, name, createStyle(textDecoration = "bold", wrapText = TRUE),
             rows = next_row + 1, cols = seq_len(ncol(values)), gridExpand = TRUE)
    numeric_columns <- which(vapply(values, is.numeric, logical(1)))
    if (length(numeric_columns) && nrow(values))
      addStyle(table_workbook, name, createStyle(numFmt = "General"),
               rows = next_row + 1 + seq_len(nrow(values)), cols = numeric_columns, gridExpand = TRUE)
    next_row <- next_row + nrow(values) + 4
  }
  freezePane(table_workbook, name, firstActiveRow = 5)
  setColWidths(table_workbook, name, cols = 1:30, widths = 24)
}
saveWorkbook(table_workbook, table_path, overwrite = TRUE)
cat("Saved Supplementary Tables S1–S6 with six worksheets.\n")
