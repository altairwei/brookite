#' Pull data table from biomart
#'
#' @param host Host url such as plants.ensembl.org
#' @param mart Mart name such as plants_mart
#' @param dataset Dataset name such as taestivum_eg_gene
#' @param dest_attrs Attributes you want to retrieve.
#'  It can be a named vector or named list.
#' @param filename File Name to save table.
#'
#' @export
pull_data_from_biomart <- function(
  host, mart, dataset,
  dest_attrs = default_biomart_attributes,
  filename = NULL) {

  # Build connections
  dataset_conn <- biomaRt::useMart(mart, dataset = dataset, host = host)

  # Get all annotated gene list from chrosomes
  all_gene_list <- biomaRt::getBM(
    attributes = c("ensembl_gene_id"),
    filters = "chromosome_name",
    values = biomaRt::keys(dataset_conn, keytype = "chromosome_name"),
    mart = dataset_conn
  ) %>% dplyr::as_tibble()
  cat(paste("Total genes: ",
    length(unique(all_gene_list[["ensembl_gene_id"]])), "\n"))

  # Get all data
  all_gene_data <- biomaRt::getBM(
      attributes = dest_attrs,
      filters = "ensembl_gene_id",
      values = all_gene_list[["ensembl_gene_id"]],
      mart = dataset_conn
  ) %>% dplyr::as_tibble()

  # Rename column if dest_attrs is a named vector.
  if (!is.null(names(dest_attrs))) {
    all_gene_data <- all_gene_data %>%
      dplyr::rename(!!!dest_attrs)
  }

  if (!is.null(filename)) {
    write.table(
      all_gene_data,
      file = filename,
      sep = "\t",
      row.name = FALSE,
      col.names = TRUE,
      quote = FALSE)
  } else {
    return(all_gene_data)
  }
}

#' Default biomart attributes to retrieve
#' @export
default_biomart_attributes <- c(
    Gene_ID = "ensembl_gene_id",
    GO_ID = "go_id",
    GO_Name = "name_1006",
    GO_Level = "namespace_1003"
)