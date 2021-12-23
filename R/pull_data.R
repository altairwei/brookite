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
    Description = "description",
    Gene_ID = "ensembl_gene_id",
    GO_ID = "go_id",
    GO_Name = "name_1006",
    GO_Level = "namespace_1003"
)

# Collapse go_id and go_name
collapse <- function(x, sep = "|") {
    if (all(x == "")) {
        return("-")
    } else {
        return(paste(unique(x), collapse = sep))
    }
}

#' Pull Annotation from BioMart
#'
#' @param genes Ensembl gene id
#' @param mart_info A list contains BioMart information:
#' \code{host}, \code{mart} and \code{dataset}
#' @param dest_attrs BioMart attributes to pull down
#' @param collapse_sep Separator for collapsed attributes
#' @return a tibble
#'
#' @export
pull_annotation <- function(
  genes, mart_info,
  dest_attrs = default_biomart_attributes,
  collapse_sep = "|"
) {
  dataset_conn <- biomaRt::useMart(
    mart_info$mart,
    dataset = mart_info$dataset,
    host = mart_info$host
  )

  anno <- biomaRt::getBM(
      attributes = dest_attrs,
      filters = "ensembl_gene_id",
      values = genes,
      mart = dataset_conn,
      quote = "\""
  )

  anno_collapse <- anno %>%
    tibble::as_tibble() %>%
    dplyr::group_by(ensembl_gene_id) %>%
    dplyr::group_modify(function(x, y) {
      df <- apply(x, 2, FUN = collapse, collapse_sep, simplify = FALSE) %>%
        tibble::as_tibble()
      names(df) <- names(default_biomart_attributes)
      df
    })

  anno_collapse
}