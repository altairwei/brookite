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
#' @param filter Filter name used to query BioMart
#' @param dest_attrs BioMart attributes to pull down
#' @param collapse_sep Separator for collapsed attributes
#' @return a tibble
#'
#' @export
pull_annotation <- function(
  genes, mart_info,
  filter = "ensembl_gene_id",
  dest_attrs = default_biomart_attributes,
  collapse_sep = "|"
) {
  dataset_conn <- biomaRt::useMart(
    mart_info$mart,
    dataset = mart_info$dataset,
    host = mart_info$host
  )

  anno <- biomaRt::getBM(
      attributes = c(filter, dest_attrs),
      filters = filter,
      values = genes,
      mart = dataset_conn,
      quote = "\""
  )

  anno_collapse <- anno %>%
    tibble::as_tibble() %>%
    dplyr::group_by(.data[[filter]]) %>%
    dplyr::group_modify(function(x, y) {
      df <- apply(x, 2, FUN = collapse, collapse_sep, simplify = FALSE) %>%
        tibble::as_tibble()
      df
    })

  names(anno_collapse) <- ifelse(
    names(anno_collapse) %in% dest_attrs,
    names(dest_attrs)[match(names(anno_collapse), dest_attrs)],
    names(anno_collapse)
  )

  anno_collapse
}

.isSingleString <- AnnotationForge:::.isSingleString
.isSingleStringOrNull <- AnnotationForge:::.isSingleStringOrNull
.isSingleStringOrNA <- AnnotationForge:::.isSingleStringOrNA
isSingleNumber <- S4Vectors::isSingleNumber

#' Copy from AnnotationForge
p_makeOrgPackage <- function(
  data, version, maintainer, author,
  outputDir = getwd(), tax_id, genus = NULL,
  species = NULL, goTable = NA,
  databaseOnly = FALSE, verbose = TRUE,
  gid_type = "eg"
) {
  ## unique names
  if (length(unique(names(data))) != length(data))
      stop("All elements of '...' must have unique names")

  ## expired names
  blackListedNames <- c("genes", "metadata")
  if (any(names(data) %in% blackListedNames))
      stop("'genes' and 'metadata' are reserved. Please choose different ",
          "names for elements of '...'.")

  ## coerce to data.frame
  data <- lapply(data, as.data.frame)

  ## drop rownames, no duplicated rows
  data <- lapply(data, function(x) {
              rownames(x) <- NULL
              if (any(duplicated(x)))
                  stop("data.frames in '...' cannot contain duplicated rows")
              x
          })

  ## unique colnames for each data.frame
  .noGID <- function(x) x[!(x %in% "GID")]
  colnamesUni <- unique(.noGID(unlist(sapply(data, colnames))))
  colnamesAll <- .noGID(unlist(sapply(data, colnames)))
  names(colnamesAll) <- NULL
  if (any(colnamesUni != colnamesAll))
      stop(paste0("data.frames should have completely unique names for all ",
                  "fields that are not the primary gene id 'GID'"))

  ## first column of each data.frame must be gene ID (GID)
  colnameGIDs <- sapply(data, function(x) colnames(x)[1])
  if (any(colnameGIDs != "GID"))
      stop("The 1st column must always be the gene ID 'GID'")

  ## check GID type
  GIDCols <- unique(sapply(data,
      function(x) class(x[["GID"]])
  ))
  if (length(GIDCols) > 1)
      stop(paste0("The type of data in the 'GID' columns must be the same ",
                  "for all data.frames"))

  ## check other arguments
  if (!.isSingleString(version))
      stop("'version' must be a single string")
  if (!.isSingleString(maintainer))
      stop("'maintainer' must be a single string")
  if (!.isSingleString(author))
      stop("'author' must be a single string")
  if (outputDir != "." && file.access(outputDir)[[1]] != 0)
      stop("Selected outputDir '", outputDir, "' does not exist.")
  if (!(isSingleNumber(tax_id) || .isSingleString(tax_id)))
      stop("'tax_id' must be a single integer")
  if (!is.integer(tax_id))
      tax_id <- as.integer(tax_id)
  if (!.isSingleStringOrNull(genus))
      stop("'genus' must be a single string or NULL")
  if (!.isSingleStringOrNull(species))
      stop("'species' must be a single string or NULL")
  ## only an NA internally - a NULL is what would have come in from outside...
  if (!.isSingleStringOrNA(goTable))
      stop("'goTable' argument needs to be a single string or NULL")
  if (!is.na(goTable) && !(goTable %in% names(data)))
      stop("'goTable' must be a name from the data.frames passed in '...'")

  ## genus and species
  if (is.null(genus))
      genus <- GenomeInfoDb:::lookup_organism_by_tax_id(tax_id)[["genus"]]
  if (is.null(species)) {
      species <- GenomeInfoDb:::lookup_organism_by_tax_id(tax_id)[["species"]]
      species <- gsub(" ", ".", species)
  }

  dbName <- paste0(
    "org.",
    paste0(toupper(substr(genus, 1, 1)), species),
    ".", gid_type
  )

  dbFileName <- file.path(outputDir, paste0(dbName, ".sqlite"))

  AnnotationForge:::makeOrgDbFromDataFrames(
    data, tax_id, genus, species, dbFileName, goTable)

  if (databaseOnly) {
      ## return the path to the database file
      file.path(outputDir, dbFileName)
  } else {
      seed <- new("AnnDbPkgSeed",
                  Package = paste0(dbName, ".db"),
                  Version = version,
                  Author = author,
                  Maintainer = maintainer,
                  PkgTemplate = "NOSCHEMA.DB",
                  AnnObjPrefix = dbName,
                  organism = paste(genus, species),
                  species = paste(genus, species),
                  biocViews = "annotation",
                  manufacturerUrl = "no manufacturer",
                  manufacturer = "no manufacturer",
                  chipName = "no manufacturer")

      AnnotationForge::makeAnnDbPkg(seed, dbFileName, dest_dir = outputDir)

      ## cleanup
      message("Now deleting temporary database file")
      file.remove(dbFileName)
      ## return the path to the dir that was just created.
      file.path(outputDir, paste0(dbName, ".db"))
  }
}

#' Make a OrgDb package from annotations available on a BioMart database
#'
#' @inheritParams AnnotationForge::makeOrgPackage
#' @param gid_type the type of central identifier.
#' @param mart which BioMart database to use. Get the list of all available
#' BioMart databases with the \code{\link[biomaRt]{listMarts}} function from
#' the biomaRt package.
#' @param dataset which dataset from BioMart.
#' @param host The host URL of the BioMart.
#' @export
#' @examples
#' makeOrgDbFromBiomart(
#'   mart = "plants_mart",
#'   dataset = "taestivum_eg_gene",
#'   host = "https://plants.ensembl.org",
#'   version = "0.0.0.9000",
#'   maintainer = "Altair Wei <altair_wei@outlook.com>",
#'   author = "Altair Wei <altair_wei@outlook.com>",
#'   outputDir = ".",
#'   tax_id = "4565",
#'   genus = "Triticum",
#'   species = "aestivum",
#'   gid_type = "iwgsc"
#' )
makeOrgDbFromBiomart <- function(
  tax_id,
  maintainer,
  author,
  genus = NULL,
  species = NULL,
  version = "0.0.0.9000",
  gid_type = "eg",
  mart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "www.ensembl.org",
  outputDir = getwd()
) {
  dataset_conn <- biomaRt::useMart(
    mart,
    dataset = dataset,
    host = host
  )

  anno <- biomaRt::getBM(
      attributes = c(
        "ensembl_gene_id", # GID
        "description", # DESCRIPTION
        "chromosome_name", # CHROMOSOME
        "external_gene_name", # GENENAME
        "go_id", # GO
        "go_linkage_type" # EVIDENCE
      ),
      filters = "chromosome_name",
      values = biomaRt::keys(dataset_conn, keytype = "chromosome_name"),
      mart = dataset_conn,
      quote = "\""
  )

  anno <- dplyr::rename(anno,
    GID = ensembl_gene_id,
    GENENAME = external_gene_name,
    DESCRIPTION = description,
    CHROMOSOME = chromosome_name,
    GO = go_id,
    EVIDENCE = go_linkage_type
  )

  # Gene data.frame： GID, GENENAME and DESCRIPTION
  gene_df <- anno %>%
    dplyr::group_by(GID) %>%
    dplyr::summarise(
      GENENAME = GENENAME[[1]],
      DESCRIPTION = DESCRIPTION[[1]]
    ) %>%
    as.data.frame()

  # Chromosome data.frame: GID, CHROMOSOME
  chrom_df <- anno %>%
    dplyr::group_by(GID) %>%
    dplyr::summarise(
      CHROMOSOME = CHROMOSOME[[1]]
    ) %>%
    as.data.frame()

  # GO data.frame: GID, GO, EVIDENCE
  go_df <- anno %>%
    dplyr::group_by(GID, GO) %>%
    dplyr::summarise(
      EVIDENCE = EVIDENCE[[1]]
    ) %>%
    dplyr::filter(GO != "", EVIDENCE != "") %>%
    as.data.frame()

  db_file <- p_makeOrgPackage(
    list(
      gene_info = gene_df,
      chromosome = chrom_df,
      go = go_df
    ),
    goTable = "go",
    tax_id = tax_id,
    genus = genus,
    species = species,
    version = version,
    author = author,
    maintainer = maintainer,
    gid_type = gid_type
  )
}