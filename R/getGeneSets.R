#' Get a collection of gene sets to perform enrichment on
#'
#' This function allows users to select libraries and specific 
#' gene.sets to form a GeneSetCollection that is a list of gene sets.
#
#' @param species The scientific name of the species of interest in 
#' order to get correct gene nomenclature
#' @param library Individual collection(s) of gene sets, e.g. c("H", "C5").
#' See \href{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp}{msigdbr}for
#' all MSigDB collections.
#' @param subcategory MSigDB sub-collection abbreviation, such as CGP or BP.
#' @param gene.sets Select gene sets or pathways, using specific names, 
#' example: pathways = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"). Will only be
#' honored if library is set, too.
#'
#' @examples 
#' GS <- getGeneSets(library = "H")
#' 
#' @export
#' 
#' @importFrom GSEABase GeneSet GeneSetCollection
#' @importFrom msigdbr msigdbr msigdbr_species
#' @importFrom stringr str_replace_all
#' 
#' @author Nick Borcherding, Jared Andrews
#' @return A list of gene sets from msigdbr.
getGeneSets <- function(species = "Homo sapiens", 
                        library = NULL, 
                        subcategory = NULL,
                        gene.sets = NULL) {
  spec <- msigdbr_species()$species_name
  if (!(species %in% spec)) {
    stop(paste0("Please select a compatible species: ", 
                paste(spec, collapse = ", ")))
  }
  
  if (!is.null(library)) {
    m_df_list <- list()
    for (x in seq_along(library)) {
      if (is.null(subcategory)) {
        tmp2 <- msigdbr(species = species, category = library[x])
      } else {
        tmp2 <- msigdbr(species = species, category = library[x], subcategory = subcategory)
      }
      m_df_list[[x]] <- tmp2
    }
    m_df <- dplyr::bind_rows(m_df_list)
    if (!is.null(gene.sets)) {
      m_df <- m_df[m_df$gs_name %in% gene.sets, ]
    }    
  } else {
    m_df <- msigdbr(species = species)
  }
  
  if (nrow(m_df) == 0) {
    warning("No gene sets found for the specified parameters.")
    return(NULL)
  }
  
  gs <- unique(m_df$gs_name)
  ls <- vector("list", length(gs))
  
  for (i in seq_along(gs)) {
    tmp <- unique(m_df$gene_symbol[m_df$gs_name == gs[i]])
    ls[[i]] <- GSEABase::GeneSet(tmp, setName = paste(gs[i]))
  }
  
  gsc <- GSEABase::GeneSetCollection(ls)
  mod.names <- stringr::str_replace_all(names(gsc), "_", "-")
  gsc <- GSEABase::geneIds(gsc)
  names(gsc) <- mod.names
  
  return(gsc)
}
