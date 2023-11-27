
#' Tissue
#'
#' @description Return tissue types stored in HCNetlas.
#'
#' @return a vector of characters, with tissue names.
#' @export
#'
#' @examples
#' Tissue()
Tissue <- function(){
  return(Tissues)
}

#' CellType
#'
#' @description Return cell types stored in HCNetlas.
#'
#' @return a vector of characters, with cell-type names.
#' @export
#'
#' @examples
#' CellType()
CellType <- function(){
  return(CellTypes)
}

#' getCelltype
#'
#' @description Return cell types of networks stored in HCNetlas in a specific tissue.
#'
#' @param tissue a character of tissue names. Should be included in output of Tissues() function.
#'
#' @return a vector of characters, with cell-type names.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' getCelltype(tissue="LNG")
getCelltype <- function(tissue = NULL){
  tryCatch(
    if (tissue %in% Tissues){
      items <- str_split(names(HCNetlas_list),'_',n=2)

      tissue_v <- unlist(lapply(items, function(x) x[1]))
      celltype_v <- unlist(lapply(items, function(x) x[2]))

      accessible_celltypes <- celltype_v[tissue_v == tissue]

      return(accessible_celltypes)
    }
    else {
      cat("Invalid tissue.\nTo see valid tissue list, please type Tissue()\n")

    },
    error = function(e) cat("Invalid tissue.\nTo see valid tissue list, please type Tissue()\n")
  )
}



#' MergeCGN.HCNetlas
#'
#' @description Merge disease CGNs and HCNetlas CGNs with corresponding tissue.
#'
#' @param disease_cgn scHumanNet output CGN list. output from SortAddLLS() function in scHumanNet package.
#' @param disease a character specifying a disease.
#' @param tissue a tissue name for selecting HCNetlas CGNs. We recomand this corresponds with tissue where disease scRNA-seq data were obtained. Should be included in output of Tissues() function.
#'
#' @return a list of CGN networks. Disease CGNs and HCNetlas healhty CGNs merged.
#' @import dplyr
#'
#' @export
#'
#' @examples
#' MergeCGN.HCNetlas(disease_cgn=sorted.net.list, disease="SLE", tissue="LNG")
MergeCGN.HCNetlas <- function(disease_cgn  = NULL,disease = NULL, tissue = NULL){
  tryCatch(
    if (tissue %in% Tissues){
      items <- str_split(names(HCNetlas_list),'_',n=2)

      tissue_v <- unlist(lapply(items, function(x) x[1]))
      celltype_v <- unlist(lapply(items, function(x) x[2]))

      control_cgn <- HCNetlas_list[tissue_v == tissue]

      names(control_cgn) <- gsub(tissue,"HCNetlas",names(control_cgn))
      names(disease_cgn) <- paste0(disease,"_",names(disease_cgn))
      merged_cgn <- c(control_cgn, disease_cgn)

      return(merged_cgn)
    }
    else {
      cat("Invalid tissue.\nTo see valid tissue list, please type Tissue()\n")
    },
    error = function(e) cat("Invalid tissue.\nTo see valid tissue list, please type Tissue()\n")
  )
}


#' MergeCGN
#'
#' @description Merge disease CGNs and healthy CGNs to generate list of networks.
#'
#' @param disease_cgn a list of disease CGNs. scHumanNet output from SortAddLLS() function in scHumanNet package.
#' @param disease a character specifying a disease.
#' @param control_cgn a list of healthy CGNs. scHumanNet output from SortAddLLS() function in scHumanNet package.
#' @param control a character specifying control.
#'
#' @return a list of CGN networks. Disease CGNs and healhty CGNs merged.
#' @export
#'
#' @examples
#' MergeCGN(disease_cgn=dis.sorted.net.list, disease="SLE", control_cgn=ctl.sorted.net.list, control="Control")
MergeCGN <- function(disease_cgn  = NULL,disease = NULL,control_cgn=NULL,control=NULL){

  names(control_cgn) <- paste0(control,"_",names(control_cgn))
  names(disease_cgn) <- paste0(disease,"_",names(disease_cgn))
  merged_cgn <- c(control_cgn, disease_cgn)

  return(merged_cgn)
}
