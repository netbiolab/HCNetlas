#' Compactness
#'
#' @description Calcuate within group conncectivity among disease genes in each networks
#'
#' @param geneset Character string that contains disease-related genes
#' @param disease Character string that states control
#' @param control Character string that states disease
#' @param net.list Output of SortAddLLS, network list to perfrom the analysis
#'
#' @return Dataframe Compactness in each networks
#' @export
#'
#' @examples
#' geneset.conn <- Compactness(net.list = merged.net.list, control = 'hcNETLAS', disease = 'SLE', geneset = sle_genes)
Compactness <- function(net.list = NULL,
                                control = NULL,
                                disease = NULL,
                                geneset = NULL){


  # select celltypes both exist in hcNETLAS and disease
  items <- names(net.list)

  celltypes_ctl <- gsub(paste0(control,"_"),"",items[grep(control,items)])
  celltypes_dis <- gsub(paste0(disease,"_"),"",items[grep(disease,items)])

  celltypes_exist <- intersect(celltypes_dis, celltypes_ctl)


  net.list_perSample <- list()
  for (sample in c(control,disease)){
    net.list_sample <- paste(sample,celltypes_exist,sep="_")
    net.list_perSample[[sample]] <- net.list[net.list_sample]
    names(net.list_perSample[[sample]]) <- gsub(paste0(sample,"_"),"",names(net.list_perSample[[sample]]))
  }

  #these are for scNET
  geneset.within_list <- list()
  for (i in 1:length(net.list_perSample)){
    sample <- names(net.list_perSample)[i]
    geneset.within_list[[sample]] <- unlist(lapply(net.list_perSample[[sample]], function(net) nrow(net[(net[,1] %in% geneset & net[,2] %in% geneset), ])))
  }

  #bind to DF
  df.net <- rbind(geneset.within_list[[1]],geneset.within_list[[2]])
  rownames(df.net) <- names(geneset.within_list)

  return(df.net)
}


