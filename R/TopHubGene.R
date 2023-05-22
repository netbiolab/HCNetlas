#' FindNeighbors
#'
#' @description Find direct neibhbors for each CGNs in corresponding top hub genes
#'
#' @param top.df Output of TopHub, a dataframe with top hubness genes
#' @param net.list List of networks. Should contain all networks from top.df
#'
#' @return List Characters of direct neighboring nodes and top hub genes in each network
#' @export
#'
#' @examples
#' direct.neighbors <- FindNeighbors(top.df=sel_top.df, net.list =merged.net.list)
FindNeighbors <- function(top.df = NULL,
                          net.list = NULL){

  if (!(all(colnames(top.df) %in% names(net.list)))){
    cat("Not all columns in top.df are represented in network list")
  } else {
    direct.neighbors <- list()

    for (i in 1:ncol(top.df)){
      cgn <- colnames(top.df)[i]
      top.nodes <- top.df[,cgn]

      network <- net.list[[cgn]]
      network_connected_to_topnodes <- network[(network$'1' %in% top.nodes | network$'2' %in% top.nodes),]
      direct.neighbors[[cgn]] <- union(network_connected_to_topnodes$'1', network_connected_to_topnodes$'2')
    }
    return (direct.neighbors)
  }
}




#' GSAplot
#'
#' @description Plot Gene set enrichment plot using EnrichR
#'
#' @param genes Vector a gene set of interest
#' @param database Character indicating which pathway to use. Should be db contained in EnrichR
#' @param title Title label of GSEA plot
#' @param top.term Number of terms that will be represented in a plot
#'
#' @return Plot with top n pathways enriched in a given gene set.
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' go <- GSAplot(top.genes,'GO_Biological_Process_2021',"SLE_Bcell_gene GOBP",10)
GSAplot <- function(genes, database, title, top.term){
  dbs <- listEnrichrDbs()
  dbs <- database
  enrichr <- enrichr(genes, dbs)
  data.gsa <- enrichr[[dbs]]

  aaa <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',1))
  bbb <- as.numeric(sapply(strsplit(data.gsa$Overlap, '/'),'[',2))

  #add column
  data.gsa$overlap_num <- aaa / bbb
  data.gsa$log10_qvalue <- -log10(data.gsa$Adjusted.P.value)

  data.gsa.f <- data.gsa[order(data.gsa$log10_qvalue, decreasing = T),]

  p <- ggplot(data.gsa.f[1:top.term,], aes(x = reorder(Term,log10_qvalue),log10_qvalue,  y = log10_qvalue,
                                           fill = overlap_num)) +
    geom_bar(stat = 'identity', width = 0.9, position = position_dodge(width = 0.1)) +
    geom_hline(yintercept = -log10(0.05), colour = 'red', linetype = 'longdash') +
    coord_flip() +
    scale_fill_gradient(name = 'overlap percentage',low = '#a1c4fd', high = '#ff9a9e') +
    theme_classic() +
    labs(
      title = title,
      y = '-log10(qvalue)',
      x = database
    )

  return(p)
}


