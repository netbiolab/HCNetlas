#' DiffPR.HCNetlas
#'
#' @description Get differnence of normalized centrality values from the output of GetCentrality()
#'
#' @param rank.df.final Output of CombinePercRank
#' @param disease Character string that states control
#' @param control Character string that states disease
#'
#' @return Dataframe of celltypes and their diffPR values of each genes in the CGNs
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' strength.list <- GetCentrality(net.list = merged.net.list)
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
#' DiffPR.HCNetlas(rank.df.final, control = 'HCNetlas', disease = 'SLE')
DiffPR.HCNetlas <- function(rank.df.final = NULL,
                            control = NULL,
                            disease = NULL){

  items <- colnames(rank.df.final)

  spl_items <- str_split(names(rank.df.final),'_',n=2)
  diagnosis_v <- unlist(lapply(spl_items, function(x) x[1]))


  #check each parameter input
  if (!(control %in% unique(diagnosis_v))){
    print(paste('condition', control, 'does not exist in input rank.df.final'))
    stop('this should contain condition represented in rank.df.final')
  }

  if (!(disease %in% unique(diagnosis_v))){
    print(paste('condition', disease, 'does not exist in input rank.df.final'))
    stop('this should contain condition represented in rank.df.final')
  }

  celltypes_ctl <- gsub(paste0(control,"_"),"",items[grep(control,items)])
  celltypes_dis <- gsub(paste0(disease,"_"),"",items[grep(disease,items)])

  celltypes.analyze <- intersect(celltypes_dis, celltypes_ctl)

  rank.list <- list()
  for (celltype in celltypes.analyze){
    celltype_condition_cols <- c(paste(control, celltype, sep = '_'), paste(disease, celltype,sep = '_'), paste(celltype, control, sep = '_'), paste(celltype, disease,sep = '_'))
    colnames.in <- celltype_condition_cols[celltype_condition_cols %in% colnames(rank.df.final)]
    df <- rank.df.final[,colnames.in]
    #its percentile rank...DIS-CTL lets higher value have higher centrality in DIS
    df$diff.rank <- df[,2] - df[,1]

    #make a list for each celltype get genes that have the most rank differential
    rank.list[[celltype]] <- df$diff.rank
    names(rank.list[[celltype]]) <- rownames(df)
    rank.list[[celltype]] <- rank.list[[celltype]][order(abs(rank.list[[celltype]]), decreasing = T)]
  }

  #check that all element of list contain same genes(all genes)
  df.final <- data.frame(matrix(nrow=length(rank.list[[1]]), ncol=2*length(rank.list)))

  for (i in seq(1,length(rank.list))){
    t <- 2*i - 1
    df.final[,t] <- names(rank.list[[i]])
    df.final[,t+1] <- rank.list[[i]]
  }

  column.name <- names(rank.list)
  final.column.names <- vector()
  for (i in seq_along(column.name)){
    celltype <- column.name[i]
    column.name2 <- paste(column.name[[i]], paste0(disease,'-',control),sep = '_')
    index = 2*i -1
    final.column.names[index] <- celltype
    final.column.names[index+1] <- column.name2
  }
  colnames(df.final) <- final.column.names

  return(df.final)

}


#' FindDiffHub.HCNetlas
#'
#' @description Calcuate statistical signficance of diffPR values calculated with DiffPR()
#'
#' @param rank.df.final Output from CombinePercRank
#' @param dis.meta disease metadata data.frame used for this analayis
#' @param celltypes column name where celltypes annotaiton are stored
#' @param disease Character string that states control
#' @param control Character string that states disease
#' @param net.list Output of SortAddLLS, network list to perfrom the analysis
#' @param centrality the centrality method used in the getCentrality. both must match! this is used to find the null distribution
#' @param q.method Method to perfrom multiple testing correction. accepted values are: c("BH","holm","hochberg","bonferroni","BY","fdr","none")
#' @param min.cells minimum cells required to perform diff network analysis default value is 500
#'
#' @return Dataframe Percentile rank of centrality of each networks, their difference, and their statistical significance
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' strength.list <- GetCentrality(net.list = merged.net.list)
#' rank.df.final <- CombinePercRank(perc.rank.list = strength.list)
#' diffPR.df.sig <- FindDiffHub.HCNetlas(rank.df.final, control = 'HCNetlas', disease='SLE', net.list = merged.net.list, centrality = "degree", q.method = "BH", min.cells = 1000, dis.meta = meta_select, celltypes = "celltypes_merged")
FindDiffHub.HCNetlas <- function(rank.df.final = NULL,
                                 disease = NULL,
                                 control = NULL,
                                 net.list = NULL,
                                 centrality = "degree",
                                 q.method = "BH",
                                 min.cells = 1000,
                                 dis.meta = NULL,
                                 celltypes = NULL){


  # select celltypes both exist in HCNetlas and disease
  items <- colnames(rank.df.final)

  spl_items <- str_split(names(rank.df.final),'_',n=2)
  diagnosis_v <- unlist(lapply(spl_items, function(x) x[1]))

  #check each parameter input
  if (!(celltypes %in% names(meta))){
    print(paste('celltypes column', celltypes, 'does not exist in metadata'))
    stop('this column should only have celltype names, not conditions...')
  }

  if (!(control %in% unique(diagnosis_v))){
    print(paste('condition', control, 'does not exist in input rank.df.final'))
    stop('this should contain condition represented in rank.df.final')
  }

  if (!(disease %in% unique(diagnosis_v))){
    print(paste('condition', disease, 'does not exist in input rank.df.final'))
    stop('this should contain condition represented in rank.df.final')
  }

  if (!((type(min.cells) == "double") && min.cells > 0)){
    stop('wrong input for minimum cell counts')
  }

  if (!(q.method %in% c("BH","holm","hochberg","bonferroni","BY","fdr","none"))){
    stop('wrong input for Q-value method')
  }

  celltypes_ctl <- gsub(paste0(control,"_"),"",items[grep(control,items)])
  celltypes_dis <- gsub(paste0(disease,"_"),"",items[grep(disease,items)])

  celltypes_exist <- intersect(celltypes_dis, celltypes_ctl)

  #set 1 core for data.table frank
  data.table::setDTthreads(threads = 1)


  final.df.list <- list()

  #check that each celltype pass minimum threshold (this only checks for disease networks since all HCNetlas CGNs were built with more than 1,000 cell counts)
  celltypes.analyze <- vector()
  for (celltype in celltypes_exist){
    disease.cells <- meta[(meta[,celltypes] == celltype),]

    if (nrow(disease.cells) >= min.cells){
      celltypes.analyze <- c(celltypes.analyze, celltype)
    }
    else{
      print(paste(celltype,':',nrow(disease.cells),'Disease cells. Lower than min.cells threshold, skipping diffHub analysis...'))
    }
  }

  for (celltype in celltypes.analyze){
    #progress bar
    disease.cells.analyze <- meta[(meta[,celltypes] == celltype),]

    print(paste("Finding", celltype, "DiffHubs between HCNetlas CGNs and Disease CGNs with ", nrow(disease.cells.analyze),"cells..."))

    #get diffPR.df of ctrl disase for each celltype
    celltype_condition_cols <- c(paste(control, celltype, sep = '_'), paste(disease, celltype,sep = '_'), paste(celltype, control, sep = '_'), paste(celltype, disease,sep = '_'))
    colnames.in <- celltype_condition_cols[celltype_condition_cols %in% colnames(rank.df.final)]
    df <- rank.df.final[,colnames.in]

    #get disease net # this is not needed in the null distribution generation
    disease.net.name <- colnames.in[grep(disease, colnames.in)]
    disease.net <- igraph::graph_from_data_frame(net.list[[disease.net.name]], directed = F)
    shuffled.weight1 <- sample(E(disease.net)$LLS) #rewire does not suporte weight..so we are going to shuffle both node and edges, while preserving topology

    #get control net
    control.net.name <- colnames.in[!(colnames.in %in% disease.net.name)]
    control.net <- igraph::graph_from_data_frame(net.list[[control.net.name]], directed = F)
    shuffled.weight2 <- sample(E(control.net)$LLS) #rewire does not suporte weight..so we are going to shuffle both node and edges, while preserving topology

    #get df where at least one column has non-zero value. we don't need to evaluate diffPR zero..
    df.f <- df[!(df[,1] ==0 & df[,2] == 0),]

    #its percentile rank...DIS-CTL lets higher value have higher centrality in Disease CGNs compare to HCNetlas CGNs
    df.f$gene <- rownames(df.f)
    df.f$DiffPR <- df.f[,2] - df.f[,1]

    #make this once for all gene comparision, we assume that null diffPR distribution is same for all genes
    #to get significance of our diffPR value, we make n random diffPR and compare
    null.distribution <- vector()
    #make two random network and append diffPR values to null.distribution
    while (length(null.distribution) < 1000000){
      # we make two random network from control..and find null diffPR values

      #random control network 1
      random.control1 <- rewire(control.net, with = each_edge(0.9, loops = F))
      #random.disease <- rewire(control.net, with = keeping_degseq(niter = vcount(control.net * 50)))
      E(random.control1)$LLS <- shuffled.weight2
      random.net.control1 <- igraph::as_data_frame(random.control1)
      random.cent.control1 <- GetCentrality(method = centrality, net.list =list(disease.net.name = random.net.control1))
      random.cent.control1 <- as.data.frame(random.cent.control1[[1]])
      random.cent.control1$gene <- rownames(random.cent.control1)

      #random control network
      random.control2 <- rewire(control.net, with = each_edge(0.9, loops = F))
      #random.control <- rewire(control.net, with = keeping_degseq(niter = vcount(control.net * 50)))
      E(random.control2)$LLS <- shuffled.weight2
      random.net.control2 <- igraph::as_data_frame(random.control2)
      random.cent.control2 <- GetCentrality(method = centrality, net.list =list(control.net.name = random.net.control2))
      random.cent.control2 <- as.data.frame(random.cent.control2[[1]])
      random.cent.control2$gene <- rownames(random.cent.control2)

      #bind two centrality
      diffPR.random <- dplyr::full_join(random.cent.control1, random.cent.control2, by='gene')
      diffPR.random[is.na(diffPR.random)] <- 0 #replace NA with 0
      diffPR.values <- as.numeric(diffPR.random[,1]) - as.numeric(diffPR.random[,3])
      null.distribution <- c(null.distribution, diffPR.values)
    }

    #hist(null.distribution, col='grey')

    #remove RPS RPL MRPS MRPL from union geneset of df.f to save time..these genes will be disregarded
    genes.all <- df.f$gene
    ribo.genes <-  genes.all[grep('^RPS[0-9]*|^RPL[0-9]*', genes.all)]
    mito.genes <- genes.all[grep('^MRPS[0-9]*|^MRPL[0-9]*', genes.all)]
    nduf.genes <- genes.all[grep('^NDUF', genes.all)]
    bad.genes <- c(ribo.genes, mito.genes, nduf.genes)

    df.f.f <- df.f[!(df.f$gene %in% bad.genes),]

    # iterate over union set of genes and get pvalue
    pb = txtProgressBar(min = 0, max = nrow(df.f.f), style = 3)
    pvalue.all <- vector()
    for (g in 1:nrow(df.f.f)){
      Sys.sleep(0.05)
      setTxtProgressBar(pb,g)
      gene <- rownames(df.f)[g]
      diffPR.gene <- df.f.f$DiffPR[g]

      #get pvalue, $frank is 10times faster
      distribution.all <- c(null.distribution, diffPR.gene)
      pvalue <- data.table::frank(-abs(distribution.all), ties.method = "min")[length(distribution.all)] / length(distribution.all) # ties are averaged..highly unlikely here
      pvalue.all[g] <- pvalue
    }
    close(pb)

    #make pvalue fo reach celltype
    df.f.f$pvalue <- pvalue.all
    df.f.f$qvalue <- p.adjust(df.f.f$pvalue, method = q.method, n = nrow(df.f.f))
    df.f.f$celltype <- rep(celltype, nrow(df.f.f))
    colnames(df.f.f) <- c("Control_scHumanNet", "Disease_scHumanNet", "gene", "diffPR", "pvalue", "qvalue", "celltype")

    final.df.list[[celltype]] <- df.f.f
  }

  names(final.df.list) <- NULL
  diffPR.df.result <- do.call("rbind", final.df.list)

  return(diffPR.df.result)
}
