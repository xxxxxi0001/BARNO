#' Use loess select high variance gene
#'
#' @param df Data frame with gene expression
#' @param p_value Threshold for variance gene 
#' @param span Loess
#' 
#' @return All high variance gene
#' 
#' @export
hvg_selection<-function (df,
                         p_value=0.2,
                         span=0.3) {
  mean_var_gene<-data.frame(
    Gene_Name=colnames(df),
    Mean=colMeans(df),
    Variance=apply(df,2,var)
  )
  mean_var_gene$Mean_log10<-log10(mean_var_gene$Mean+1)
  mean_var_gene$Variance_log10<-log10(mean_var_gene$Variance+1)

  loess_model<-loess(Variance_log10~Mean_log10, data=mean_var_gene, span=span)
  
  mean_var_gene$Expected_Variance<-loess_model$fitted

  mean_var_gene$Residual<-mean_var_gene$Variance_log10-mean_var_gene$Expected_Variance

  mean_var_gene$scaled_residual<-scale(mean_var_gene$Residual)
  mean_var_gene$p_value_cal<-pnorm(mean_var_gene$scaled_residual,lower.tail = F)

  hvg_list<-subset(mean_var_gene, p_value_cal < p_value)
  hvg_gene_names<-hvg_list$Gene_Name

  return(hvg_gene_names)
}

#' Use kNN impute missing value
#'
#' @param df Data frame with gene expression
#' @param hvg_list High variance gene calculate before
#' @param k The number of k you choose for knn imputation
#' 
#' @return Imputated expression data frame
#' 
#' @export
pca_knn_imputation<-function (df,
                              hvg_list,
                              k=5) {
  df_hvg<-df[,hvg_list]
  pca_results<-prcomp(df_hvg)
  pca_df<-as.data.frame(pca_results$x[,1:min(30,ncol(pca_results$x))])
  message("Before Imputation:", sum(df_hvg == 0), " zeros.")
  df_distance<-as.matrix(dist(pca_df, method="euclidean"))

  k_num<-k

  for (i in 1:nrow(df_hvg)) {

    zeros<-which(df_hvg[i,]==0)
    if (length(zeros) !=0) {

      cell_name<-rownames(df_hvg)[i]

      cell_distance<-df_distance[cell_name,]

      nearest_k<-sort(cell_distance)[2:(k+1)]

      neighbor_name<-names(nearest_k)

      neighbor<-df_hvg[neighbor_name, zeros, drop = FALSE]
      df_hvg[i,zeros]<-colMeans(neighbor)
    }
  }
  
  message("After Imputation: ", sum(df_hvg == 0), " zeros.")
  
  return(df_hvg)
}

#' Use Silhouette find best k for k mean clustering
#'
#' @param df_umap UMAP data frame
#' @param test_k Number of k you wanna test
#' @param nstart K-mean clustering parameter
#' 
#' @return Best k
#' 
#' @export
find_best_k_kmean<-function(df_umap, 
                            test_k=10,
                            nstart=10) {
  best_k <- 2
  max_avg_sj<--Inf
  df_distance<-as.matrix(dist(df_umap, method="euclidean"))
  for (k in 2:test_k) {
    kmean_result<-kmeans(df_umap, centers = k, nstart = nstart)
    clusters<-kmean_result$cluster
    sj_all<-numeric(nrow(df_umap))
    # for each cell
    for (i in 1:nrow(df_umap)) {
      my_cluster<-clusters[i]
      cluster_all<-which(clusters == my_cluster)
      cluster_all<-cluster_all[cluster_all!=i]
      if (length(cluster_all)>0) {
        aj<-mean(df_distance[i,cluster_all])
      }
      else{
        aj<-0
      }
      bj<-Inf
      for (j in 1:k) {
        if (j != my_cluster) {
          other_cluster<-which(clusters==j)
          if (length(other_cluster)>0){
            mean_distance_bj<-mean(df_distance[i,other_cluster])
            if(mean_distance_bj<bj){
              bj<-mean_distance_bj
            }
          }
        }
      }
      sj_all[i]<-(bj-aj)/max(aj,bj)
    }
    avg_sj<-mean(sj_all)
    if (avg_sj>max_avg_sj){
      max_avg_sj<-avg_sj
      best_k<-k
    }
  }
  message("Through Silhouette, the best k is ",best_k," and the Silhouette is ",round(max_avg_sj,2))
  return(best_k)
}

#' Calculate pseudotime for cell progression
#'
#' @param RDS_file RDS file with tpm 
#' @param nstart K-mean clustering parameter
#' 
#' @return Pseudotime
#' 
#' @export
load_pseudotime<-function(RDS_file,
                          nstart=10){
  
  # load necessary packages for pseudotime
  if (!requireNamespace("BiocManager")) {
    stop("Need package BiocManager")
  }
  if (!requireNamespace("slingshot")) {
    stop("Need package slingshot")
  }
  if (!requireNamespace("RcisTarget")) {
    stop("Need package RcisTarget")
  }
  if (!requireNamespace("umap")) {
    stop("Loading package umap")
  }

  data("motifAnnotations_hgnc", package="RcisTarget")
  allTFs <- unique(motifAnnotations$TF)
  df<-readRDS(RDS_file)
  df_tpm<-df$tpm
  if (nrow(df_tpm) != length(df$cells)) {
    df_tpm <- t(df_tpm) 
  }
  hvg_gene<-hvg_selection(df_tpm)

  genes_to_run <- intersect(unique(c(hvg_gene, allTFs)), colnames(df_tpm))
  message("High Variance Gene Number (with TF): ", length(genes_to_run))
  df_imputated<-pca_knn_imputation(df_tpm,genes_to_run)
  df_imputated_scaled <- scale(df_imputated)
  suppressMessages({umap_result<-umap(df_imputated_scaled)})
  umap_matrix <- umap_result$layout
  best_k<-find_best_k_kmean(umap_matrix)
  kmean_result<-kmeans(umap_matrix, centers = best_k, nstart = nstart)
  df_run<-data.frame(
    UMAP1=umap_matrix[,1],
    UMAP2=umap_matrix[,2],
    clusters = as.character(kmean_result$cluster),
    treated=df$treated
  )
  
  rownames(df_run) <- rownames(df_imputated)
  slingshot_result<-slingshot(data = umap_matrix, clusterLabels = df_run$clusters, start.clus = '1')
  pt_matrix <- slingshot::slingPseudotime(slingshot_result)
  pt_vector <- rowMeans(pt_matrix, na.rm = TRUE)
  df_run$pseudotime<-pt_vector
  df_run <- df_run[!is.na(df_run$pseudotime), ]
  
  message(length(df$cells)-nrow(df_run), " cells are droped in pseudotime calculation.")
  
  return(list(
    meta_data=df_run,
    imputed_matrix=df_imputated,
    slingshot_result=slingshot_result
  ))
}

#' Calculate state score as pseudo-pseudotime for cell progression
#'
#' @param RDS_file RDS file with tpm 
#' @param gene_list A list of gene that represent your selection of cell type
#' 
#' @return Pseudo-pseudotime
#' 
#' @export
load_state_score<-function(RDS_file,
                           gene_list){
  if (!requireNamespace("BiocManager")) {
    stop("Need package BiocManager")
  }
  if (!requireNamespace("Seurat")) {
    stop("Need package Seurat")
  }
  if (!requireNamespace("RcisTarget")) {
    stop("Need package RcisTarget")
  }
  data("motifAnnotations_hgnc", package="RcisTarget")
  allTFs <- unique(motifAnnotations$TF)
  df<-readRDS(RDS_file)
  df_tpm<-df$tpm

  if (nrow(df_tpm) != length(df$cells)) {
    df_tpm <- t(df_tpm) 
  }

  seurat_object<-CreateSeuratObject(counts=t(df_tpm))
  seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
  target_features <- list(toupper(gene_list))
  seurat_object<-AddModuleScore(seurat_object, features = target_features, name = "Target_Score")
  scores<-seurat_object$Target_Score1
  alt_pseudotime<-(max(scores)-scores)/(max(scores)-min(scores))
  hvg_gene<-hvg_selection(df_tpm)
  genes_to_run <- intersect(unique(c(hvg_gene, allTFs)), colnames(df_tpm))
  message("High Variance Gene Number (with TF): ", length(genes_to_run))
  df_imputated<-pca_knn_imputation(df_tpm,genes_to_run)
  
  df_run<-data.frame(
    treated=df$treated,
    pseudotime=alt_pseudotime,
    row.names=rownames(df_tpm)
  )
  
  return(list(
    meta_data=df_run,
    imputed_matrix=df_imputated,
    raw_tpm = df$tpm
  ))
}

#' Find peak of loess
#'
#' @param loess_fit The loess construct before
#' @param x Gene expression
#' @param y Seurat Score
#' 
#' @return Peak of loess
#' 
#' @export
find_peak<-function(loess_fit,
                    x,
                    y){
  x_pred<-seq(min(x),max(x),length.out=length(x))
  y_pred<-predict(loess_fit,newdata=data.frame(x=x_pred))
  peak_index<-which.max(y_pred)
  T_peak<-x_pred[peak_index]
  return(T_peak)
}

#' Caculate weight based on pseudotime
#'
#' @param meta_data Data frame that has all calculated value we calculate before
#' @param imputed_matrix Imputed gene matrix
#' @return Peak information & correlation informaiton
#' @export
dynamic_weight_initialization<-function(meta_data,
                                        imputed_matrix){
  
  genes <- colnames(imputed_matrix)
  results<-data.frame(
    gene=genes,
    T_peak_treated=NA,
    T_peak_untreated=NA,
    R_treated=NA,
    R_untreated=NA
  )

  groups_in_data<-unique(meta_data$treated)
  for (i in 1:ncol(imputed_matrix)){
    gene_expression<-imputed_matrix[,i]
    for (group in groups_in_data){
      cells<-rownames(meta_data)[which(meta_data$treated == group & !is.na(meta_data$pseudotime))]
      if(length(cells)<5) {next}
      x<-as.numeric(meta_data[cells,"pseudotime"])
      y<-as.numeric(gene_expression[cells])
      fit<-loess(y~x,span=0.75,control = loess.control(surface="direct"))
      T_peak<-find_peak(fit,x,y)
      correlation<-cor(x,y,method="spearman")
      if (group == FALSE) {
        column_peak <- "T_peak_untreated"
        column_correlation <- "R_untreated"
      } else if (group == TRUE) {
        column_peak <- "T_peak_treated"
        column_correlation <- "R_treated"
      }
      results[i, column_peak] <- T_peak
      results[i, column_correlation] <- correlation
    }
  }
  return(results)
}

#' Caculate penalty score based on expression's dispersion
#'
#' @param tf_expression Data frame that has all calculated value we calculate before
#' @param pca_coordinates Imputed gene matrix
#' @param global_disperse Global dispersion
#' @param quantile_por Threshold set for highly expression cluster group
#' @param threshold Threshold set for portion
#' @param alpha Penalty degree
#' @param amplify Amplification degree of portion
#' 
#' @return Penalty score of batch effect
#' 
#' @export
penalty_score_batch_effect<-function(tf_expression,
                                     pca_coordinates,
                                     global_disperse,
                                     quantile_por=0.9,
                                     threshold=0.5,
                                     alpha=2,
                                     amplify=3){
  true_cell_index<-which(tf_expression>quantile(tf_expression,quantile_por))
  if (length(true_cell_index)<3) {
    return(0)
  }
  cell_coordinate<-pca_coordinates[true_cell_index, ,drop=FALSE]
  cell_center<-colMeans(cell_coordinate)
  cell_disperse<-mean(sqrt(rowSums(sweep(cell_coordinate,2,cell_center,"-")^2)))
  portion<-cell_disperse/global_disperse
  portion<-portion^amplify
  if(portion>threshold){
    return(1)
  }else{
    return((portion/threshold)^alpha)
  }
}

#' Caculate penalty score for each tf
#'
#' @param genie_matrix The GENIE3 weight matrix
#' @param df_tpm Expression matrix (original)
#' @param quantile_por Threshold set for highly expression cluster group
#' @param threshold Threshold set for portion
#' @param alpha Penalty degree
#' @param amplify Amplification degree of portion
#' 
#' @return Penalty score of batch effect
#' 
#' @export
generate_penalty_score<-function(genie_matrix,
                                 df_tpm,
                                 quantile_por=0.9,
                                 threshold=0.5,
                                 alpha=2,
                                 amplify=3) {
  tfs<-rownames(genie_matrix)
  df_tpm<-log2(df_tpm+1)
  gene_variance<-apply(df_tpm,1,var)
  top_2000<-names(sort(gene_variance,decreasing=TRUE))[1:2000]
  top_2000_df<-df_tpm[top_2000,]
  pca_result<-prcomp(t(top_2000_df),center=TRUE, scale.=TRUE)
  pca_coordinates<-pca_result$x[,1:5]
  global_center<-colMeans(pca_coordinates)
  global_disperse<-mean(sqrt(rowSums(sweep(pca_coordinates,2,global_center,"-")^2)))
  penalty_score<-numeric(length(tfs))
  names(penalty_score)<-tfs
  for (tf in tfs) {
    tf_expression<-df_tpm[tf,]
    penalty_score[tf]<-penalty_score_batch_effect(tf_expression,pca_coordinates,global_disperse,quantile_por,threshold,alpha,amplify)
  }
  penalty_df<-data.frame(
    TF=names(penalty_score),
    panelty_score=penalty_score
  )
  return(penalty_df)
}

#' Caculate weight by calling previous function one by one
#'
#' @param RDS_file RDS file with treated T/F, tpm
#' @param gene_list List of gene selected for scoring
#' @param genie_file If genie calculated before, import file
#' @param batch_effect Whether select slingshot or seurat
#' @param alpha Penalty degree
#' @param nTree Number of tree for genie3
#' @param nCores Number of core for genie3
#' 
#' @return Weight matrix and other additional info might be useful
#' 
#' @export
TR_weight<-function(RDS_file,
                    gene_list,
                    genie_file=NULL,
                    batch_effect=TRUE,
                    alpha=0.5,
                    nTree=500,
                    nCores=1){
  
  message("[1/4] TR_Weight: Loading Pseudotime...")
  
  if (batch_effect) {
    result <- load_state_score(RDS_file,gene_list)
  }else {result<-load_pseudotime(RDS_file)}
  
  meta_data<-result$meta_data
  imputed_matrix<-result$imputed_matrix
  
  message("[2/4] TR_Weight: Loading GENIE3 weight matrix...")
  if (is.null(genie_file)) {
    if (!requireNamespace("BiocManager")) {
      stop("Need Package BiocManager")
    }
    if (!requireNamespace("GENIE3")) {
      stop("Need Package GENIE3")
    }
    message("Preparing run GENIE3, this may takes a while")
    expression_matrix<-as.matrix(t(imputed_matrix))
    data("motifAnnotations_hgnc", package="RcisTarget")
    allTFs <- unique(motifAnnotations$TF)
    inputTFs<-intersect(allTFs, rownames(expression_matrix))
    
    genie_matrix<-GENIE3(
      exprMatrix = expression_matrix,
      regulators=inputTFs,
      nTrees = nTree,
      nCores = nCores,
      verbose=TRUE)
    saveRDS(genie_matrix,file="genie_matrix.rds")
  }
  else {genie_matrix<-readRDS(genie_file)}
  
  message("[3/4] TR_Weight: Initializing dynamix weight matrix...")

  dynamic_matrix<-dynamic_weight_initialization(meta_data,imputed_matrix)
  
  message("[4/4] TR_Weight: Constructing Dynamic Weight Matrix...")

  tf_names <- rownames(genie_matrix)
  target_names <- colnames(genie_matrix)

  tf_peak_Treated<-dynamic_matrix$T_peak_treated[match(tf_names, dynamic_matrix$gene)]
  target_peak_Treated<-dynamic_matrix$T_peak_treated[match(target_names, dynamic_matrix$gene)]
  target_correlation_Treated<-dynamic_matrix$R_treated[match(target_names, dynamic_matrix$gene)]

  tf_peak_Untreated<-dynamic_matrix$T_peak_untreated[match(tf_names, dynamic_matrix$gene)]
  target_peak_Untreated<-dynamic_matrix$T_peak_untreated[match(target_names, dynamic_matrix$gene)]
  target_correlation_Untreated<-dynamic_matrix$R_untreated[match(target_names, dynamic_matrix$gene)]

  delta_T_treated<-outer(tf_peak_Treated, target_peak_Treated, function(tf, target) target - tf)
  delta_T_untreated<-outer(tf_peak_Untreated, target_peak_Untreated, function(tf, target) target - tf)

  a<-alpha
  penalty_Treated<-ifelse(delta_T_treated >= 0, 1, exp(a * delta_T_treated))
  penalty_Untreated<-ifelse(delta_T_untreated >= 0, 1, exp(a * delta_T_untreated))
  penalty_Treated[is.na(penalty_Treated)]<-0
  penalty_Untreated[is.na(penalty_Untreated)]<-0

  cor_Treated<-matrix(rep(ifelse(is.na(target_correlation_Treated), 0, target_correlation_Treated), each = length(tf_names)), nrow = length(tf_names))
  cor_Untreated<-matrix(rep(ifelse(is.na(target_correlation_Untreated), 0, target_correlation_Untreated), each = length(tf_names)), nrow = length(tf_names))

  weight_treated<-genie_matrix * penalty_Treated * cor_Treated
  weight_untreated<-genie_matrix * penalty_Untreated * cor_Untreated
  weight_treated<-scale(weight_treated)
  weight_untreated<-scale(weight_untreated)
  TR_weight_matrix<-weight_treated-weight_untreated

  if(batch_effect){
    message("Punishing TF with Batch Effect...")
    df_tpm<-result$raw_tpm
    penalty_df<-generate_penalty_score(genie_matrix,df_tpm)
    penalty_aligned<-penalty_df$panelty_score[match(rownames(TR_weight_matrix), penalty_df$TF)]
    TR_weight_matrix<-sweep(TR_weight_matrix,1,penalty_aligned, "*")
  }
  distance_matrix<-1/(abs(TR_weight_matrix)+1e-6)
  
  return(list(
    TR_weight_matrix=TR_weight_matrix,
    distance_matrix=distance_matrix,
    meta_data=meta_data,
    imputed_matrix=imputed_matrix,
    penalty_df=penalty_df
  )
  )
}

#' Get top expressed TF based on previously calculated weight
#'
#' @param weight_matrix Weight matrix you calculated before
#' 
#' @return Top expressed TF
#' 
#' @export
top_TRweight_genes<-function(weight_matrix){
  tr_matrix <- weight_matrix$TR_weight_matrix
  tr_df <- as.data.frame(as.table(tr_matrix))
  colnames(tr_df) <- c("TF", "Target", "TR_Weight")
  top_pos<-tr_df[order(tr_df$TR_Weight, decreasing=TRUE), ]
  top_neg<-tr_df[order(tr_df$TR_Weight, decreasing=FALSE), ]
  pos<-tr_df[tr_df$TR_Weight>0,]
  tf_pos_summary<-aggregate(TR_Weight ~ TF, data = pos, sum)
  tf_pos_summary<-tf_pos_summary[order(tf_pos_summary$TR_Weight, decreasing = TRUE), ]
  neg<-tr_df[tr_df$TR_Weight<0,]
  tf_neg_summary<-aggregate(TR_Weight ~ TF, data = neg, sum)
  tf_neg_summary<-tf_neg_summary[order(tf_neg_summary$TR_Weight, decreasing = FALSE), ]
  
  return(list(
    top_positive_genes=top_pos,
    top_negative_genes=top_neg,
    top_positive_summary=tf_pos_summary,
    top_negative_summary=tf_neg_summary
  ))
}