#' Select co-factor for selected tf
#'
#' @param weight_matrix Weight Matrix's name
#' @param TF Selected TF for km construction
#' @param top_percentage Treshold for co-factor
#' 
#' @return co-factor
#' 
#' @export
select_co_factor<-function(weight_matrix,
                           TF,
                           top_percentage=0.9){
  expression<-weight_matrix$imputed_matrix
  meta_data<-weight_matrix$meta_data
  treated_index<-which(meta_data$treated==TRUE)
  top_expression<-quantile(expression[treated_index,TF],top_percentage)
  high_expression_cells<-rownames(expression)[treated_index][expression[treated_index,TF]>top_expression]
  low_expression_cells<-setdiff(rownames(expression)[treated_index], high_expression_cells)

  mean_high<-colMeans(expression[high_expression_cells, ])
  mean_low<-colMeans(expression[low_expression_cells, ])

  df_fc<-data.frame(
    Gene=colnames(expression),
    Log2FC=log2(mean_high+1)-log2(mean_low+1)
  )

  df_fc<-df_fc[order(df_fc$Log2FC, decreasing=TRUE),]

  return(df_fc)
}

#' Construct km-curve
#'
#' @param signature_genes Genes selected for km construction
#' @param patient_RDS RDS with survival & patient info
#' @param TF Selected TF for km construction
#' @param cut Default median, otherwise optimal
#' 
#' @return km-curve
#' 
#' @export
km_curve<-function(signature_genes,
                   patient_RDS,
                   TF,
                   cut="median") {
  if (!requireNamespace("survival")) {
    stop("Need package survival")
  }
  if (!requireNamespace("survminer")) {
    stop("Need package survminer")
  }

  real<-readRDS(patient_RDS)
  real_expression<-real$tpm[intersect(signature_genes, rownames(real$tpm)), ]
  normalized_expression<-t(scale(t(log2(real_expression + 1))))

  survival_df<-data.frame(
    time=as.numeric(real$survival[, 1]),
    status=as.numeric(real$survival[, 2]),
    score=colMeans(normalized_expression)
  )
  if (cut!="median"){
    res.cut<-surv_cutpoint(survival_df, time = "time", event = "status", variables = "score")
    survival_df <- surv_categorize(res.cut)
    km_fit <- survfit(Surv(time, status) ~ score, data = survival_df)
    
  } else {
    median_score<-median(survival_df$score, na.rm=TRUE)
    survival_df$score_group<-ifelse(survival_df$score >= median_score, "High", "Low")
    survival_df$score_group<-factor(survival_df$score_group, levels = c("High", "Low"))

    km_fit<-survfit(Surv(time, status) ~ score_group, data = survival_df)
  }
  
  
  km_curve_visualization(km_fit,survival_df,TF)
  
  return(km_fit)
}