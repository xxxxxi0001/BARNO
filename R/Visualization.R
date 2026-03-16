#' Visualize a directed regulatory cascade
#'
#' @param directed_mst A directed igraph object representing the regulatory cascade.
#' @param TF Parent transcription factor used as the root node.
#' @param pdf_width As parameter name
#' @param pdf_height As parameter name
#' @param png_width As parameter name
#' @param png_height As parameter name
#' @param png_resolution As parameter name
#' @param line_transparency As parameter name
#' @param line_thickness As parameter name
#' @param parent_node_size As parameter name
#' @param intermediate_node_size As parameter name
#' @param lead_node_size As parameter name
#' @param text_size As parameter name
#' @param text_type As parameter name
#' @param text_centrality As parameter name
#' @param text_diagnol As parameter name
#' @param node_color As parameter name
#' @param node_frame_color As parameter name
#' @param arrow_size As parameter name
#' @param aspect_ratio As parameter name
#' 
#' @return GRN graph
#' @export
directed_mst_visualization<-function(directed_mst,
                                     TF,
                                     pdf_width=8,
                                     pdf_height=7,
                                     png_width=2400,
                                     png_height=2100,
                                     png_resolution=300,
                                     line_transparency=1,
                                     line_thickness=1,
                                     parent_node_size=18,
                                     intermediate_node_size=10,
                                     lead_node_size=5,
                                     text_size=0.65,
                                     text_type=2,
                                     text_centrality=-0.5,
                                     text_diagnol=-pi/2,
                                     node_color="white",
                                     node_frame_color="grey",
                                     arrow_size=0.6,
                                     aspect_ratio=0.6){
  
  all_nodes<-igraph::V(directed_mst)$name
  out_degrees<-igraph::degree(directed_mst,mode="out")
  intermediate_tfs<-all_nodes[out_degrees > 0 & all_nodes != TF]

  igraph::E(directed_mst)$color<-ifelse(igraph::E(directed_mst)$TR_Weight > 0, 
                                adjustcolor("#CC3311", alpha.f=line_transparency), 
                                adjustcolor("#0077BB", alpha.f=line_transparency))
  igraph::E(directed_mst)$width <- (abs(igraph::E(directed_mst)$TR_Weight) / 5) * line_thickness

  igraph::V(directed_mst)$level <- ifelse(igraph::V(directed_mst)$name == TF, "Root",
                                  ifelse(igraph::V(directed_mst)$name %in% intermediate_tfs, "Middle", "Leaf"))

  igraph::V(directed_mst)$size <- ifelse(igraph::V(directed_mst)$level == "Root", parent_node_size, 
                                 ifelse(igraph::V(directed_mst)$level == "Middle", intermediate_node_size, lead_node_size))
  options<-c("PDF","PNG","No Export")
  choice<-menu(options,title="How do you want to save result?")
  if (choice==1) {
    pdf(paste0(TF,"_MST_Cascade.pdf"), width=pdf_width, height=pdf_height)
    plot(directed_mst,
         layout = igraph::layout_with_sugiyama(directed_mst)$layout,
         vertex.label = igraph::V(directed_mst)$name,
         vertex.label.cex = text_size,
         vertex.label.font=text_type,
         vertex.label.dist = text_centrality,
         vertex.label.degree = text_diagnol,
         vertex.color = node_color,
         vertex.frame.color = node_frame_color,
         edge.color = igraph::E(directed_mst)$color,
         edge.width = igraph::E(directed_mst)$width,
         edge.arrow.size = arrow_size,
         asp = aspect_ratio,
         main = paste0(TF, " Multi-Layer Regulatory Cascade"),
         sub = "Hierarchy: Parent Transcript Factor -> Intermediate Transcript Factor -> Target Gene")
    dev.off()
    message("Exported")
  }
  else if(choice==2){
    png(paste0(TF,"_MST_Cascade.png"), width=png_width, height=png_height, res=png_resolution)
    plot(directed_mst,
         layout = igraph::layout_with_sugiyama(directed_mst)$layout,
         vertex.label = igraph::V(directed_mst)$name,
         vertex.label.cex = text_size,
         vertex.label.font=text_type,
         vertex.label.dist = text_centrality,
         vertex.label.degree = text_diagnol,
         vertex.color = node_color,
         vertex.frame.color = node_frame_color,
         edge.color = igraph::E(directed_mst)$color,
         edge.width = igraph::E(directed_mst)$width,
         edge.arrow.size = arrow_size,
         asp = aspect_ratio,
         main = paste0(TF, " Multi-Layer Regulatory Cascade"),
         sub = "Hierarchy: Parent Transcript Factor -> Intermediate Transcript Factor -> Target Gene")
    dev.off()
    message("Exported")
  }else{
    plot(directed_mst,
         layout = igraph::layout_with_sugiyama(directed_mst)$layout,
         vertex.label = igraph::V(directed_mst)$name,
         vertex.label.cex = text_size,
         vertex.label.font=text_type,
         vertex.label.dist = text_centrality,
         vertex.label.degree = text_diagnol,
         vertex.color = node_color,
         vertex.frame.color = node_frame_color,
         edge.color = igraph::E(directed_mst)$color,
         edge.width = igraph::E(directed_mst)$width,
         edge.arrow.size = arrow_size,
         asp = aspect_ratio,
         main = paste0(TF, " Multi-Layer Regulatory Cascade"),
         sub = "Hierarchy: Parent Transcript Factor -> Intermediate Transcript Factor -> Target Gene")
  }
}

#' Construct km-curve
#'
#' @param km_fit The km curve constructed earlier
#' @param df Survival df
#' @param Title Usually TF
#' @param pdf_width As parameter name
#' @param pdf_height As parameter name
#' @param png_width As parameter name
#' @param png_height As parameter name
#' @param png_resolution As parameter name
#' 
#' @return visualized km-curve
#' 
#' @importFrom survminer ggsurvplot
#' @export
km_curve_visualization<-function(km_fit,
                                 df,
                                 Title,
                                 pdf_width=10,
                                 pdf_height=6.5,
                                 png_width=2400,
                                 png_height=1600,
                                 png_resolution=300){
  plot<-survminer::ggsurvplot(km_fit,
                   data = df,
                   pval = TRUE, 
                   risk.table = TRUE,
                   conf.int = TRUE,
                   palette = c("#CC3311", "#0077BB"),
                   title = paste0(Title),
                   legend.labs = c("High (Identity Switch)", "Low (Normal)"),
                   xlab = "Time (Years)",
                   risk.table.title="Patients at risk")
  options<-c("PDF","PNG","No Export")
  choice<-menu(options,title="How do you want to save result?")
  if (choice==1) {
    pdf(paste0(Title,"_Survival_Analysis.pdf"), width=pdf_width, height=pdf_height)
    print(plot)
    dev.off()
    message("PDF Exported")
  }
  else if(choice==2){
    png(paste0(Title,"_Survival_Analysis.png"), width=png_width, height=png_height, res=png_resolution)
    print(plot)
    dev.off()
    message("PNG Exported")
  }
  print(plot)
}

#' Construct PCA plot
#'
#' @param target_gene Gene name you want to see distribution
#' @param RDS_file RDS file with tpm
#' @param PCA Number of PCA, 10 is sufficient
#' @param cut Quantile you set
#' @param background_dot_transparency As parameter name
#' @param top10_dot_transparency As parameter name
#' @param background_dot_size As parameter name
#' @param top10_dot_size As parameter name
#' @param eclipse_coverage As parameter name
#' @param eclipse_line_width As parameter name
#' @param axis_title_size As parameter name
#' @param plot_title_size As parameter name
#' @param pdf_width As parameter name
#' @param pdf_height As parameter name
#' @param png_width As parameter name
#' @param png_height As parameter name
#' @param png_resolution As parameter name
#' 
#' @return visualized pca plot
#' 
#' @export
pca_plot_visualization<-function(target_gene,
                                 RDS_file,
                                 PCA=10,
                                 cut=0.9,
                                 background_dot_transparency=0.5,
                                 top10_dot_transparency=0.8,
                                 background_dot_size=1,
                                 top10_dot_size=1.5,
                                 eclipse_coverage=0.95,
                                 eclipse_line_width=1.2,
                                 axis_title_size=12,
                                 plot_title_size=14,
                                 pdf_width=10,
                                 pdf_height=6.5,
                                 png_width=2400,
                                 png_height=1600,
                                 png_resolution=300){
  
  df<-readRDS(RDS_file)
  df_tpm<-df$tpm
  df_log<-t(log2(df_tpm + 1))
  pca_result<-stats::prcomp(df_log, center=TRUE, scale.=TRUE, rank.=PCA)

  pca_df <- data.frame(
    Cell_Name = rownames(df_log),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Expression = as.numeric(df$tpm[target_gene, ]))

  threshold<-stats::quantile(pca_df$Expression, cut)

  pca_df$Group<-ifelse(pca_df$Expression>threshold, "Top 10% Expressed Cells", "Background")

  pca_df$Group <- factor(pca_df$Group, levels=c("Background", "Top 10% Expressed Cells"))
  pca_df<-pca_df[order(pca_df$Group), ]
  
  plot<-ggplot2::ggplot(pca_df, ggplot2::aes(x=PC1, y=PC2, color=Group)) +
    ggplot2::geom_point(ggplot2::aes(alpha=Group, size=Group)) +

    ggplot2::scale_color_manual(values=c("Background"="grey", "Top 10% Expressed Cells" = "#CC3311")) +

    ggplot2::scale_alpha_manual(values=c("Background"=background_dot_transparency, "Top 10% Expressed Cells"=top10_dot_transparency)) +

    ggplot2::scale_size_manual(values=c("Background"=background_dot_size, "Top 10% Expressed Cells"=top10_dot_size)) +
    ggplot2::stat_ellipse(data=subset(pca_df, Group=="Top 10% Expressed Cells"),
                 type="norm", level=eclipse_coverage, linewidth=eclipse_line_width, color="black", linetype="dashed") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face="bold", size=axis_title_size),
      plot.title = ggplot2::element_text(face="bold", size=plot_title_size)
    ) +
    ggplot2::labs(
      title = paste("PCA Spatial of", target_gene, "in 2 Dimension"),
      subtitle = "(We used 5 dimension in real cacluation, this is for visualziation)",
      x = "PC 1", 
      y = "PC 2"
    )

  options<-c("PDF","PNG","No Export")
  choice<-menu(options,title="How do you want to save result?")
  if (choice==1) {
    pdf(paste0(target_gene,"_PC_plot.pdf"), width=pdf_width, height=pdf_height)
    print(plot)
    dev.off()
    message("PDF Exported")
  }
  else if(choice==2){
    png(paste0(target_gene,"_PC_plot.png"), width=png_width, height=png_height, res=png_resolution)
    print(plot)
    dev.off()
    message("PNG Exported")
  } else{
    print(plot)
  }
}

#' Visualize gene expression across patient samples using violin plots
#'
#' @param target_genes Character vector of gene names to visualize.
#' @param RDS_file Path to an RDS file containing TPM expression data and sample information.
#' @param Title Title of the plot.
#' @param Subtitle Subtitle of the plot.
#' @param output_name File name prefix used for exporting the figure.
#' @param pdf_width Width of the exported PDF.
#' @param pdf_height Height of the exported PDF.
#' @param png_width Width of the exported PNG.
#' @param png_height Height of the exported PNG.
#' @param png_resolution Resolution of the exported PNG.
#'
#' @return A violin plot showing gene expression across patient samples.
#' @export
violin_plot_visualization <- function(target_genes,
                                      RDS_file,
                                      Title,
                                      Subtitle,
                                      output_name,
                                      pdf_width = 8,
                                      pdf_height = 5,
                                      png_width = 2400,
                                      png_height = 1500,
                                      png_resolution = 300) {
  
  df <- readRDS(RDS_file)
  df_tpm <- df$tpm
  
  valid_genes <- intersect(target_genes, rownames(df_tpm))
  
  if (length(valid_genes) == 0) {
    stop("None of the target_genes were found in the expression matrix.")
  }
  
  df_expression <- as.data.frame(log2(t(df_tpm[valid_genes, , drop = FALSE] + 1)))
  df_expression$Cell_Name <- rownames(df_expression)
  
  df_plot <- data.frame(
    Cell_Name = colnames(df_tpm),
    Sample = df$samples
  )
  
  df_plot <- merge(df_plot[, c("Cell_Name", "Sample")], df_expression, by = "Cell_Name")
  
  long_plot_data <- reshape2::melt(
    df_plot,
    id.vars = c("Cell_Name", "Sample"),
    measure.vars = valid_genes,
    variable.name = "Gene",
    value.name = "Expression"
  )
  
  violin_plot <- ggplot2::ggplot(
    long_plot_data,
    ggplot2::aes(x = Sample, y = Expression, fill = Sample)
  ) +
    ggplot2::geom_violin(
      scale = "width",
      trim = FALSE,
      alpha = 0.8,
      color = "black",
      linewidth = 0.3
    ) +
    ggplot2::geom_boxplot(
      width = 0.1,
      fill = "white",
      color = "black",
      outlier.shape = NA,
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~ Gene, ncol = 1, scales = "free_y") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 10, color = "black"),
      axis.text.y = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 13, face = "bold", color = "black"),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(face = "italic"),
      strip.text = ggplot2::element_text(face = "bold.italic", size = 13),
      strip.background = ggplot2::element_rect(fill = "grey95", color = "black", linewidth = 1),
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linetype = "dashed")
    ) +
    ggplot2::labs(
      title = Title,
      subtitle = Subtitle,
      x = "Patient Sample ID",
      y = "Expression Level log2(TPM + 1)"
    )
  
  options <- c("PDF", "PNG", "No Export")
  choice <- menu(options, title = "How do you want to save result?")
  
  if (choice == 1) {
    pdf(paste0(output_name, ".pdf"), width = pdf_width, height = pdf_height)
    print(violin_plot)
    dev.off()
    message("PDF Exported")
  } else if (choice == 2) {
    png(paste0(output_name, ".png"), width = png_width, height = png_height, res = png_resolution)
    print(violin_plot)
    dev.off()
    message("PNG Exported")
  } else {
    print(violin_plot)
  }
}