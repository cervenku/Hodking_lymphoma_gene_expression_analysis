#' Calculate hiearchical clustering and plot dendrogram.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   method_distance: Distance metric used for clustering. See ?dist.
#'                    Can be also correlation metrix ("pearson", "kendall", "spearman"). See ?cor.
#'   method_clustering: Clustering method. See ?hclust.
#'   color_by: Vector of discrete values to color samples by.
#'             Length must be the same as number of columns in 'm'.
#'   color_by_lab: Name of the color legend.
#'   title: Main plot title.
#'
#' Returns:
#'   dendrogram object
plot_hc <- function(
  m,
  method_distance = "euclidean",
  method_clustering = "complete",
  color_by = NULL,
  color_by_lab = "Group",
  title = "Hierarchical Clustering"
) {
  m2 <- data.frame(t(m))
  if(method_distance %in% c("pearson", "kendall", "spearman")){
    d <- as.dist(1-cor(t(m2),method=method_distance))
  }
  else{
    d <- dist(m2, method = method_distance)
  }
  
  dend <- as.dendrogram(hclust(d, method_clustering))
  
  cols <- unique(as.vector(color_by))
  cols_vec <- rep(NA, length(as.vector(color_by)))
  j<-1
  for (i in as.vector(color_by)) {
    cols_vec[j] <- match(i, cols) + 1
    j <- j + 1
  }
  labels_order <- order.dendrogram(dend)
  ordered_cols <- rep(NA, length(as.vector(color_by)))
  j<- 1
  for (i in as.vector(labels_order)) {
    ordered_cols[j] <- cols_vec[i]
    j <- j + 1
  }
  
  dend <- dend %>% set("labels_col", ordered_cols)
  plot(dend, main = title)
  legend(
    "topright", title = color_by_lab, legend = cols,
    pch = 19, col = 2:(length(cols)+1), bty = "n"
  )
}

#' Calculate hiearchical clustering and plot dendrogram using the dendextend package.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   method_distance: Distance metric used for clustering. See ?dist.
#'                    Can be also correlation metrix ("pearson", "kendall", "spearman"). See ?cor.
#'   method_clustering: Clustering method. See ?hclust.
#'   color_by: Vector of discrete values to color samples by.
#'             Length must be the same as number of columns in 'm'.
#'   color_by_lab: Name of the color legend.
#'   title: Main plot title.
#'
#' Returns:
#'   dendrogram object
plot_hc2 <- function(
  m,
  method_distance = "euclidean",
  method_clustering = "complete",
  color_by = NULL,
  color_by_lab = "Group",
  title = "Hierarchical Clustering"
) {
}

#' Select features (rows) with the highest variance across the samples (columns).
#' This will be useful for visualization (mainly heatmaps) of large expression matrices.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   n_top_features: Number of top features.
#'
#' Returns:
#'   Subset of 'm' with 'n_top_features' with the highest variance across the samples.
select_var_features <- function(m, n_top_features) {
  variances <- apply(X=m, MARGIN=1, FUN=var)
  sorted <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[1:n_top_features]
  return(m[sorted, ])
}

#' Using the ggplot2 package, plot the first PC or first to three PCs of samples in expression matrix.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   sample_data: Dataframe describing samples.
#'   plot_type: "single" for PC1 vs. PC2, "multi" for combinations of PC1-3 and their cumulative explained variance.
#'   n_top_features: Number of top features with the highest variance across the samples.
#'   color_by: Column name in sample_data to use for point coloring.
#'   shape_by: Column noneame in sample_data to use for point shape.
#'   label_by: Column name in sample_data to use for point labels.
#'   point_size: Point size (numeric).
#'   text_size: Label text size (numeric).
#'   center: Whether to center PCA. See ?prcomp.
#'   scale.: Whether to scale PCA. See ?prcomp.
#'
#' Returns:
#'   list(pca = <prcomp object>, pca_df = <combined dataframe of sample_data and PCs>, plot = <ggplot2 or patchwork object (depends on plot_type)>)
plot_pca <- function(
  m,
  sample_data,
  plot_type = c("single", "multi"),
  n_top_features = Inf,
  color_by = NULL,
  shape_by = NULL,
  label_by = NULL,
  point_size = 2,
  text_size = 2.5,
  center = TRUE,
  scale. = TRUE
) {
  
  if(!is.infinite(n_top_features)){
    m <- select_var_features(m, n_top_features)
  }
  m_pca <- prcomp(as.data.frame(t(m)),scale = scale., center = center)
  df_pca <- as.data.frame(m_pca$x)
  df_combined <- cbind(sample_data, df_pca)
  var <- as.vector(summary(m_pca)$importance[2,])
  var <- ceiling(var*100)
  
  if("single" %in% plot_type){
    p <- ggplot(df_combined, aes(x = PC1, y = PC2)) +
      geom_point(aes_string(shape=shape_by, color=color_by), size=point_size) + xlab(paste("PC1 (", var[1], "%)", sep ="")) + 
      ylab(paste("PC2 (", var[2], "%)", sep="")) 
    if(!is.null(label_by)){
      p <- p +geom_text(aes_string(label=label_by),hjust=0, vjust=0, size=text_size)
    }
  }
  else{
    cum_var <- as.vector(summary(m_pca)$importance[3,])
    df_cum <- data.frame(PC = 1:length(cum_var), Cum_var = cum_var)
    p1 <- ggplot(df_combined, aes(x = PC1, y = PC2)) +
      geom_point(aes_string(shape=shape_by, color=color_by), size=point_size) + xlab(paste("PC1 (", var[1], "%)", sep ="")) + 
      ylab(paste("PC2 (", var[2], "%)", sep=""))
    p2 <- ggplot(df_combined, aes(x = PC1, y = PC3)) +
      geom_point(aes_string(shape=shape_by, color=color_by), size=point_size) + xlab(paste("PC1 (", var[1], "%)", sep ="")) + 
      ylab(paste("PC3 (", var[3], "%)", sep=""))
    p3 <- ggplot(df_combined, aes(x = PC2, y = PC3)) +
      geom_point(aes_string(shape=shape_by, color=color_by), size=point_size) + xlab(paste("PC2 (", var[2], "%)", sep ="")) + 
      ylab(paste("PC3 (", var[3], "%)", sep=""))
    p4 <- ggplot(data=df_cum, aes(x=PC, y=cum_var)) +
      geom_bar(stat="identity") + ylab("Cumulative % of var. explained")
    p <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
    
  }
  res <- list("pca" = m_pca, "pca_df" = df_combined, "plot" = p)
  return(res)
}

#' Using the GGally::ggpairs() function, plot grid of PCA plots (PC1 vs. PC2, PC1 vs. PC3, PC2 vs. PC3, etc).
#' When n_components == 2, use normal ggplot2.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   sample_data: Dataframe describing samples.
#'   output_file: File to save plot in.
#'   n_components: Number of PCs to plot.
#'   n_top_features: Number of top features with the highest variance across the samples.
#'   color_by: Column name in sample_data to use for point coloring.
#'   color_legend_lab: Name of the color legend.
#'   shape_by: Column name in sample_data to use for point shape.
#'   shape_legend_lab: Name of the shape legend.
#'   label_by: Column name in sample_data to use for point labels.
#'   point_size: Point size (numeric).
#'   text_size: Label text size (numeric).
#'   title: Plot title.
#'   subtitle: Plot subtitle.
#'   center: Whether to center PCA. See ?prcomp.
#'   scale.: Whether to scale PCA. See ?prcomp.
#'
#' Returns:
#'   list(pca = <prcomp object>, pca_df = <combined dataframe of sample_data and PCs>, plot = <ggplot2 object>)
plot_pca_ggpairs <- function(
  m,
  sample_data,
  n_components = 5,
  n_top_features = Inf,
  color_by = NULL,
  color_legend_lab = NULL,
  shape_by = NULL,
  shape_legend_lab = NULL,
  label_by = NULL,
  point_size = 2,
  text_size = 2.5,
  title = NULL,
  subtitle = NULL,
  center = TRUE,
  scale. = TRUE
) {
  m_pca <- prcomp(as.data.frame(t(m)),scale = scale., center = center)
  df_pca <- as.data.frame(m_pca$x)
  df_combined <- cbind(df_pca, sample_data)
  
  
  if(n_components == 2){
    p <- plot_pca(m,sample_data, plot_type = "single",
                  n_top_features = n_top_features,
                  color_by = color_by,
                  shape_by = shape_by,
                  label_by = label_by,
                  point_size = point_size,
                  text_size = text_size,
                  center = center,
                  scale. = scale
    )$plot
  }
  else{
    p <- GGally::ggpairs(df_combined, columns = 1:n_components,
                         ggplot2::aes_string(colour=color_by, shape=shape_by),
                         upper = list(continuous = "points", discrete = "points", na = "points"),
                         diag=list(continuous = "blankDiag", discrete = "blankDiag", na = "blankDiag"),
                         legend = 3, progress=FALSE)
  }
  return(list(pca = m_pca, pca_df = df_combined, plot = p + labs(title = title, subtitle = subtitle)))
  
}

#' Plot heatmap using the ComplexHeatmap package.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   z_score: If TRUE, calculate row z-score.
#'   column_annotation: Dataframe used for annotation of columns.
#'   row_annotation: Dataframe used for annotation of rows.
#'   title: Heatmap title.
#'   legend_title: Heatmap color legend title.
#'   show_row_names: If TRUE, show rownames in the heatmap.
#'   show_col_names: If TRUE, show colnames in the heatmap.
#'   color_palette: Function to generate colors for annotations.
#'   color_mapping: Named list of named vectors to map colors to variable levels.
#'                  See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#simple-annotation
#'
#' Returns:
#'   ComplexHeatmap object
plot_heatmap <- function(
  m,
  n_top_features = Inf,
  z_score = FALSE,
  column_annotation = NULL,
  row_annotation = NULL,
  title = "",
  legend_title = "Values",
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  color_palette = scales::hue_pal(),
  color_mapping = NULL
) {
  if(!is.infinite(n_top_features)){
    m <- select_var_features(m, n_top_features)
  }
  #if(z_score){
   # m <- scale(as.matrix(m))
  #}
  if(!(is.null(row_annotation))){
    m <- as.data.frame(m)
    m$rowNames <- row.names(m)
    m <- as.data.frame(m) %>% dplyr::filter(., rowNames %in% row.names(as.data.frame(row_annotation)))
    m <- subset(m, select = -rowNames)
    p_heatmap <- Heatmap(
      as.matrix(m),
      name = legend_title,
      column_title = title,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      top_annotation = HeatmapAnnotation(df = as.data.frame(column_annotation)),
      right_annotation = rowAnnotation(df = as.data.frame(row_annotation))
    )
  }
  else{
    p_heatmap <- Heatmap(
      as.matrix(m),
      name = legend_title,
      column_title = title,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      top_annotation = HeatmapAnnotation(df = column_annotation)
    )
  }
  
  return(p_heatmap)
  
  
}

#' Create a heatmap using the heatmaply package.
#'
#' Args:
#'   m: Expression matrix (rows are features, columns are samples).
#'   z_score: If TRUE, calculate row z-score.
#'   column_annotation: Dataframe used for annotation of columns.
#'   row_annotation: Dataframe used for annotation of rows.
#'   title: Heatmap title.
#'   legend_title: Heatmap color legend title.
#'
#' Returns:
#'   plotly object
plot_heatmaply <- function(
  m,
  z_score = FALSE,
  column_annotation = NULL,
  row_annotation = NULL,
  main = NULL,
  legend_title = NULL,
  showticklabels = c(TRUE, TRUE)
) {
  m <- as.data.frame(m)
  m$rowNames <- row.names(m)
  m <- as.data.frame(m) %>% dplyr::filter(., rowNames %in% row.names(as.data.frame(row_annotation)))
  m <- subset(m, select = -rowNames)
  m <- m %>% dplyr::select(., row.names(row_annotation))
  
  if(is.null(column_annotation)){
    p <- heatmaply_cor(
      m,
      main = main, 
      showticklabels = showticklabels,
      key.title = legend_title,
      row_side_colors = row_annotation
    )
    
  }
  else{
    p <- heatmaply_cor(
      m,
      main = main, 
      showticklabels = showticklabels,
      key.title = legend_title,
      row_side_colors = row_annotation,
      col_side_colors = column_annotation
    )
  }
  
  return(p)
  
}

#' Using the ggpubr::ggboxplot() function, plot boxplots of gene expression.
#'
#' Args:
#'   plot_data: data.frame (long format)
#'   x: Column to divide x-axis values to (e.g. sample groups).
#'   y: Column to compute boxplots on y-axis.
#'   facet_by: One or two columns used for facetting.
#'   feature: Name of a feature from which boxplots will be made.
#'            Data will be filtered based on facet_by.
#'            E.g. if facet_by = "gene" and feature = "CD24", only boxplots for "CD24" will be made.
#'   color_by: Column to use for boxplot and point coloring.
#'   x_lab: Name of x-axe.
#'   y_lab: Name of y-axe.
#'   main: Main plot title.
#'   add: Add something more to boxplots.
#'        Allowed values are one or the combination of:
#'        "none", "dotplot", "jitter", "boxplot", "point", "mean",
#'        "mean_se", "mean_sd", "mean_ci", "mean_range", "median",
#'        "median_iqr", "median_mad", "median_range".
#'        See ?ggpubr::ggboxplot
#'   point_size: Size of points inside boxplots.
#'   outlier_shape: Which point shape to use for outliers.
#'   do_t_test: Whether to do the t-test and display p-values inside the plot.
#'
#' Returns:
#'   ggplot2 object
plot_boxplots <- function(
  plot_data,
  x,
  y,
  facet_by,
  feature = NULL,
  color_by = NULL,
  x_lab = x,
  y_lab = y,
  main = NULL,
  add = "jitter",
  point_size = 2,
  outlier_shape = 0,
  do_t_test = TRUE
) {
  if(!is.null(feature)){
    f <- feature
    plot_data <- plot_data %>% dplyr::filter(!!sym(facet_by) == f)
  }
  p<-ggplot(plot_data, aes_string(x=x, y=y, color=color_by)) +
    geom_boxplot(outlier.shape = outlier_shape) + geom_jitter(aes(size = point_size)) + facet_wrap(reformulate(facet_by,"."), nrow = 2, ncol = 2) +
    xlab(x_lab) + ylab(y_lab) + ggtitle(main) + theme(legend.position="top")
  
  if(do_t_test){
    p <- p + stat_compare_means(method = "t.test")
  }
  
  return(p)
  
}

#' Compute the M value of CP values.
#'
#' Args:
#'  gene: Name of gene to compute the M value for.
#'  cp: Matrix or dataframe of CP values. Rows are genes and columns are samples.
#'
#' Returns:
#'  M value (numeric)
compute_m <- function(gene, cp) {
  # Ajk = {(aji/aki)}
  # Vjk = st.dev(Ajk)
  # Mj = sum(Vjk)/(n-1)
  
  others <- rownames(cp)[!(rownames(cp) == gene)]
  a_jk <- rep(NA, ncol(cp))
  v_j <- rep(NA, length(others))
  v_index <- 1
  
  for(k in others){
    in_A <- cp[c(gene, k),]
    index <- 1
    for(i in colnames(in_A)){
      a_jk[index] <- (in_A[gene, i]) - (in_A[k, i])
      index <- index + 1
    }
    v_j[v_index] <- sd(a_jk)
    v_index <- v_index + 1
  }
  m_j <- sum(v_j)/length(others)
  return(m_j)
  
}

#' For a single gene, test for statistical significance of difference in group means.
#'
#' Args:
#'   gene: Name of gene to test.
#'   gene_data: Dataframe in long format.
#'   gene_col: Column with genes.
#'   value_col: Column with values to test.
#'   group_col: Column with sample groups. There must be exactly two different groups.
#'   test: Statistical test to perform. It must have the same interface as t.test()
#'   verbose: Whether to print test results.
#'
#' Returns:
#'   htest object
test_gene <- function(gene, gene_data, gene_col, value_col, group_col, test = t.test, verbose = TRUE) {
  g <- gene
  gene_df <- gene_data %>% dplyr::filter(!!sym(gene_col) == g)
  group_col
  value_col
  
  exp1 <- expr(!!ensym(value_col) ~ !!ensym(group_col))
  
  test_res <- test(formula = eval(exp1), data = gene_df)
  if(verbose){
    print(test_res)
  }
  return(test_res)
}

#' For all genes in the input dataframe, test for statistical significance of difference in group means.
#'
#' Args:
#'   gene_data: Dataframe in long format.
#'   gene_col: Column with genes.
#'   value_col: Column with values to test.
#'   group_col: Column with sample groups. There must be exactly two different groups.
#'   test: Statistical test to perform. It must have the same interface as t.test()
#'
#' Returns:
#'   tibble object
test_gene_table <- function(gene_data, gene_col, value_col, group_col, test = t.test) {
  genes <- unique(as.matrix(gene_data)[,gene_col])
  htest_list <- lapply(genes, function(x) test_gene(x, gene_data, gene_col, value_col, group_col, test, verbose = FALSE))
  
  p_val <- sapply(1:length(htest_list), function(x) htest_list[[x]]$p.value)
  
  df <- data.frame(gene = genes, p_value = p_val)
  
  df <- dplyr::mutate(
    df,
    significance = case_when(
      is.na(p_value) ~ NA_character_,
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    ),
    test_res = htest_list
  )
  
  return(as_tibble(df))
}

#' Return asterisks according to p-values.
#'
#' Args:
#'   p_value: Vector of p-values.
#'
#' Returns:
#'  character vector
asterisk <- function(p_value) {
}
