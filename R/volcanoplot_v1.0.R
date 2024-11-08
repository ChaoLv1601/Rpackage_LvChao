#' Volcano Plot Function
#'
#' Generates a volcano plot with labeled genes.
#'
#' This function takes a data frame containing gene names, log2 fold change, and p-values,
#' and generates a volcano plot. It labels a specified number of upregulated and 
#' downregulated genes.
#'
#' @param data Data frame containing gene names, log2 fold change, and p-values.
#' @param logFC_col Column name for log2 fold change (default is "log2FC").
#' @param p_value_col Column name for p-values (default is "p_value").
#' @param gene_col Column name for gene names (default is "Gene").
#' @param logFC_threshold Threshold for log2 fold change (default is 1).
#' @param p_value_threshold Threshold for p-value (default is 0.05).
#' @param n_labels Number of up/down genes to label (default is 5).
#' @param seed Random seed for reproducibility (default is 2024).
#' @return A ggplot object representing the volcano plot.
#' @import ggplot2
#' @import ggrepel
#' @export
#' @examples
#' # Assuming 'volcano_data' is a data frame with Gene, log2FC, and p_value columns
#' volcanoplot_v1.0(volcano_data)

volcanoplot_v1.0 <- function(data, logFC_col = "log2FC", p_value_col = "p_value", 
                             gene_col = "Gene", logFC_threshold = 1, 
                             p_value_threshold = 0.05, n_labels = 5, seed = 2024) {
  
  # 设置随机种子
  set.seed(seed)
  
  # 判断输入的数据框是否包含指定的列
  if(!all(c(logFC_col, p_value_col, gene_col) %in% colnames(data))) {
    stop("The data frame must contain the specified columns for logFC, p-value, and gene names.")
  }
  
  # 计算基因的状态（上调、下调、未变）
  data$change <- ifelse(data[[logFC_col]] > logFC_threshold & data[[p_value_col]] < p_value_threshold, 
                        "up", 
                        ifelse(data[[logFC_col]] < -logFC_threshold & data[[p_value_col]] < p_value_threshold, 
                               "down", 
                               "unchange"))
  
  # 随机选择要标注的上调和下调基因
  up_genes <- sample(data[[gene_col]][data$change == "up"], n_labels)
  down_genes <- sample(data[[gene_col]][data$change == "down"], n_labels)
  
  # 创建火山图
  p <- ggplot(data, aes_string(x = logFC_col, y = paste0("-log10(", p_value_col, ")"))) + 
    geom_point(alpha = 0.4, size = 3.5, aes(color = change)) + 
    ylab("-log10(Pvalue)") + 
    scale_color_manual(values = c("blue4", "grey", "red3")) +  
    theme_bw() + 
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = 4, color = "black", lwd = 0.8) + 
    geom_hline(yintercept = -log10(p_value_threshold), linetype = 4, color = "black", lwd = 0.8)
  
  # 添加标签
  p <- p + geom_text_repel(data = subset(data, data[[gene_col]] %in% c(up_genes, down_genes)), 
                           aes_string(label = gene_col), size = 3, 
                           box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50')
  
  return(p)
}

