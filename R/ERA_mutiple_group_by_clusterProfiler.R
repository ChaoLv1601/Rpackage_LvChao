#' Perform Multi-group Enrichment Analysis Using clusterProfiler
#'
#' This function performs Gene Ontology (GO) and KEGG enrichment analysis for multiple groups on a given gene list.
#' The results are processed, simplified, and saved in CSV files.
#'
#' @param data_file A character string specifying the path to the input data file containing "gene" and "group" columns.
#' @param species A character string specifying the species, either "human" or "mouse". Default is "human".
#' @param save_path A character string specifying the path where results will be saved. Default is current directory.
#' @param qvalueCutoff A numeric value specifying the cutoff for the q-value. Default is 0.05.
#' @param pAdjustMethod A character string specifying the method to adjust p-values, such as "BH". Default is "BH".
#' @param simplify_cutoff A numeric value specifying the cutoff for simplifying GO terms. Default is 0.7.
#' @param save_raw A logical value indicating whether to save the raw results. Default is TRUE.
#'
#' @return A list containing the processed results for GO categories ("BP", "CC", "MF") and KEGG pathways,
#' and a comparison of the number of terms before and after simplification.
#'
#' @importFrom clusterProfiler compareCluster simplify
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom crayon red
#' @importFrom dplyr mutate
#'
#' @export
ERA_multi_group_by_clusterProfiler <- function(data_file, gene_col = "gene", group_col = "group", species = "human", save_path = ".", qvalueCutoff = 0.05, pAdjustMethod = "BH", simplify_cutoff = 0.7, save_raw = TRUE) {

  # 加载必要包
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(crayon))
# 检查数据框是否提供
  if (!is.data.frame(data)) {
    stop("The input data must be a data frame.")
  }  
  if (species == "human") {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    orgDb <- org.Hs.eg.db
    kegg_species <- "hsa"
  } else if (species == "mouse") {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    orgDb <- org.Mm.eg.db
    kegg_species <- "mmu"
  } else {
    stop("Unsupported species. Please choose 'human' or 'mouse'.")
  }


  # 读取输入数据文件并检查指定的列
  if (!(gene_col %in% colnames(data)) || !(group_col %in% colnames(data))) {
    stop("The specified gene or group column is not found in the data.")
  }  
  
  # 获取分组信息
  groups <- unique(data$group)

  # 打印默认参数值（带颜色）
  message(red("Default Parameters:"))
  message(red("Species: "), species)
  message(red("Save Path: "), save_path)
  message(red("Q-value Cutoff: "), qvalueCutoff)
  message(red("P-value Adjustment Method: "), pAdjustMethod)
  message(red("Simplify Cutoff: "), simplify_cutoff)
  message(red("Gene Column: "), gene_col)
  message(red("Group Column: "), group_col)
  message(red("Groups in Data: "), paste(groups, collapse = ", "))
  
  # 将基因符号转换为 ENTREZID
  gene_groups <- split(data$gene, data$group)
  message("gene_groups is")
  print(str(gene_groups))

  gene_list_entrez <- lapply(gene_groups, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})
  gene_id_groups <- lapply(gene_groups, function(genes) bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgDb))
  message("gene_id_groups is")
  print(str(gene_id_groups))
  # GO富集分析（BP, CC, MF）
  go_terms <- c("BP", "CC", "MF")
  go_results_list <- list()
  
  # 定义数据处理函数
process_go_results <- function(go_results, orgDb) {
  # 提取geneID并转换为基因符号
  entrez_ids <- unlist(strsplit(as.character(go_results$geneID), "/"))
  symbol_map <- mapIds(orgDb, keys = unique(entrez_ids), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

  go_results <- go_results %>%
    mutate(
      SYMBOL = sapply(strsplit(as.character(geneID), "/"), function(x) paste(symbol_map[x], collapse = ",")),
      GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )

  return(go_results)
}


for (ont in go_terms) {
    message("Performing GO enrichment for ", ont, " ...")
    
    # 执行多组富集分析
    go_results <- compareCluster(
      geneCluster = gene_list_entrez,  # 使用转换后的ENTREZID基因集合
      fun = "enrichGO",
      OrgDb = orgDb,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = pAdjustMethod,
      qvalueCutoff = qvalueCutoff
    )
    
  if (!is.null(go_results) && length(go_results@compareClusterResult) > 0) {
    # 获取未简化的GO结果
    go_results_raw_df <- as.data.frame(go_results@compareClusterResult)
    
    # 保存原始的GO结果（未简化）
    if (save_raw) {
      # 对原始结果进行处理
      go_results_raw_df <- process_go_results(go_results_raw_df, orgDb)
      write.csv(go_results_raw_df, file = file.path(save_path, paste0("rawdata_go_results_", ont, ".csv")), row.names = FALSE)
    }
    
    # 简化GO结果
    go_results_simplified <- simplify(go_results, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min)
    go_results_df <- as.data.frame(go_results_simplified@compareClusterResult)
    
    # 保存简化后的GO结果
    write.csv(go_results_df, file = file.path(save_path, paste0("multi_group_go_results_", ont, ".csv")), row.names = FALSE)
    
    # 将结果添加到结果列表中
    go_results_list[[ont]] <- go_results_df
  } else {
    message("No significant results found for ", ont, ". Skipping simplification.")
  }
  }

  # KEGG富集分析
  message("Performing KEGG enrichment...")
  kegg_results <- compareCluster(
    geneCluster = gene_id_groups,
    fun = "enrichKEGG",
    organism = kegg_species,
    pAdjustMethod = pAdjustMethod,
    qvalueCutoff = qvalueCutoff
  )

  # 保存KEGG结果
  kegg_results_df <- as.data.frame(kegg_results)
  write.csv(kegg_results_df, file = file.path(save_path, "multi_group_kegg_results.csv"), row.names = FALSE)

  return(list(GO_results = go_results_list, KEGG_results = kegg_results_df))
}

