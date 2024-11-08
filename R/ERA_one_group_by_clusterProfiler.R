#' Perform Enrichment Analysis Using clusterProfiler
#'
#' This function performs Gene Ontology (GO) and KEGG enrichment analysis on a given gene list. 
#' It supports both human and mouse species. The results are processed, simplified, and saved in CSV files.
#'
#' @param gene_list A character vector of gene symbols.
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
#' @details 
#' This function performs enrichment analysis for biological process (BP), cellular component (CC), 
#' and molecular function (MF) categories in GO, as well as KEGG pathway enrichment. The results 
#' are simplified by removing redundant terms, and both raw and simplified results are saved in CSV files.
#' 
#' @examples
#' \dontrun{
#' gene_list <- c("TP53", "BRCA1", "EGFR")
#' result <- ERA_one_group_by_clusterProfiler(gene_list = gene_list, species = "human")
#' }
#' 
#' @importFrom clusterProfiler enrichGO enrichKEGG simplify
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom dplyr mutate
#' 
#' @export
ERA_one_group_by_clusterProfiler <- function(gene_list, species = "human", save_path = ".", qvalueCutoff = 0.05, pAdjustMethod = "BH", simplify_cutoff = 0.7, save_raw = TRUE) {
    
  # 打印默认参数值
  message("Default Parameters:")
  message("Species: ", species)
  message("Save Path: ", save_path)
  message("Q-value Cutoff: ", qvalueCutoff)
  message("P-value Adjustment Method: ", pAdjustMethod)
  message("Simplify Cutoff: ", simplify_cutoff)  

  # 加载必要包
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(dplyr))
  
  # 根据物种选择数据库
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
  
  # 将基因符号转换为 ENTREZID
  gene_id <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgDb)
  
  # 定义 GO 术语
  go_terms <- c("BP", "CC", "MF")
  results_list <- list()
  comparison_df <- data.frame(Term = character(), Before = integer(), Significant = integer(), After = integer(), stringsAsFactors = FALSE)
  
  # 帮助函数用于处理富集分析结果
  process_go_results <- function(go_results) {
    entrez_ids <- unlist(strsplit(as.character(go_results$geneID), "/"))
    symbol_map <- mapIds(orgDb, keys = unique(entrez_ids), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    go_results <- go_results %>%
      mutate(
        SYMBOL = sapply(strsplit(as.character(geneID), "/"), function(x) paste(symbol_map[x], collapse = ",")),
        GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
      )
    return(go_results)
  }
 
  # Perform GO enrichment analysis
for (term in go_terms) {
  message("Performing ", term, " enrichment...")
  go_results <- enrichGO(gene = gene_id$ENTREZID, OrgDb = orgDb, keyType = "ENTREZID", ont = term, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)

  # Total results count
  before_count <- nrow(go_results@result)

  # Count results with padj < 0.05 (Significant entries)
  significant_count <- sum(go_results@result$p.adjust < 0.05, na.rm = TRUE)

  go_results_raw_df <- as.data.frame(go_results@result)

  # Process raw results
  go_results_raw_df <- process_go_results(go_results_raw_df)

  if (save_raw) {
    write.csv(go_results_raw_df, file = file.path(save_path, paste0("rawdata_go_results_", term, ".csv")), row.names = FALSE)
  }

  # Simplify GO terms
  go_results <- simplify(go_results, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min)
  after_count <- nrow(go_results@result)

  # Print before_count, significant_count, and after_count for each iteration
  message("Term: ", term)
  message("Before count: ", before_count)
  message("Significant count (padj < 0.05): ", significant_count)
  message("After count (after simplification): ", after_count)

  # Add counts to comparison_df
  comparison_df <- rbind(comparison_df, data.frame(Term = term, Before = before_count, Significant = significant_count, After = after_count))

  # Print comparison_df in each iteration
  print(comparison_df)


  go_results_df <- process_go_results(as.data.frame(go_results@result))
  go_results_df$Category <- term
  results_list[[term]] <- go_results_df
}

# Save comparison_df to CSV
write.csv(comparison_df, file = file.path(save_path, "go_results_comparison.csv"), row.names = FALSE)

# Perform KEGG enrichment analysis
message("Performing KEGG enrichment...")
kegg_results <- enrichKEGG(gene = gene_id$ENTREZID, organism = kegg_species, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff)

if (!is.null(kegg_results) && nrow(kegg_results@result) > 0) {
  kegg_results_df <- as.data.frame(kegg_results@result)
  
  # Process KEGG results (using the same function as for GO results)
  kegg_results_df <- process_go_results(kegg_results_df)
  
  # Count significant entries
  significant_count_kegg <- sum(kegg_results_df$p.adjust < 0.05, na.rm = TRUE)
  
  # Save raw KEGG results
  if (save_raw) {
    write.csv(kegg_results_df, file = file.path(save_path, "rawdata_kegg_results.csv"), row.names = FALSE)
  }

  # Update comparison_df
  comparison_df <- rbind(comparison_df, data.frame(Term = "KEGG", Before = nrow(kegg_results_df), Significant = significant_count_kegg, After = nrow(kegg_results_df)))
  
  print(comparison_df)

  results_list[["KEGG"]] <- kegg_results_df
} else {
  message("No significant KEGG pathways found.")
}

# 保存最终结果
for (term in c(go_terms, "KEGG")) {
  if (!is.null(results_list[[term]])) {
    write.csv(results_list[[term]], file = file.path(save_path, paste0("results_", term, ".csv")), row.names = FALSE)
  }
}

return(list(results = results_list, comparison = comparison_df))

}
