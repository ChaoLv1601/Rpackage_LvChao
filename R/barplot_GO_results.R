#' Plot Gene Ontology Results in a Barplot
#'
#' This function generates a bar plot to visualize Gene Ontology (GO) enrichment results.
#' The bars are sorted by the log10-transformed p-values, and gene symbols are added as labels.
#' Only the top 15 most significant GO terms are displayed by default.
#' If any GO term has more than 10 associated genes (symbols), all symbols will not be displayed.
#'
#' @param data A data frame containing GO results.
#' @param description_col A string for the name of the column containing GO term descriptions.
#' @param logP_col A string for the name of the column containing log10-transformed p-values. If provided, padj will not be used.
#' @param padj_col A string for the name of the column containing adjusted p-values.
#' @param symbol_col A string for the name of the column containing gene symbols.
#' @param title A character string for the plot title. Default is "Gene Annotation Results".
#' @param title A character string for the plot title. Default is "Gene Annotation Results".
#' @param fill_low The color for the low end of the gradient fill. Default is "lightyellow".
#' @param fill_high The color for the high end of the gradient fill. Default is "darkorange".
#' @param top_n An integer specifying how many top GO terms to display. Default is 15.
#'
#' @return A list containing:
#'   - `plot`: A ggplot2 object representing the GO enrichment barplot.
#'   - `data`: A data frame with the columns Description, Symbol, and LogP.
#'
#' @import ggplot2
#' @import crayon
#' @export
#'
#' @examples
#' # Example usage:
#' # result <- barplot_GO_results(Description = data$Description, LogP = NULL, Symbol = data$SYMBOL, padj = data$p.adjust)
#' # View(result$data)  # To view the output data frame
#'
barplot_GO_results <- function(data, description_col, padj_col, symbol_col, logP_col = NULL,
                               title = "Gene Annotation Results", 
                               fill_low = "lightyellow", fill_high = "darkorange", 
                               top_n = 10) {
  
  # Extract relevant columns from data
  Description <- data[[description_col]]
  Symbol <- data[[symbol_col]]

  # Determine LogP value
  if (!is.null(logP_col)) {
    LogP <- data[[logP_col]]  
    message(red("Using LogP directly.")) # Use provided LogP values directly
  } else if (!is.null(padj_col)) {
    padj <- data[[padj_col]]
    LogP <- -log10(padj)  # Calculate LogP from padj
    message(red("Using adjusted p-values (padj) to calculate LogP."))
  } else {
    stop("Either logP_col or padj_col must be provided.")
  }

  # Create a data frame to hold the provided data
  data <- data.frame(Description = Description, LogP = LogP, Symbol = Symbol)
  
  # Count the number of associated genes for each Symbol
  gene_counts <- sapply(strsplit(as.character(data$Symbol), ","), length) 
  
  # Check if any GO term has more than 10 associated genes
  if (any(gene_counts > 10)) {
    data$Symbol <- ""  # Set all Symbols to "" if any gene count exceeds 10
    message(red("Some GO terms have more than 10 associated genes. Symbols will not be displayed."))
  }

  # Sort the data by LogP and select the top_n terms
  data <- data[order(data$LogP, decreasing = T), ][1:min(top_n, nrow(data)), ]
  
  cat(blue(paste(capture.output(print(data)), collapse = "\n")))

  # Generate the bar plot
  ggplot(data, aes(x = reorder(Description, LogP, decreasing = F), y = LogP, fill = LogP)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_gradient(low = fill_low, high = fill_high) +
    geom_text(aes(x = as.numeric(reorder(Description, LogP, decreasing = F)) - 0.4, y = 0.1, label = Symbol), 
              hjust = 0, size = 3.5) +
    coord_flip() +
    geom_text(aes(x = reorder(Description, LogP, decreasing = F), y = 0.1, label = Description), 
              hjust = 0, size = 3.5, color = "black") +
    labs(x = NULL, y = expression(log[10](P)), title = title) +
    theme_classic() +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          plot.title = element_text(hjust = 0.5, size = 20), 
          axis.text.x = element_text(size = 15), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.box.background = element_blank(), 
          legend.background = element_blank())
}

