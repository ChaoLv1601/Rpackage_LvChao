#' View Matrix Subset
#'
#' Extracts a subset of a matrix, including the first 5 rows and specific columns.
#'
#' This function extracts a subset of a matrix. If the number of columns in the matrix 
#' is less than 6, all columns are included. Otherwise, the function selects the first 
#' 3 columns and the last 3 columns. The function returns the first 5 rows of these 
#' selected columns.
#'
#' @param mat A matrix from which to extract a subset.
#' @return A matrix containing the first 5 rows and selected columns.
#' @examples
#' # Create a sample matrix
#' example_matrix <- matrix(1:20, nrow = 5, ncol = 4)
#' # View subset of the matrix
#' result <- view_matrix_subset(example_matrix)
#' print(result)
#' @export

view_matrix_subset <- function(mat) {
  # 获取矩阵的列数
  num_cols <- ncol(mat)
  
  # 判断列数是否小于6
  if (num_cols < 6) {
    selected_cols <- 1:num_cols  # 选择所有列
  } else {
    # 提取前3列和后3列
    selected_cols <- c(1:3, (num_cols-2):num_cols)
  }
  
  # 提取前5行
  selected_data <- mat[1:5, selected_cols]
  
  # 返回结果
  return(selected_data)
}

# 示例用法
# 假设 your_matrix 是你要查看的矩阵
# result <- view_matrix_subset(your_matrix)
# print(result)

