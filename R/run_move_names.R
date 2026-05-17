#' Move Row or Column Names into the Data Area
#'
#' @description
#' This function moves existing row names into a specified column or column names 
#' into a specified row. By default, it prioritizes moving row names to the first 
#' column. If both a column index and a row index are provided, row name 
#' movement takes precedence and a warning is issued.
#'
#' @param x A dataframe, tibble, or matrix.
#' @param col_index Numeric. The column index where to place the row names. Defaults to 1.
#' @param row_index Numeric or NULL. The row index where to place the column names. 
#'   If not NULL and `col_index` is also active, this will be ignored with a warning.
#' @param col_name String. The name of the new column created when moving row names. 
#'   Defaults to "row_names".
#' @param remove Logical. If TRUE (default), removes the original row/column names 
#'   after moving them into the data.
#'
#' @return A data frame with the names moved into the data area. When moving row 
#'   names to a column, original data types of other columns are preserved.
#'
#' @importFrom methods is
#' 
#' @author John Lennon L. Calorio
#' 
#' @export
#'
#' @examples
#' # 1. Move row names to the first column (default)
#' run_move_names(mtcars[1:5, 1:3])
#'
#' # 2. Move row names to a specific column index
#' run_move_names(mtcars[1:5, 1:3], col_index = 2)
#'
#' # 3. Move column names to the first row
#' # Note: col_index must be NULL to process row_index per requirements
#' run_move_names(mtcars[1:5, 1:3], col_index = NULL, row_index = 1)
#'
#' # 4. Demonstrate priority warning
#' run_move_names(mtcars[1:5, 1:3], col_index = 1, row_index = 1)
run_move_names <- function(
  x, 
  col_index = 1, 
  row_index = NULL, 
  col_name = "row_names",
  remove = TRUE
) {
  
  # --- 1. Validation ---
  if (!is.data.frame(x) && !is.matrix(x)) {
    stop("Input 'x' must be a dataframe, tibble, or matrix. Received: ", class(x)[1])
  }
  
  # --- 2. Conflict Handling & Priority Logic ---
  # If both are provided, prioritize col_index (row names)
  if (!is.null(col_index) && !is.null(row_index)) {
    warning("Ambiguity detected: Both 'col_index' and 'row_index' provided. ",
            "Prioritizing movement of row names to column. 'row_index' will be ignored.")
    row_index <- NULL
  }

  # --- 3. Move Row Names to Column ---
  if (!is.null(col_index)) {
    # Extract current row names
    rn <- rownames(x)
    
    # Check if row names are just default numeric indices
    if (is.null(rn) || all(rn == as.character(seq_len(nrow(x))))) {
      warning("The input does not appear to have explicit row names to move.")
    }

    res <- as.data.frame(x)
    
    # Bound check for col_index
    target_pos <- as.integer(col_index)
    if (target_pos < 1 || target_pos > (ncol(res) + 1)) {
      stop("Invalid 'col_index' (", target_pos, "). Must be between 1 and ", ncol(res) + 1)
    }

    # Create the new column as a data frame to preserve names and types during cbind
    new_col <- stats::setNames(data.frame(rn, stringsAsFactors = FALSE), col_name)

    # Insert row names at target column index
    if (target_pos == 1) {
      res <- cbind(new_col, res)
    } else if (target_pos > ncol(res)) {
      res <- cbind(res, new_col)
    } else {
      res <- cbind(res[, 1:(target_pos - 1), drop = FALSE],
                   new_col,
                   res[, target_pos:ncol(res), drop = FALSE])
    }

    if (remove) {
      rownames(res) <- NULL
    }
    
    return(res)
  }

  # --- 4. Move Column Names to Row ---
  if (!is.null(row_index)) {
    cn <- colnames(x)
    
    if (is.null(cn)) {
      stop("The input does not have column names to move.")
    }

    # Moving headers to a row usually requires character conversion to avoid type coercion issues
    res_mat <- as.matrix(x)
    target_pos <- as.integer(row_index)
    
    if (target_pos < 1 || target_pos > (nrow(res_mat) + 1)) {
      stop("Invalid 'row_index' (", target_pos, "). Must be between 1 and ", nrow(res_mat) + 1)
    }

    # Insert column names as a new row
    if (target_pos == 1) {
      res_mat <- rbind(cn, res_mat)
    } else if (target_pos > nrow(res_mat)) {
      res_mat <- rbind(res_mat, cn)
    } else {
      res_mat <- rbind(res_mat[1:(target_pos - 1), , drop = FALSE],
                       cn,
                       res_mat[target_pos:nrow(res_mat), , drop = FALSE])
    }

    res <- as.data.frame(res_mat, stringsAsFactors = FALSE)

    if (remove) {
      # Set generic column names (V1, V2...) to prevent duplicate/confusing headers
      colnames(res) <- paste0("V", seq_len(ncol(res)))
    }
    
    return(res)
  }

  return(x)
}
