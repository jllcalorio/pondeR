#' Recode and Factorize a Column
#'
#' This function recodes values in a specified column of a data frame according
#' to a provided mapping, then converts the column to a factor with ordered levels.
#' The recoded values appear first in the factor levels, followed by any remaining
#' unique values in alphabetical order.
#'
#' @param df A data frame containing the column to be recoded and factorized.
#' @param column A character string specifying the name of the column to recode.
#'   The column must exist in the data frame.
#' @param conversion_map A named list or vector where names are the original values
#'   to be replaced and values are the new values. All elements must be atomic
#'   (character, numeric, or logical).
#'
#' @return A data frame identical to the input except that the specified column
#'   has been recoded according to the conversion map and converted to a factor
#'   with ordered levels.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Converts the specified column to character to handle mixed data types
#'   \item Applies the recoding based on the conversion map
#'   \item Creates factor levels with recoded values first, followed by other
#'         unique values in alphabetical order
#'   \item Converts the column to a factor with the specified level ordering
#' }
#'
#' @examples
#' # Basic usage with character recoding
#' df <- data.frame(
#'   id = 1:5,
#'   status = c("A", "B", "C", "A", "D")
#' )
#' map <- c("A" = "Active", "B" = "Inactive", "C" = "Pending")
#' result <- recode_and_factorize(df, "status", map)
#' levels(result$status)  # "Active" "Inactive" "Pending" "D"
#'
#' # Usage with numeric to character conversion
#' df2 <- data.frame(
#'   score = c(1, 2, 3, 1, 4),
#'   grade = c(90, 85, 78, 92, 65)
#' )
#' grade_map <- c("1" = "Excellent", "2" = "Good", "3" = "Fair")
#' result2 <- recode_and_factorize(df2, "score", grade_map)
#'
#' @export
recode_and_factorize <- function(df, column, conversion_map) {
  # Input validation
  if (!is.data.frame(df)) {
    stop("Argument 'df' must be a data frame.")
  }

  if (!is.character(column) || length(column) != 1) {
    stop("Argument 'column' must be a single character string.")
  }

  if (!column %in% names(df)) {
    stop("Column '", column, "' not found in the data frame. ",
         "Available columns: ", paste(names(df), collapse = ", "))
  }

  if (is.null(conversion_map) || length(conversion_map) == 0) {
    stop("Argument 'conversion_map' must be a non-empty named list or vector.")
  }

  if (is.null(names(conversion_map)) || any(names(conversion_map) == "")) {
    stop("All elements in 'conversion_map' must be named.")
  }

  # Check that conversion_map values are atomic
  if (!all(sapply(conversion_map, function(x) is.atomic(x) && length(x) == 1))) {
    stop("All values in 'conversion_map' must be atomic (single character, numeric, or logical values).")
  }

  # Create a working copy to avoid modifying the original data frame
  df_copy <- df

  # Handle empty data frame or column with all NA values
  if (nrow(df_copy) == 0) {
    warning("Data frame is empty. Returning original data frame.")
    return(df_copy)
  }

  if (all(is.na(df_copy[[column]]))) {
    warning("Column '", column, "' contains only NA values. Converting to factor with NA level.")
    df_copy[[column]] <- factor(df_copy[[column]])
    return(df_copy)
  }

  # Convert the column to character to handle mixed data types
  # Store original column for potential error reporting
  original_column <- df_copy[[column]]

  tryCatch({
    df_copy[[column]] <- as.character(df_copy[[column]])
  }, error = function(e) {
    stop("Failed to convert column '", column, "' to character: ", e$message)
  })

  # Apply the conversions using the provided map
  # Use vectorized approach for better performance
  conversion_keys <- names(conversion_map)
  conversion_values <- as.character(unname(conversion_map))

  tryCatch({
    for (i in seq_along(conversion_keys)) {
      key <- conversion_keys[i]
      value <- conversion_values[i]

      # Handle potential NA values in the key matching
      matches <- !is.na(df_copy[[column]]) & df_copy[[column]] == key
      df_copy[[column]][matches] <- value
    }
  }, error = function(e) {
    stop("Error during recoding: ", e$message)
  })

  # Identify the recoded levels in the order they appear in conversion_map
  ordered_levels <- conversion_values

  # Get all other unique values in the column (excluding NA values for sorting)
  all_unique <- unique(df_copy[[column]])
  other_levels <- all_unique[!all_unique %in% ordered_levels & !is.na(all_unique)]

  # Sort other levels alphabetically
  if (length(other_levels) > 0) {
    other_levels <- sort(other_levels)
  }

  # Combine the ordered levels with the sorted other levels
  # Include NA if present in the data
  all_levels <- c(ordered_levels, other_levels)
  if (any(is.na(df_copy[[column]]))) {
    all_levels <- c(all_levels, NA)
  }

  # Convert the column to a factor with the defined levels
  tryCatch({
    df_copy[[column]] <- factor(df_copy[[column]], levels = all_levels)
  }, error = function(e) {
    stop("Failed to convert column '", column, "' to factor: ", e$message)
  })

  # Verify that the factorization was successful
  if (!is.factor(df_copy[[column]])) {
    stop("Failed to create factor for column '", column, "'.")
  }

  return(df_copy)
}
