# pondeR: Statistical Analysis Made Easier

`pondeR` is an R package providing a curated collection of functions for generating publication-ready descriptive and inferential statistics tables and figures — with minimal coding overhead. It is built upon the robust `gtsummary` framework and leverages the power of the `tidyverse`.

This package aims to make statistical analysis more accessible for everyone, ensuring transparency in analytical methods and providing a consistent template for results.

## Key Features

`pondeR` offers a comprehensive suite of tools for various stages of data analysis and general statistical reporting:

*   **Summarizing Data & Checking Distributions**: Generate publication-ready summary tables, calculate anthropometric indices according to WHO standards, and test for normality or skewness.
*   **Comparing Groups & Identifying Associations**: Provides methods for detecting differences between cohorts, calculating fold changes, and testing statistical associations. Includes plotting functions like volcano plots and heatmaps.
*   **Predictive Modeling & Classification**: Offers tools for regularized regression (Ridge, LASSO, Elastic Net), Linear Mixed-Effect Models, and Logistic regression analysis with performance evaluation (AUC/AUROC) and automated bias correction options.
*   **Multivariate Exploration and Discrimination**: Features dimensionality reduction and visualization tools (PCA, PLS, PCoA) for simplified exploration of complex datasets, including scree and score plots.
*   **Plotting Functions**: A variety of plotting utilities, including distribution plots, mean maps, and multi-Y plots.
*   **Export Functions**: Easily export data frames to Excel/CSV files and figures to PNGs.
*   **Utilities & Helpers**: A collection of utility functions for data management, relabelling, and randomization.

## Installation

You can install the development version of `pondeR` directly from GitHub:

```R
# Install pak if you haven't already
install.packages("pak")

# Install the GitHub package
pak::pak("jllcalorio/pondeR")
```

## Quick Start

Here's a quick example of how to generate a descriptive summary table using `pondeR` (assuming you have a data frame named `my_data`):

```R
library(pondeR)
library(dplyr) # Often useful with gtsummary-based packages

# Example data (replace with your actual data)
my_data <- tibble::tibble(
  age = rnorm(100, 30, 5),
  gender = sample(c("Male", "Female"), 100, replace = TRUE),
  treatment = sample(c("A", "B", "C"), 100, replace = TRUE),
  score = runif(100, 0, 100)
)

# Generate a summary table
my_data %>%
  run_summarytable(
    # Specify variables and grouping if needed
    include = c(age, gender, score),
    by = treatment
  )
```
*(Note: The exact arguments for `run_summarytable` may vary. Please refer to the package documentation for detailed usage.)*

## Contributing

We welcome contributions to `pondeR`! Please review our Code of Conduct before contributing.

## Reporting Bugs

If you encounter any bugs or have feature requests, please open an issue on our GitHub Issues page.

## License

This project is licensed under the MIT License.

## Acknowledgements

Huge thanks to The R Project for Statistical Computing, Positron, RStudio, the developers of the gtsummary, dplyr, ggplot2, everything in the tidyverse package, and all the dependencies I fail to mention. Your work makes `pondeR` possible.

## Future Development

More functions are continuously being added to `pondeR` to further enhance its capabilities and ease of use. Stay tuned for updates!
