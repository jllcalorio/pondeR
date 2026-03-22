# Transform Metabolomics Data

Applies various transformation methods to stabilize variance and reduce
heteroscedasticity in metabolomics data.

## Usage

``` r
run_transform(
  x,
  method = "log2",
  metadata = NULL,
  groups = "Group",
  qc_types = c("QC", "SQC", "EQC"),
  num_cores = 1,
  verbose = TRUE
)
```

## Arguments

- x:

  Matrix or data frame. Numeric data with samples in rows and features
  in columns.

- method:

  Character. Transformation method to apply. Options:

  - `"log2"`: Log base 2 transformation

  - `"log10"`: Log base 10 transformation

  - `"sqrt"`: Square root transformation

  - `"cbrt"`: Cube root transformation

  - `"vsn"`: Variance Stabilizing Normalization (requires 'vsn' package)

  - `"glog"`: Generalized logarithm transformation (requires 'pmp'
    package)

- metadata:

  Data frame. Sample metadata with number of rows equal to nrow(x).
  Required for glog transformation to identify QC samples. Default:
  NULL.

- groups:

  Character. Name of column in `metadata` containing group labels.
  Required for glog transformation. Default: "Group".

- qc_types:

  Character vector. Group labels identifying QC samples for glog.
  Default: c("QC", "SQC", "EQC").

- num_cores:

  Integer or "max". Number of CPU cores to use for VSN transformation.
  Use "max" for automatic detection (uses all available - 2). Default:
  1.

- verbose:

  Logical. Print progress messages. Default: TRUE.

## Value

A list of class "run_transform" containing:

- data:

  Matrix or data frame of transformed data (same class as input `x`)

- method_used:

  Character describing transformation method applied

- shift_applied:

  Numeric value added to data before transformation (if applicable)

- parameters:

  List of parameters used

## Details

**Transformation Methods:**

- **log2 / log10**: Logarithmic transformations compress large values
  and expand small values. Useful for data spanning multiple orders of
  magnitude. If data contains zeros or negative values, a shift is
  applied first.

- **sqrt**: Square root transformation. Moderate variance stabilization.
  Works with zero values but requires non-negative data.

- **cbrt**: Cube root transformation. Similar to sqrt but can handle
  negative values. Provides gentler compression than log.

- **vsn**: Variance Stabilizing Normalization. Model-based method from
  microarray analysis that calibrates variance across intensity range.
  Requires 'vsn' package. Can be computationally intensive; use
  `num_cores` for parallel processing.

- **glog**: Generalized logarithm transformation. Hybrid log
  transformation that handles low values better than standard log.
  Requires 'pmp' package and metadata with QC sample labels.

**When to Use Each Method:**

- **log2/log10**: Most common for metabolomics. Use when variance
  increases with mean (heteroscedasticity).

- **sqrt**: When variance is proportional to mean.

- **vsn**: When you need sophisticated variance modeling. Recommended
  for large datasets with clear intensity-dependent variance.

- **glog**: Alternative to log when you have many low-intensity values.
  Good for LC-MS data with detection limits.

## References

Huber, W., von Heydebreck, A., Sueltmann, H., Poustka, A., & Vingron, M.
(2002). Variance stabilization applied to microarray data calibration
and to the quantification of differential expression. Bioinformatics,
18(Suppl 1), S96-S104.
[doi:10.1093/bioinformatics/18.suppl_1.s96](https://doi.org/10.1093/bioinformatics/18.suppl_1.s96)

Parsons, H.M., Ludwig, C., Günther, U.L., & Viant, M.R. (2007). Improved
classification accuracy in 1- and 2-dimensional NMR metabolomics data
using the variance stabilising generalised logarithm transformation. BMC
Bioinformatics, 8, 234.
[doi:10.1186/1471-2105-8-234](https://doi.org/10.1186/1471-2105-8-234)

## Author

John Lennon L. Calorio

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data
set.seed(123)
x <- matrix(abs(rnorm(100 * 50, mean = 100, sd = 50)),
            nrow = 100, ncol = 50)
colnames(x) <- paste0("Feature", 1:50)

# Log transformation
result1 <- run_transform(x, method = "log2")

# VSN with parallel processing
result2 <- run_transform(x, method = "vsn", num_cores = 4)

# Glog transformation (requires metadata)
metadata <- data.frame(
  Sample = paste0("S", 1:100),
  Group = rep(c("Control", "Treatment", "QC"), c(40, 40, 20))
)
result3 <- run_transform(x, method = "glog", metadata = metadata)
} # }
```
