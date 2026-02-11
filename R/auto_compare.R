#' Automatic Statistical Comparison with Comprehensive Analysis
#'
#' This function performs an automatic comparison of a numeric outcome variable across
#' two or more groups. It intelligently selects the appropriate statistical test
#' (parametric or non-parametric) based on checks for normality (of data or residuals) and homogeneity of
#' variances, calculates effect sizes, and performs post-hoc tests when appropriate.
#'
#' @details
#' The function follows a decision-making process to determine the most appropriate test:
#'
#' \strong{Test Family Selection:}
#' If `paired = FALSE`, it compares independent groups. If there are 2 groups, it
#' selects a two-sample test family. If there are more than 2 groups, it selects
#' an ANOVA-family test.
#'
#' \strong{Assumption Checking (only for \code{test_type = "auto"}):}
#' \itemize{
#'   \item \strong{Normality:} Assessed using the Shapiro-Wilk test for small samples
#'     (n < 50). For larger samples (n >= 50), uses the Lilliefors (Kolmogorov-Smirnov)
#'     test from the \code{nortest} package. For two groups, each group is tested.
#'     For more than two groups, the normality of model residuals is tested.
#'   \item \strong{Homogeneity of Variance:} For independent tests, assessed using
#'     Levene's test from the \code{car} package. If \code{car} is not installed,
#'     falls back to Bartlett's test.
#' }
#'
#' \strong{Test Selection Logic (\code{test_type = "auto"}):}
#' \itemize{
#'   \item \strong{Independent Two-Sample:}
#'     \itemize{
#'       \item If normality is violated -> Mann-Whitney U test (Wilcoxon rank-sum).
#'       \item If normality is OK but variance is violated -> Welch's t-test.
#'       \item If both are OK -> Student's t-test.
#'     }
#'   \item \strong{Paired Two-Sample:}
#'     \itemize{
#'       \item If normality is violated -> Wilcoxon signed-rank test.
#'       \item If normality is OK -> Paired t-test.
#'     }
#'   \item \strong{Independent ANOVA (>2 groups):}
#'     \itemize{
#'       \item If normality is violated -> Kruskal-Wallis test.
#'       \item If normality is OK but variance violated -> Welch's ANOVA.
#'       \item If both OK -> One-way ANOVA.
#'     }
#' }
#'
#' \strong{Post-Hoc Tests:}
#' For significant results with >2 groups, appropriate post-hoc tests are automatically performed:
#' \itemize{
#'   \item Tukey HSD for one-way ANOVA
#'   \item Games-Howell for Welch's ANOVA
#'   \item Dunn's test for Kruskal-Wallis
#' }
#'
#' \strong{Effect Sizes:}
#' The function automatically calculates appropriate effect sizes:
#' \itemize{
#'   \item Cohen's d for t-tests
#'   \item Rank-biserial correlation for Mann-Whitney U
#'   \item Eta-squared for ANOVA
#'   \item Epsilon-squared for Kruskal-Wallis
#' }
#'
#' @param data A data frame containing the variables for analysis.
#' @param outcome A character string specifying the name of the numeric outcome variable.
#' @param group A character string specifying the name of the grouping factor variable.
#' @param paired A logical indicating whether the observations are paired. Default is FALSE.
#' @param test_type A character string specifying the test strategy. Must be one of
#'   \code{"auto"} (default), \code{"parametric"}, or \code{"nonparametric"}.
#' @param normality_method Method for assessing normality: "shapiro" (Shapiro-Wilk test),
#'   "lilliefors" (Lilliefors (Kolmogorov-Smirnov) test), or "auto" (uses Shapiro-Wilk for n < 50,
#'   Lilliefors for n >= 50). Default is "auto".
#' @param alpha_normality Significance level for normality tests. Default is 0.05.
#' @param alpha_variance Significance level for variance homogeneity tests. Default is 0.05.
#' @param test_alpha Significance level for the final statistical test. Default is 0.05.
#' @param min_n_threshold Minimum sample size required per group. Default is 3.
#' @param calculate_effect_size Logical indicating whether to calculate effect sizes. Default is TRUE.
#' @param perform_posthoc Logical indicating whether to perform post-hoc tests for >2 groups. Default is TRUE.
#' @param p_adjust_method The p-value adjustment method for post-hoc comparisons: "holm", "hochberg",
#'   "hommel", "bonferroni", "BH", "BY", "fdr", "none". Default is "BH".
#' @param group_order Character vector specifying the order of groups. If NULL, uses factor levels
#'   in alphabetical order.
#' @param verbose A logical indicating whether to print detailed messages. Default is TRUE.
#' @param num_cores Integer specifying the number of cores to use for parallel processing, or the character string "max".
#'   If "max", the function uses the number of available cores minus 2 (to reserve system resources), with a minimum of 1.
#'   Default is 1. Supports both Windows and POSIX (Mac/Linux) systems.
#'
#' @return An object of class "auto_compare" which is a list containing:
#'   \item{test_used}{Character string of the statistical test performed.}
#'   \item{test_result}{List containing results of the statistical test.}
#'   \item{effect_size}{List containing effect size estimate, confidence interval, magnitude, and interpretation.}
#'   \item{posthoc_result}{Data frame of post-hoc pairwise comparisons (NULL if not applicable).}
#'   \item{assumptions}{Data frame summarizing assumption check results.}
#'   \item{data_summary}{Data frame with comprehensive descriptive statistics per group.}
#'   \item{outcome}{Name of the outcome variable.}
#'   \item{group_var}{Name of the group variable.}
#'   \item{paired}{The paired setting used.}
#'   \item{test_type}{The test_type strategy used.}
#'   \item{test_alpha}{The significance level used.}
#'   \item{warnings}{Character vector of any warnings generated.}
#'   \item{parameters}{List of all parameters used in the analysis.}
#'   \item{raw_data}{Data frame containing the variables used.}
#'
#' @importFrom car leveneTest
#' @importFrom effectsize cohens_d rank_biserial eta_squared epsilon_squared interpret
#' @importFrom rstatix tukey_hsd games_howell_test dunn_test
#' @importFrom moments skewness kurtosis
#' @importFrom nortest lillie.test
#' @importFrom tibble tibble
#' @importFrom stats aov as.formula bartlett.test complete.cases kruskal.test lm
#' @importFrom stats median oneway.test quantile sd shapiro.test t.test wilcox.test
#' @importFrom utils install.packages
#'
#' @references Aaron Schlege, \emph{Games-Howell Post-Hoc Test}. \url{https://rpubs.com/aaronsc32/games-howell-test}.
#' @references Albers, C., & Lakens, D. (2018). \emph{When power analyses based on pilot data are biased: Inaccurate effect size estimators and follow-up bias}. Journal of Experimental Social Psychology. URL \url{https://doi.org/10.1016/j.jesp.2017.09.004}
#' @references Algina, J., Keselman, H. J., & Penfield, R. D. (2006). \emph{Confidence intervals for an effect size when variances are not equal}. Journal of Modern Applied Statistical Methods, 5(1), 2. URL \url{https://jmasm.com/index.php/jmasm/article/view/222} FULL PDF \url{https://digitalcommons.wayne.edu/cgi/viewcontent.cgi?article=1256&context=jmasm}
#' @references Allen, R. (2017). \emph{Statistics and Experimental Design for Psychologists: A Model Comparison Approach}. World Scientific Publishing Company.
#' @references Bache S, Wickham H (2022). \emph{_magrittr: A Forward-Pipe Operator for R_}. doi:10.32614/CRAN.package.magrittr <https://doi.org/10.32614/CRAN.package.magrittr>, R package version 2.0.3, <https://CRAN.R-project.org/package=magrittr>.
#' @references Bartlett, M. S. (1937). \emph{Properties of sufficiency and statistical tests}. Proceedings of the Royal Society of London Series A 160, 268–282. doi:10.1098/rspa.1937.0109.
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
#' @references Benjamini, Y., and Hochberg, Y. (1995). \emph{Controlling the false discovery rate: a practical and powerful approach to multiple testing}. Journal of the Royal Statistical Society Series B, 57, 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x.
#' @references Benjamini, Y., and Yekutieli, D. (2001). \emph{The control of the false discovery rate in multiple testing under dependency}. Annals of Statistics, 29, 1165–1188. doi:10.1214/aos/1013699998.
#' @references Ben-Shachar M, Lüdecke D, Makowski D (2020). \emph{effectsize: Estimation of Effect Size Indices and Standardized Parameters}. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815
#' @references Chambers, J. M., Freeny, A and Heiberger, R. M. (1992) \emph{Analysis of variance; designed experiments}. Chapter 5 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' @references Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}. Wadsworth & Brooks/Cole.
#' @references Cliff, N. (1993). \emph{Dominance statistics: Ordinal analyses to answer ordinal questions}. Psychological bulletin, 114(3), 494.
#' @references Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences (2nd Ed.)}. New York: Routledge.
#' @references Cureton, E. E. (1956). \emph{Rank-biserial correlation}. Psychometrika, 21(3), 287-290.
#' @references Dallal, G.E. and Wilkinson, L. (1986): \emph{An analytic approximation to the distribution of Lilliefors' test for normality}. The American Statistician, 40, 294–296.
#' @references David F. Bauer (1972). \emph{Constructing confidence sets using rank statistics}. Journal of the American Statistical Association 67, 687–690. doi:10.1080/01621459.1972.10481279.
#' @references Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021, May 7). \emph{Why Hedges’ g*s based on the non-pooled standard deviation should be reported with Welch's t-test}. doi:10.31234/osf.io/tu6mp
#' @references Dunn, O. J. (1964) \emph{Multiple comparisons using rank sums Technometrics}, 6(3):241-252.
#' @references Eric Langford (2006) \emph{Quartiles in Elementary Statistics}, Journal of Statistics Education 14 3; doi:10.1080/10691898.2006.11910589
#' @references Fox, J. (2016) \emph{Applied Regression Analysis and Generalized Linear Models}, Third Edition. Sage.
#' @references Fox J, Weisberg S (2019). \emph{_An R Companion to Applied Regression_}, Third edition. Sage, Thousand Oaks CA. <https://www.john-fox.ca/Companion/>.
#' @references Glass, G. V. (1965). \emph{A ranking variable analogue of biserial correlation: Implications for short-cut item analysis}. Journal of Educational Measurement, 2(1), 91-95.
#' @references Gross J, Ligges U (2015). \emph{_nortest: Tests for Normality_}. doi:10.32614/CRAN.package.nortest <https://doi.org/10.32614/CRAN.package.nortest>, R package version 1.0-4, <https://CRAN.R-project.org/package=nortest>.
#' @references Hedges, L. V. & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. Orlando, FL: Academic Press.
#' @references Hochberg, Y. (1988). \emph{A sharper Bonferroni procedure for multiple tests of significance}. Biometrika, 75, 800–803. doi:10.2307/2336325
#' @references Holm, S. (1979). \emph{A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics}, 6, 65–70. \url{https://www.jstor.org/stable/4615733}.
#' @references Hommel, G. (1988). \emph{A stagewise rejective multiple test procedure based on a modified Bonferroni test}. Biometrika, 75, 383–386. doi:10.2307/2336190.
#' @references Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis: Correcting error and bias in research findings}. Sage.
#' @references Hyndman, R. J. and Fan, Y. (1996) \emph{Sample quantiles in statistical packages}, American Statistician 50, 361–365. doi:10.2307/2684934.
#' @references Kassambara A (2023). \emph{_rstatix: Pipe-Friendly Framework for Basic Statistical Tests_}. doi:10.32614/CRAN.package.rstatix <https://doi.org/10.32614/CRAN.package.rstatix>, R package version 0.7.2, <https://CRAN.R-project.org/package=rstatix>.
#' @references Kelley, T. (1935) \emph{An unbiased correlation ratio measure}. Proceedings of the National Academy of Sciences. 21(9). 554-559.
#' @references Kerby, D. S. (2014). \emph{The simple difference formula: An approach to teaching nonparametric correlation}. Comprehensive Psychology, 3, 11-IT.
#' @references King, B. M., & Minium, E. W. (2008). \emph{Statistical reasoning in the behavioral sciences}. John Wiley & Sons Inc.
#' @references Komsta L, Novomestky F (2022). \emph{_moments: Moments, Cumulants, Skewness, Kurtosis and Related Tests_}. doi:10.32614/CRAN.package.moments <https://doi.org/10.32614/CRAN.package.moments>, R package version 0.14.1, <https://CRAN.R-project.org/package=moments>.
#' @references Makkonen, L. and Pajari, M. (2014) \emph{Defining Sample Quantiles by the True Rank Probability}, Journal of Probability and Statistics; Hindawi Publ.Corp. doi:10.1155/2014/326579
#' @references Myles Hollander and Douglas A. Wolfe (11973). \emph{Nonparametric Statistical Methods}. New York: John Wiley & Sons. Pages 27–33 (one-sample), 68–75 (two-sample). Or second edition (1999).
#' @references Myles Hollander and Douglas A. Wolfe (1973), \emph{Nonparametric Statistical Methods}. New York: John Wiley & Sons. Pages 115–120.
#' @references Müller K, Wickham H (2025). \emph{_tibble: Simple Data Frames_}. doi:10.32614/CRAN.package.tibble <https://doi.org/10.32614/CRAN.package.tibble>, R package version 3.3.0, <https://CRAN.R-project.org/package=tibble>.
#' @references Olejnik, S., & Algina, J. (2003). \emph{Generalized eta and omega squared statistics: measures of effect size for some common research designs}. Psychological methods, 8(4), 434.
#' @references Patrick Royston (1982). \emph{Algorithm AS 181: The W test for Normality}. Applied Statistics, 31, 176–180. doi:10.2307/2347986.
#' @references Patrick Royston (1982). \emph{An extension of Shapiro and Wilk's W test for normality to large samples}. Applied Statistics, 31, 115–124. doi:10.2307/2347973.
#' @references Patrick Royston (1995). \emph{Remark AS R94: A remark on Algorithm AS 181: The W test for normality}. Applied Statistics, 44, 547–551. doi:10.2307/2986146.
#' @references Ruxton, G.D., and Beauchamp, G. (2008) \emph{‘Time for some a priori thinking about post hoc testing’}, Behavioral Ecology, 19(3), pp. 690-693. doi: 10.1093/beheco/arn020. In-text citations: (Ruxton and Beauchamp, 2008)
#' @references Sangseok Lee, Dong Kyu Lee. \emph{What is the proper way to apply the multiple comparison test?}. Korean J Anesthesiol. 2018;71(5):353-360.
#' @references Sarkar, S. (1998). \emph{Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture}. Annals of Statistics, 26, 494–504. doi:10.1214/aos/1028144846.
#' @references Sarkar, S., and Chang, C. K. (1997). \emph{The Simes method for multiple hypothesis testing with positively dependent test statistics}. Journal of the American Statistical Association, 92, 1601–1608. doi:10.2307/2965431.
#' @references Shaffer, J. P. (1995). \emph{Multiple hypothesis testing. Annual Review of Psychology}, 46, 561–584. doi:10.1146/annurev.ps.46.020195.003021. (An excellent review of the area.)
#' @references Steiger, J. H. (2004). \emph{Beyond the F test: Effect size confidence intervals and tests of close fit in the analysis of variance and contrast analysis}. Psychological Methods, 9, 164-182.
#' @references Stephens, M.A. (1974): \emph{EDF statistics for goodness of fit and some comparisons}. Journal of the American Statistical Association, 69, 730–737.
#' @references Thode Jr., H.C. (2002): \emph{Testing for Normality}. Marcel Dekker, New York.
#' @references Tomczak, M., & Tomczak, E. (2014). \emph{The need to report effect size estimates revisited}. An overview of some recommended measures of effect size.
#' @references Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). \emph{“Welcome to the tidyverse.”} _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
#' @references Wicklin, R. (2017) \emph{Sample quantiles: A comparison of 9 definitions}; SAS Blog. https://blogs.sas.com/content/iml/2017/05/24/definitions-sample-quantiles.htm
#' @references Wright, S. P. (1992). \emph{Adjusted P-values for simultaneous inference}. Biometrics, 48, 1005–1013. doi:10.2307/2532694. (Explains the adjusted P-value approach.)
#' @references R Core Team (2025). \emph{_R: A Language and Environment for Statistical Computing_}. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' # --- Basic Usage with 3+ Groups ---
#' # Using the classic 'iris' dataset. Sepal.Length is numeric, Species is a factor.
#' # This will likely trigger the One-way ANOVA -> Tukey HSD pathway.
#' data(iris)
#' res_iris <- auto_compare(data = iris, outcome = "Sepal.Length", group = "Species")
#' print(res_iris)
#' summary(res_iris)
#'
#' # --- Two-Sample Independent Comparison ---
#' # Using 'mtcars' to compare miles per gallon (mpg) between automatic (0) and manual (1) cars.
#' # The variances are unequal, so this will trigger the Welch's t-test.
#' data(mtcars)
#' mtcars$am <- factor(mtcars$am, labels = c("automatic", "manual"))
#' res_mtcars <- auto_compare(data = mtcars, outcome = "mpg", group = "am")
#' print(res_mtcars)
#'
#' # --- Paired-Sample Comparison ---
#' # Using the 'sleep' dataset, which is designed for a paired t-test.
#' # It compares the extra hours of sleep for the same subjects under two drugs.
#' data(sleep)
#' # We specify the subject ID is implicit in the row order for a paired test.
#' res_sleep <- auto_compare(data = sleep, outcome = "extra", group = "group", paired = TRUE)
#' summary(res_sleep)
#'
#' # --- Forcing a Non-Parametric Test ---
#' # Use 'iris' data again, but force the non-parametric Kruskal-Wallis test.
#' data(iris)
#' res_iris_np <- auto_compare(
#'   data = iris,
#'   outcome = "Sepal.Length",
#'   group = "Species",
#'   test_type = "nonparametric"
#' )
#' print(res_iris_np)
#'
#' # --- Forcing a Parametric Test & Custom P-Value Adjustment ---
#' # Use the 'ToothGrowth' data with 'dose' as the grouping factor.
#' # We force an ANOVA (even if normality is slightly off) and use Bonferroni correction for post-hocs.
#' data(ToothGrowth)
#' # Convert dose to a factor to be treated as a group
#' ToothGrowth$dose <- as.factor(ToothGrowth$dose)
#' res_tooth <- auto_compare(
#'   data = ToothGrowth,
#'   outcome = "len",
#'   group = "dose",
#'   test_type = "parametric",
#'   p_adjust_method = "bonferroni"
#' )
#' # Look at the post-hoc results and the adjusted p-values
#' print(res_tooth$posthoc_result)
#'
#' # --- Demonstrating Custom Group Ordering ---
#' # Use 'PlantGrowth' and specify a non-alphabetical order for the groups.
#' data(PlantGrowth)
#' custom_order <- c("trt1", "ctrl", "trt2")
#' res_plant <- auto_compare(
#'   data = PlantGrowth,
#'   outcome = "weight",
#'   group = "group",
#'   group_order = custom_order,
#'   verbose = FALSE # Suppress messages for cleaner output
#' )
#' # The summary and plots will respect the custom order
#' print(res_plant$data_summary)
#'
#' # --- Disabling Features for Quicker Analysis ---
#' # If you only need the main test result, you can disable effect sizes and post-hoc tests.
#' # Using 'mtcars' again with the 'cyl' variable as a group.
#' data(mtcars)
#' mtcars$cyl <- as.factor(mtcars$cyl)
#' res_mtcars_fast <- auto_compare(
#'   data = mtcars,
#'   outcome = "mpg",
#'   group = "cyl",
#'   calculate_effect_size = FALSE,
#'   perform_posthoc = FALSE
#' )
#' print(res_mtcars_fast)
#'
#' # --- Generating Plots ---
#' # The function includes a plot method to visualize the results.
#' \dontrun{
#' # Run an analysis first
#' data(PlantGrowth)
#' res_plot <- auto_compare(data = PlantGrowth, outcome = "weight", group = "group")
#'
#' # Generate the plots
#' plot(res_plot)
#' # The plot will show:
#' # 1. A boxplot with means
#' # 2. A Q-Q plot for checking normality
#' # 3. A heatmap of post-hoc significant differences
#' }
#'
#' @export
auto_compare <- function(
    data,
    outcome,
    group,
    paired                = FALSE,
    test_type             = c("auto", "parametric", "nonparametric"),
    normality_method      = c("auto", "shapiro", "lilliefors"),
    alpha_normality       = 0.05,
    alpha_variance        = 0.05,
    test_alpha            = 0.05,
    min_n_threshold       = 3,
    calculate_effect_size = TRUE,
    perform_posthoc       = TRUE,
    p_adjust_method       = "BH",
    group_order           = NULL,
    verbose               = TRUE,
    num_cores             = 1
) {

  # --- 1. Input Validation and Preparation ---
  test_type        <- match.arg(test_type)
  normality_method <- match.arg(normality_method)
  warnings_list    <- character(0)

  # Validate num_cores and check OS
  avail_cores                         <- parallel::detectCores()
  if (is.na(avail_cores)) avail_cores <- 1

  # Handle "max" or Integer input
  if (identical(num_cores, "max")) {
    # Use max cores - 2, but ensure at least 1 core is used
    num_cores <- max(1, avail_cores - 2)
    if (verbose) message(sprintf("Parallel processing: Using 'max' setting (%d cores active, 2 reserved).", num_cores))
  } else if (is.numeric(num_cores)) {
    # Validate integer range
    if (num_cores < 1 || num_cores > avail_cores || num_cores %% 1 != 0) {
      stop(sprintf("'num_cores' must be an integer between 1 and %d, or the string 'max'.", avail_cores))
    }
    num_cores <- as.integer(num_cores)
    if (num_cores > 1 && verbose) message(sprintf("Parallel processing: Using %d cores.", num_cores))
  } else {
    stop("Invalid 'num_cores' argument. Must be an integer or 'max'.")
  }

  # Validata data
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!all(c(outcome, group) %in% names(data))) {
    stop(paste0("Outcome variable '", outcome, "' or group variable '", group, "' not found in data."))
  }

  y   <- data[[outcome]]
  grp <- as.factor(data[[group]])

  if (!is.numeric(y)) {
    stop("Outcome variable must be numeric.")
  }

  # Apply group order if specified
  if (!is.null(group_order)) {
    if (!all(group_order %in% levels(grp))) {
      stop("group_order contains levels not present in the data")
    }
    grp <- factor(grp, levels = group_order)
  }

  # Remove missing values
  complete_cases <- complete.cases(y, grp)
  if (sum(!complete_cases) > 0) {
    msg           <- sprintf("Removing %d rows with missing values.", sum(!complete_cases))
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }
  y   <- y[complete_cases]
  grp <- droplevels(grp[complete_cases])

  if (nlevels(grp) < 2) {
    stop("Grouping factor must have at least 2 levels after removing missing values.")
  }

  groups   <- levels(grp)
  n_groups <- length(groups)
  group_ns <- table(grp)
  min_n    <- min(group_ns)

  if (verbose) {
    message("\n=== Auto Compare Analysis ===")
    message(sprintf("Outcome: %s | Group: %s | Paired: %s", outcome, group, paired))
    message(sprintf("Number of groups: %d", n_groups))
    message(sprintf("Sample sizes: %s", paste(names(group_ns), "=", group_ns, collapse = ", ")))
    message(sprintf("Test Type: %s", test_type))
  }

  # --- 2. Initialize Assumption Checks and Determine Test Family ---
  assumptions <- data.frame(
    check            = character(),
    result           = character(),
    p_value          = numeric(),
    decision         = character(),
    stringsAsFactors = FALSE
  )

  # Check minimum sample size
  assumptions <- rbind(assumptions, data.frame(
    check            = "Minimum sample size",
    result           = sprintf("Min n = %d", min_n),
    p_value          = NA,
    decision         = ifelse(min_n >= min_n_threshold, "PASS", "WARNING"),
    stringsAsFactors = FALSE
  ))

  if (min_n < min_n_threshold) {
    msg           <- sprintf("Warning: Minimum sample size (%d) is below threshold (%d)", min_n, min_n_threshold)
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }

  # Decide test family
  if (n_groups == 2) {
    test_family <- ifelse(paired, "paired_two_sample", "independent_two_sample")
  } else if (n_groups > 2) {
    test_family <- ifelse(paired, "repeated_measures_anova", "independent_anova")
  }

  # For paired tests, ensure equal sample sizes
  if (paired && n_groups == 2) {
    if (group_ns[1] != group_ns[2]) {
      stop("For paired tests, groups must have equal sample sizes after removing missing values.")
    }
  }

  # --- 3. Helper Function: Check Normality ---
  check_normality <- function(data_vector, group_label = NULL, method = normality_method) {
    n <- length(data_vector)

    if (n < 3) {
      return(list(violated = TRUE, p_value = NA, method = "Insufficient data"))
    }

    # Determine which method to use
    use_method <- if (method == "auto") {
      if (n < 50) "shapiro" else "lilliefors"
    } else {
      method
    }

    if (use_method == "shapiro") {
      shapiro_result <- shapiro.test(data_vector)
      p_val          <- shapiro_result$p.value
      violated       <- p_val < alpha_normality
      return(list(
        violated = violated,
        p_value  = p_val,
        method   = "Shapiro-Wilk test",
        label    = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"
      ))
    } else if (use_method == "lilliefors") {
      # Use Lilliefors (Kolmogorov-Smirnov) test
      if (!requireNamespace("nortest", quietly = TRUE)) {
        warning("Package 'nortest' not available. Falling back to Shapiro-Wilk test.")
        shapiro_result <- shapiro.test(data_vector)
        p_val          <- shapiro_result$p.value
        violated       <- p_val < alpha_normality
        return(list(
          violated = violated,
          p_value  = p_val,
          method   = "Shapiro-Wilk test (fallback)",
          label    = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"
        ))
      }

      lillie_result <- nortest::lillie.test(data_vector)
      p_val         <- lillie_result$p.value
      violated      <- p_val < alpha_normality

      return(list(
        violated = violated,
        p_value  = p_val,
        method   = "Lilliefors (K-S) test",
        label    = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"
      ))
    }
  }

  # --- 4. Assumption Checking (only if test_type is "auto") ---
  normality_violated <- FALSE
  variance_violated  <- FALSE

  if (test_type == "auto") {
    # --- Normality Checks ---
    if (n_groups == 2) {
      # Parallel wrapper for group-wise normality
      run_norm_check <- function(g) {
        group_data <- y[grp == g]
        # Call the helper function defined in Section 3
        check_normality(group_data, g, normality_method)
      }

      # Run checks (Parallel or Sequential)
      norm_results <- list()

      if (num_cores > 1) {
        if (.Platform$OS.type == "windows") {
          # Windows: Use PSOCK cluster
          cl <- parallel::makeCluster(num_cores)
          on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

          # Export necessary variables/functions to workers
          parallel::clusterExport(cl,
                                  varlist = c("y", "grp", "check_normality", "normality_method", "alpha_normality"),
                                  envir = environment())

          # Load packages on workers
          parallel::clusterEvalQ(cl, {
            if (requireNamespace("nortest", quietly = TRUE)) library(nortest)
          })

          norm_results <- parallel::parLapply(cl, groups, run_norm_check)
        } else {
          # Mac/Linux: Use Forking
          norm_results <- parallel::mclapply(groups, run_norm_check, mc.cores = num_cores)
        }
      } else {
        # Sequential
        norm_results <- lapply(groups, run_norm_check)
      }

      # Process results
      for (norm_check in norm_results) {
        # Error handling for failed workers
        if (inherits(norm_check, "try-error")) {
          if (verbose) warning("A parallel worker failed. Defaulting assumption to OK.")
          next
        }

        assumptions <- rbind(assumptions, data.frame(
          check = norm_check$label, result = norm_check$method,
          p_value = norm_check$p_value, decision = ifelse(norm_check$violated, "VIOLATED", "OK"),
          stringsAsFactors = FALSE
        ))
        if (norm_check$violated) normality_violated <- TRUE
      }
    } else {
      # ANOVA: check normality of residuals (Single vector, no parallelization needed)
      fit <- lm(y ~ grp)
      norm_check <- check_normality(fit$residuals, NULL, normality_method)
      norm_check$label <- "Normality of residuals"
      assumptions <- rbind(assumptions, data.frame(
        check = norm_check$label, result = norm_check$method,
        p_value = norm_check$p_value, decision = ifelse(norm_check$violated, "VIOLATED", "OK"),
        stringsAsFactors = FALSE
      ))
      if (norm_check$violated) normality_violated <- TRUE
    }

    # Homogeneity of variance check (for independent tests only)
    if (!paired) {
      variance_check_result <- tryCatch(
        {
          if (requireNamespace("car", quietly = TRUE)) {
            car::leveneTest(y ~ grp)$`Pr(>F)`[1]
          } else {
            if (verbose) message("Package 'car' not found. Using Bartlett's test for variance homogeneity.")
            bartlett.test(y ~ grp)$p.value
          }
        },
        error = function(e) {
          msg            <- "Error in variance test. Assuming homogeneity."
          warnings_list <<- c(warnings_list, msg)
          if (verbose) message(msg)
          return(1)
        }
      )

      p_val          <- variance_check_result
      decision       <- ifelse(p_val < alpha_variance, "VIOLATED", "OK")
      test_used_name <- if (requireNamespace("car", quietly = TRUE)) "Levene's test" else "Bartlett's test"

      assumptions <- rbind(assumptions, data.frame(
        check            = "Homogeneity of variance",
        result           = test_used_name,
        p_value          = p_val,
        decision         = decision,
        stringsAsFactors = FALSE
      ))

      if (decision == "VIOLATED") variance_violated <- TRUE
    }
  }

  # --- 5. Select and Perform Test ---
  test_used       <- NULL
  test_result     <- NULL
  parametric_used <- TRUE
  fitted_model    <- NULL # Initialize to store fitted model for effect size/post-hoc

  if (test_type == "auto") {
    if (test_family == "independent_two_sample") {
      if (normality_violated) {
        test_used   <- "Mann-Whitney U test (Wilcoxon rank-sum)"
        test_result <- wilcox.test(y ~ grp)
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's t-test"
        test_result <- t.test(y ~ grp, var.equal = FALSE)
      } else {
        test_used   <- "Student's t-test"
        test_result <- t.test(y ~ grp, var.equal = TRUE)
      }
    } else if (test_family == "paired_two_sample") {
      group1_data <- y[grp == groups[1]]
      group2_data <- y[grp == groups[2]]
      if (normality_violated) {
        test_used       <- "Wilcoxon signed-rank test"
        test_result     <- wilcox.test(group1_data, group2_data, paired = TRUE)
        parametric_used <- FALSE
      } else {
        test_used   <- "Paired t-test"
        test_result <- t.test(group1_data, group2_data, paired = TRUE)
      }
    } else if (test_family == "independent_anova") {
      if (normality_violated) {
        test_used       <- "Kruskal-Wallis test"
        test_result     <- kruskal.test(y ~ grp)
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's ANOVA"
        test_result <- oneway.test(y ~ grp, var.equal = FALSE)
        # Note: fitted_model remains NULL for Welch's ANOVA
      } else {
        test_used    <- "One-way ANOVA"
        fitted_model <- aov(y ~ grp) # Store the model
        summary_aov  <- summary(fitted_model)[[1]]
        test_result  <- list(
          statistic = c("F value" = summary_aov$`F value`[1]),
          parameter = c("df1" = summary_aov$Df[1], "df2" = summary_aov$Df[2]),
          p.value   = summary_aov$`Pr(>F)`[1],
          method    = test_used,
          data.name = deparse(substitute(y))
        )
      }
    }
  } else if (test_type == "parametric") {
    if (test_family == "independent_two_sample") {
      test_used   <- "Student's t-test (forced)"
      test_result <- t.test(y ~ grp, var.equal = TRUE)
    } else if (test_family == "paired_two_sample") {
      group1_data <- y[grp == groups[1]]
      group2_data <- y[grp == groups[2]]
      test_used   <- "Paired t-test (forced)"
      test_result <- t.test(group1_data, group2_data, paired = TRUE)
    } else if (test_family == "independent_anova") {
      test_used    <- "One-way ANOVA (forced)"
      fitted_model <- aov(y ~ grp) # Store the model
      summary_aov  <- summary(fitted_model)[[1]]
      test_result  <- list(
        statistic = c("F value" = summary_aov$`F value`[1]),
        parameter = c("df1" = summary_aov$Df[1], "df2" = summary_aov$Df[2]),
        p.value   = summary_aov$`Pr(>F)`[1],
        method    = test_used,
        data.name = deparse(substitute(y))
      )
    }
  } else if (test_type == "nonparametric") {
    parametric_used <- FALSE
    if (test_family == "independent_two_sample") {
      test_used   <- "Mann-Whitney U test (forced)"
      test_result <- wilcox.test(y ~ grp)
    } else if (test_family == "paired_two_sample") {
      group1_data <- y[grp == groups[1]]
      group2_data <- y[grp == groups[2]]
      test_used   <- "Wilcoxon signed-rank test (forced)"
      test_result <- wilcox.test(group1_data, group2_data, paired = TRUE)
    } else if (test_family == "independent_anova") {
      test_used   <- "Kruskal-Wallis test (forced)"
      test_result <- kruskal.test(y ~ grp)
    }
  }

  if (test_family == "repeated_measures_anova") {
    stop("Repeated measures ANOVA is not yet implemented. This requires a subject ID and different data structuring.")
  }

  # --- 6. Calculate Effect Size ---
  effect_size_result <- NULL

  if (calculate_effect_size && !is.null(test_result)) {
    if (requireNamespace("effectsize", quietly = TRUE)) {
      effect_size_result <- tryCatch({
        if (test_family == "independent_two_sample") {
          if (parametric_used) {
            es <- effectsize::cohens_d(y ~ grp, pooled_sd = !variance_violated, verbose = FALSE)
            interp <- effectsize::interpret(es, rules = "cohen1988")
            list(
              estimate       = es$Cohens_d,
              ci_low         = es$CI_low,
              ci_high        = es$CI_high,
              # magnitude      = effectsize::interpret(es, rules = "cohen1988")[[2]],
              magnitude      = as.character(interp$Interpretation),
              metric         = "Cohen's d",
              interpretation = sprintf("Effect size: %s (%.3f)",
                                       # effectsize::interpret(es, rules = "cohen1988")[[2]],
                                       as.character(interp$Interpretation),
                                       es$Cohens_d)
            )
          } else {
            es <- effectsize::rank_biserial(y ~ grp, verbose = FALSE)
            interp <- effectsize::interpret(es, rules = "funder2019")
            list(
              estimate       = es$r_rank_biserial,
              ci_low         = es$CI_low,
              ci_high        = es$CI_high,
              # magnitude      = effectsize::interpret(es, rules = "funder2019")[[2]],
              magnitude      = as.character(interp$Interpretation),
              metric         = "Rank-biserial correlation",
              interpretation = sprintf("Effect size: %s (%.3f)",
                                       # effectsize::interpret(es, rules = "funder2019")[[2]],
                                       as.character(interp$Interpretation),
                                       es$r_rank_biserial)
            )
          }
        } else if (test_family == "paired_two_sample") {
          if (parametric_used) {
            es <- effectsize::cohens_d(y[grp == groups[1]], y[grp == groups[2]], paired = TRUE, verbose = FALSE)
            interp <- effectsize::interpret(es, rules = "cohen1988")
            list(
              estimate       = es$Cohens_d,
              ci_low         = es$CI_low,
              ci_high        = es$CI_high,
              # magnitude      = effectsize::interpret(es, rules = "cohen1988")[[2]],
              magnitude      = as.character(interp$Interpretation),
              metric         = "Cohen's d (paired)",
              interpretation = sprintf("Effect size: %s (%.3f)",
                                       # effectsize::interpret(es, rules = "cohen1988")[[2]],
                                       as.character(interp$Interpretation),
                                       es$Cohens_d)
            )
          } else {
            es <- effectsize::rank_biserial(y[grp == groups[1]], y[grp == groups[2]], paired = TRUE, verbose = FALSE)
            interp <- effectsize::interpret(es, rules = "funder2019")
            list(
              estimate       = es$r_rank_biserial,
              ci_low         = es$CI_low,
              ci_high        = es$CI_high,
              # magnitude      = effectsize::interpret(es, rules = "funder2019")[[2]],
              magnitude      = as.character(interp$Interpretation),
              metric         = "Rank-biserial correlation (paired)",
              interpretation = sprintf("Effect size: %s (%.3f)",
                                       # effectsize::interpret(es, rules = "funder2019")[[2]],
                                       as.character(interp$Interpretation),
                                       es$r_rank_biserial)
            )
          }
        } else if (test_family == "independent_anova") {
          if (parametric_used) {
            # Handle Welch's ANOVA where fitted_model is NULL
            if (is.null(fitted_model)) {
              es_model       <- aov(y ~ grp)
              warnings_list <<- c(warnings_list, "Effect size (Eta-squared) calculated from a standard ANOVA model, as the test used was Welch's ANOVA which does not provide a model object for effect size calculation.")
            } else {
              es_model <- fitted_model
            }
            es <- effectsize::eta_squared(es_model, verbose = FALSE)
            interp <- effectsize::interpret(es, rules = "field2013")
            list(
              estimate       = es$Eta2,
              ci_low         = es$CI_low,
              ci_high        = es$CI_high,
              # magnitude      = effectsize::interpret(es, rules = "field2013")[[2]],
              magnitude      = as.character(interp$Interpretation),
              metric         = "Eta-squared",
              interpretation = sprintf("Effect size: %s (%.3f)",
                                       # effectsize::interpret(es, rules = "field2013")[[2]],
                                       as.character(interp$Interpretation),
                                       es$Eta2)
            )
          } else {
            # Non-parametric ANOVA (Kruskal-Wallis)
            es <- effectsize::rank_epsilon_squared(
              x       = y,
              groups  = grp,
              verbose = FALSE
            )
            interp <- effectsize::interpret(es, rules = "field2013")
            list(
              estimate       = es$rank_epsilon_squared,
              ci_low         = if(!is.null(es$CI_low)) es$CI_low else NA,
              ci_high        = if(!is.null(es$CI_high)) es$CI_high else NA,
              # magnitude      = effectsize::interpret(es, rules = "field2013")[[2]],
              magnitude      = as.character(interp$Interpretation),
              metric         = "Rank epsilon-squared",
              interpretation = sprintf("Effect size: %s (%.3f)",
                                       # effectsize::interpret(es, rules = "field2013")[[2]],
                                       as.character(interp$Interpretation),
                                       es$rank_epsilon_squared)
            )
          }
        }
      }, error = function(e) {
        msg            <- paste("Could not calculate effect size:", e$message)
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message(msg)
        NULL
      })
    } else {
      msg           <- "Package 'effectsize' not available. Skipping effect size calculation."
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    }
  }

  # --- 7. Perform Post-Hoc Tests ---
  posthoc_result <- NULL

  if (perform_posthoc && n_groups > 2 && !is.null(test_result) && test_result$p.value < test_alpha) {
    if (requireNamespace("rstatix", quietly = TRUE)) {
      posthoc_result <- tryCatch({
        # Create tibble with consistent column names
        test_data <- tibble::tibble(y = y, grp = grp)
        ph <- NULL # Initialize ph

        if (test_used == "One-way ANOVA" || test_used == "One-way ANOVA (forced)") {
          # Tukey HSD
          ph              <- rstatix::tukey_hsd(test_data, y ~ grp)
          ph$posthoc_test <- "Tukey HSD"
          # Add placeholder columns for consistency
          ph$statistic    <- NA_real_
          ph$p            <- NA_real_

          # === START FIX ===
          # Manually add n1 and n2 for Tukey using a robust lookup

          # Create a named vector for fast lookup
          n_lookup        <- as.numeric(group_ns)
          names(n_lookup) <- names(group_ns)

          # Add n1 and n2 by looking up the group names
          # This is more robust than merge()
          ph$n1 <- n_lookup[ph$group1]
          ph$n2 <- n_lookup[ph$group2]
          # === END FIX ===

        } else if (test_used == "Welch's ANOVA") {
          # Games-Howell
          ph              <- rstatix::games_howell_test(test_data, y ~ grp)
          ph$posthoc_test <- "Games-Howell"
          # Add placeholder columns for consistency
          ph$p <- NA_real_

        } else if (test_used == "Kruskal-Wallis test" || test_used == "Kruskal-Wallis test (forced)") {
          # Dunn's test
          ph              <- rstatix::dunn_test(test_data, y ~ grp, p.adjust.method = p_adjust_method)
          ph$posthoc_test <- "Dunn's test"
        }

        # --- Standardize Columns ---

        # Add formatted p-values - check which columns exist
        if ("p" %in% names(ph)) {
          ph$p.format <- ifelse(ph$p < 0.001, "p < 0.001",
                                sprintf("p = %.3f", ph$p))
        } else {
          ph$p.format <- NA_character_
        }

        if ("p.adj" %in% names(ph)) {
          ph$p.adj.signif <- ifelse(ph$p.adj < 0.001, "***",
                                    ifelse(ph$p.adj < 0.01, "**",
                                           ifelse(ph$p.adj < 0.05, "*", "ns")))
          ph$p.adj.format <- ifelse(ph$p.adj < 0.001, "p < 0.001",
                                    sprintf("p = %.3f", ph$p.adj))
          ph$significant  <- ph$p.adj < test_alpha
        } else {
          # Handle cases like Tukey where only p.adj exists
          # Note: rstatix::tukey_hsd() returns 'p.adj'
          if ("p.adj" %in% names(ph)) {
            ph$p.adj.signif <- ifelse(ph$p.adj < 0.001, "***",
                                      ifelse(ph$p.adj < 0.01, "**",
                                             ifelse(ph$p.adj < 0.05, "*", "ns")))
            ph$p.adj.format <- ifelse(ph$p.adj < 0.001, "p < 0.001",
                                      sprintf("p = %.3f", ph$p.adj))
            ph$significant  <- ph$p.adj < test_alpha
          } else {
            ph$p.adj.signif <- NA_character_
            ph$p.adj.format <- NA_character_
            ph$significant  <- NA
          }
        }

        # Add summary statistics for each comparison
        summary_stats <- data.frame(
          group            = groups,
          mean             = tapply(y, grp, mean),
          median           = tapply(y, grp, median),
          stringsAsFactors = FALSE
        )

        # Merge with group1
        ph <- merge(ph, summary_stats, by.x = "group1", by.y = "group", all.x = TRUE)
        names(ph)[names(ph) %in% c("mean", "median")] <- c("mean1", "median1")

        # Merge with group2
        ph <- merge(ph, summary_stats, by.x = "group2", by.y = "group", all.x = TRUE)
        names(ph)[names(ph) %in% c("mean", "median")] <- c("mean2", "median2")

        # Add interpretation
        if ("p.adj.signif" %in% names(ph)) {
          if (parametric_used) {
            ph$interpretation <- ifelse(ph$p.adj.signif == "ns",
                                        "No significant difference between groups",
                                        ifelse(ph$mean1 < ph$mean2,
                                               paste0(ph$group1, " has lower mean than ", ph$group2),
                                               paste0(ph$group1, " has higher mean than ", ph$group2)))
          } else {
            ph$interpretation <- ifelse(ph$p.adj.signif == "ns",
                                        "No significant difference between groups",
                                        ifelse(ph$median1 < ph$median2,
                                               paste0(ph$group1, " has lower median than ", ph$group2),
                                               paste0(ph$group1, " has higher median than ", ph$group2)))
          }
        } else {
          ph$interpretation <- NA_character_
        }

        # --- Final Column Selection and Reordering ---
        final_cols <- c(
          "group1", "group2", "n1", "n2", "statistic", "p", "p.adj",
          "posthoc_test", "p.format", "p.adj.format", "p.adj.signif", "significant",
          "mean1", "mean2", "median1", "median2", "interpretation"
        )

        # Ensure all columns exist, adding NA for any that might be missing
        for (col in final_cols) {
          if (!col %in% names(ph)) {
            ph[[col]] <- NA
          }
        }

        # Select and reorder
        ph <- ph[, final_cols]

        ph

      }, error = function(e) {
        msg            <- paste("Could not perform post-hoc test:", e$message)
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message(msg)
        NULL
      })
    } else {
      msg           <- "Package 'rstatix' not available. Skipping post-hoc tests."
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    }
  }

  # --- 8. Comprehensive Data Summary ---
  data_summary <- data.frame(
    group     = groups,
    n         = as.numeric(group_ns),
    mean      = tapply(y, grp, mean),
    sd        = tapply(y, grp, sd),
    se        = tapply(y, grp, function(x) sd(x) / sqrt(length(x))),
    median    = tapply(y, grp, median),
    q25       = tapply(y, grp, function(x) quantile(x, 0.25)),
    q75       = tapply(y, grp, function(x) quantile(x, 0.75)),
    min       = tapply(y, grp, min),
    max       = tapply(y, grp, max),
    row.names = NULL
  )

  # Add skewness and kurtosis if moments package available
  if (requireNamespace("moments", quietly = TRUE)) {
    data_summary$skewness <- tapply(y, grp, moments::skewness)
    data_summary$kurtosis <- tapply(y, grp, moments::kurtosis)
  }

  # --- 9. Print Results ---
  if (verbose) {
    message(sprintf("\n--- Test Results ---"))
    message(sprintf("Test Selected: %s", test_used))
    message(sprintf("Test statistic: %.4f", test_result$statistic))

    if (test_result$p.value < 0.001) {
      message(sprintf("P-value: < 0.001"))
    } else {
      message(sprintf("P-value: %.4f", test_result$p.value))
    }

    if (!is.null(effect_size_result)) {
      message(sprintf("%s", effect_size_result$interpretation))
    }

    if (test_result$p.value < test_alpha) {
      message(sprintf("Interpretation: Significant difference detected (p < %.2f) ***", test_alpha))
    } else {
      message(sprintf("Interpretation: No significant difference detected (p >= %.2f)", test_alpha))
    }

    if (!is.null(posthoc_result)) {
      message(sprintf("\n--- Post-Hoc Comparisons (%s) ---", unique(posthoc_result$posthoc_test)))
      sig_comparisons <- sum(posthoc_result$significant, na.rm = TRUE)
      total_comparisons <- nrow(posthoc_result)
      message(sprintf("Significant pairwise differences: %d out of %d", sig_comparisons, total_comparisons))
    }
  }

  # --- 10. Return Object ---
  result <- list(
    test_used      = test_used,
    test_result    = test_result,
    effect_size    = effect_size_result,
    posthoc_result = posthoc_result,
    assumptions    = assumptions,
    parametric     = parametric_used,
    data_summary   = data_summary,
    outcome        = outcome,
    group_var      = group,
    paired         = paired,
    test_type      = test_type,
    test_alpha     = test_alpha,
    warnings       = warnings_list,
    parameters     = list(
      test_type             = test_type,
      normality_method      = normality_method,
      alpha_normality       = alpha_normality,
      alpha_variance        = alpha_variance,
      test_alpha            = test_alpha,
      min_n_threshold       = min_n_threshold,
      calculate_effect_size = calculate_effect_size,
      perform_posthoc       = perform_posthoc,
      p_adjust_method       = p_adjust_method,
      group_order           = group_order,
      paired                = paired
    ),
    raw_data       = tibble::tibble(outcome = y, group = grp)
  )

  class(result) <- "auto_compare"
  return(result)
}

#' Print method for auto_compare objects
#' @param x An object of class "auto_compare".
#' @param ... Further arguments passed to or from other methods.
#' @export
print.auto_compare <- function(x, ...) {
  cat("\n=== Automatic Statistical Comparison ===\n\n")
  cat(sprintf("Test Used: %s\n", x$test_used))
  cat(sprintf("Parametric: %s\n\n", x$parametric))

  cat("Test Results:\n")
  cat(sprintf("  Statistic: %.4f\n", x$test_result$statistic))

  if (x$test_result$p.value < 0.001) {
    cat("  P-value: < 0.001 ***\n")
  } else {
    cat(sprintf("  P-value: %.4f", x$test_result$p.value))
    if (x$test_result$p.value < 0.001) cat(" ***")
    else if (x$test_result$p.value < 0.01) cat(" **")
    else if (x$test_result$p.value < 0.05) cat(" *")
    cat("\n")
  }

  if (!is.null(x$effect_size)) {
    cat(sprintf("  %s\n", x$effect_size$interpretation))
  }

  cat("\n")

  cat("Data Summary:\n")
  print(x$data_summary[, c("group", "n", "mean", "sd", "median")], row.names = FALSE)

  if (!is.null(x$posthoc_result)) {
    cat(sprintf("\n\nPost-Hoc Tests (%s):\n", unique(x$posthoc_result$posthoc_test)))
    sig_only <- x$posthoc_result[x$posthoc_result$significant %in% TRUE,
                                 c("group1", "group2", "p.adj", "p.adj.signif")]
    if (nrow(sig_only) > 0) {
      cat("Significant comparisons:\n")
      print(sig_only, row.names = FALSE)
    } else {
      cat("No significant pairwise differences found.\n")
    }
  }

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat(sprintf("  - %s\n", w))
    }
  }
}

#' Summary method for auto_compare objects
#' @param object An object of class "auto_compare".
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.auto_compare <- function(object, ...) {
  cat("\n=== Detailed Statistical Comparison Summary ===\n\n")
  cat(sprintf("Outcome Variable: %s\n", object$outcome))
  cat(sprintf("Grouping Variable: %s\n", object$group_var))
  cat(sprintf("Paired Analysis: %s\n", object$paired))
  cat(sprintf("Test Strategy: %s (alpha = %.2f)\n", object$test_type, object$test_alpha))
  cat(sprintf("Test Selected: %s\n\n", object$test_used))

  cat("--- Descriptive Statistics ---\n")
  print(object$data_summary, row.names = FALSE)

  cat("\n--- Assumption Checks ---\n")
  print(object$assumptions, row.names = FALSE)

  cat("\n--- Statistical Test Results ---\n")
  print(object$test_result)

  if (!is.null(object$effect_size)) {
    cat("\n--- Effect Size ---\n")
    cat(sprintf("Metric: %s\n", object$effect_size$metric))
    cat(sprintf("Estimate: %.3f", object$effect_size$estimate))
    if (!is.na(object$effect_size$ci_low)) {
      cat(sprintf(" [95%% CI: %.3f, %.3f]", object$effect_size$ci_low, object$effect_size$ci_high))
    }
    cat("\n")
    cat(sprintf("Magnitude: %s\n", object$effect_size$magnitude))
  }

  if (!is.null(object$posthoc_result)) {
    cat(sprintf("\n--- Post-Hoc Pairwise Comparisons (%s) ---\n", unique(object$posthoc_result$posthoc_test)))
    cat(sprintf("P-value adjustment method: %s\n\n", object$parameters$p_adjust_method))

    # Build list of columns to display based on what's available
    cols_to_show <- c("group1", "group2", "n1", "n2", "statistic", "p", "p.adj", "p.adj.signif", "interpretation")

    # Filter out any columns that might not be in the result (though they should be)
    cols_to_show <- cols_to_show[cols_to_show %in% names(object$posthoc_result)]

    print(object$posthoc_result[, cols_to_show], row.names = FALSE)
  }

  cat("\n--- Interpretation ---\n")
  if (object$test_result$p.value < object$test_alpha) {
    cat(sprintf("Result: SIGNIFICANT difference detected between groups (p < %.2f).\n", object$test_alpha))
    if (!is.null(object$effect_size)) {
      cat(sprintf("The effect size is %s (%s = %.3f).\n",
                  tolower(object$effect_size$magnitude),
                  object$effect_size$metric,
                  object$effect_size$estimate))
    }
  } else {
    cat(sprintf("Result: NO significant difference detected between groups (p >= %.2f).\n", object$test_alpha))
  }

  if (length(object$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in object$warnings) {
      cat(sprintf("  - %s\n", w))
    }
  }

  invisible(object)
}
