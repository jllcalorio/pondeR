#' Automatic Statistical Comparison with Comprehensive Analysis
#'
#' @title Automatic Statistical Comparison with Comprehensive Analysis
#'
#' @description
#' Performs an automatic comparison of a numeric outcome variable across two or more
#' groups. Intelligently selects the appropriate statistical test (parametric or
#' non-parametric) based on checks for normality (of data or residuals) and homogeneity
#' of variances, calculates effect sizes, and performs post-hoc tests when appropriate.
#'
#' When \code{outcome} is a character vector of length > 1 and \code{summary_table = TRUE},
#' an additional \code{$summary_table} element is returned alongside the per-outcome list,
#' providing a single tidy data frame that consolidates key results across all outcomes —
#' useful for high-throughput screening (e.g., metabolomics, proteomics).
#'
#' When \code{subgroup} is also specified and \code{summary_table = TRUE}, additional
#' named elements are appended to the returned list, one per subgroup variable, each
#' containing a named list of per-level summary tables. For example, with
#' \code{subgroup = "Gender"}, the result will include \code{$summary_table_Gender},
#' and inside it \code{$summary_table_Gender$Male} and \code{$summary_table_Gender$Female}.
#'
#' @details
#' The function follows a decision-making process to determine the most appropriate test:
#'
#' \strong{Test Family Selection:}
#' If \code{paired = FALSE}, it compares independent groups. If there are 2 groups, it
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
#' For significant results with >2 groups, appropriate post-hoc tests are automatically
#' performed:
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
#' \strong{Non-Parametric Interpretation (Stochastic Dominance):}
#' For non-parametric tests, group comparisons are interpreted using stochastic
#' dominance rather than median comparisons. The \code{interpretation} column in
#' \code{posthoc_result} (and the 2-group case) reads "X tends to have larger/smaller
#' values than Y", reflecting that the Wilcoxon / Mann-Whitney family tests whether
#' one group's values tend to be systematically higher or lower than another's.
#'
#' @param x A data frame containing the variables for analysis.
#' @param outcome A character string or character vector specifying the name(s) of the
#'   numeric outcome variable(s). When a vector is supplied, the function loops over each
#'   outcome and returns a named list of \code{"run_diff"} objects. If
#'   \code{summary_table = TRUE}, an additional \code{$summary_table} element is appended
#'   to the returned list (see \code{summary_table}).
#' @param group A character string specifying the name of the grouping factor variable.
#' @param paired A logical indicating whether the observations are paired. Default is
#'   \code{FALSE}.
#' @param test_type A character string specifying the test strategy. Must be one of
#'   \code{"auto"} (default), \code{"parametric"}, or \code{"nonparametric"}.
#' @param normality_method Method for assessing normality: \code{"shapiro"}
#'   (Shapiro-Wilk test), \code{"lilliefors"} (Lilliefors (Kolmogorov-Smirnov) test),
#'   or \code{"auto"} (uses Shapiro-Wilk for n < 50, Lilliefors for n >= 50).
#'   Default is \code{"auto"}.
#' @param alpha_normality Significance level for normality tests. Default is 0.05.
#' @param alpha_variance Significance level for variance homogeneity tests. Default is
#'   0.05.
#' @param alpha_sphericity Significance level for sphericity tests. Default is 0.05.
#' @param type_sosquares Integer specifying the type of sums of squares for the
#'   ANOVA. Accepted values are \code{1}, \code{2} (default), or \code{3}.
#'   Type 2 yields identical results to type 1 for balanced designs but is more
#'   appropriate for unbalanced designs. Type 3 matches output from commercial
#'   software such as SPSS. Only used for RM and mixed ANOVA designs. See
#'   \code{\link[rstatix]{anova_test}} for details.
#' @param rm_effect_size Character vector specifying which effect size(s) to
#'   compute for RM and mixed ANOVA designs. Accepted values are \code{"ges"}
#'   (generalized eta-squared), \code{"pes"} (partial eta-squared, default), or
#'   \code{c("ges", "pes")} for both. When both are requested, \code{"ges"} is
#'   used as the primary metric reported in \code{$effect_size}. Only used for
#'   RM and mixed ANOVA designs. See \code{\link[rstatix]{anova_test}}.
#' @param correction Character string specifying the sphericity correction
#'   method applied to within-subject p-values when sphericity is violated.
#'   Passed to \code{\link[rstatix]{get_anova_table}}. Accepted values are
#'   \code{"auto"} (default; applies Greenhouse-Geisser when sphericity is
#'   violated), \code{"GG"} (always apply Greenhouse-Geisser), \code{"HF"}
#'   (always apply Huynh-Feldt), or \code{"none"} (uncorrected). Only used for
#'   RM and mixed ANOVA designs.
#' @param test_alpha Significance level for the final statistical test. Default is 0.05.
#' @param min_n_threshold Minimum sample size required per group. Default is 3.
#' @param calculate_effect_size Logical indicating whether to calculate effect sizes.
#'   Default is \code{TRUE}.
#' @param perform_posthoc Logical indicating whether to perform post-hoc tests for >2
#'   groups. Default is \code{TRUE}.
#' @param p_adjust_method The p-value adjustment method for post-hoc comparisons:
#'   \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"},
#'   \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}. Default is \code{"BH"}.
#' @param group_order Character vector specifying the order of groups. If \code{NULL},
#'   uses factor levels in alphabetical order.
#' @param subgroup Character vector or \code{NULL}. Column name(s) in \code{x} to use for
#'   subgroup analysis. Each must be a categorical column (factor or character) in
#'   \code{x}. When specified, the full analysis is repeated independently for each level
#'   of each subgroup column. Results are returned under \code{$subgroup_analysis} in
#'   the output, nested as \code{subgroup_var -> level -> run_diff result}.
#'   Default is \code{NULL}.
#' @param summary_table Logical. Only relevant when \code{outcome} is a character vector
#'   of length > 1. If \code{TRUE}, a tidy summary data frame is constructed across all
#'   outcomes and appended to the returned list as \code{$summary_table}. When
#'   \code{subgroup} is also specified, additional per-subgroup summary table lists are
#'   appended as \code{$summary_table_<subgroup_var>}, each containing one data frame per
#'   level (e.g., \code{$summary_table_Gender$Male}). The data frame contains one row per
#'   outcome with the following columns:
#'   \code{outcome}, \code{test_used}, \code{statistic}, \code{df}, \code{p_value},
#'   \code{effect_size_metric}, \code{effect_size_estimate}, \code{effect_size_ci_low},
#'   \code{effect_size_ci_high}, \code{effect_size_magnitude}, \code{significant},
#'   \code{n_significant_posthoc} (number of significant post-hoc pairs; \code{NA} for
#'   2-group comparisons), and \code{posthoc_pairs} (total post-hoc pairs tested;
#'   \code{NA} for 2-group comparisons). Default is \code{FALSE}.
#' @param verbose A logical indicating whether to print detailed messages. Default is
#'   \code{TRUE}.
#' @param num_cores Integer specifying the number of cores to use for parallel
#'   processing, or the character string \code{"max"}. If \code{"max"}, the function
#'   uses the number of available cores minus 2 (to reserve system resources), with a
#'   minimum of 1. Default is 1. Supports both Windows and POSIX (Mac/Linux) systems.
#'
#' @return When \code{outcome} is a single string, returns a list of class
#'   \code{"run_diff"} containing:
#'   \item{test_used}{Character string of the statistical test performed.}
#'   \item{test_result}{List containing results of the statistical test.}
#'   \item{effect_size}{List containing effect size estimate, confidence interval,
#'     magnitude, and interpretation.}
#'   \item{posthoc_result}{Data frame of post-hoc pairwise comparisons (\code{NULL} if
#'     not applicable). The \code{interpretation} column uses stochastic dominance
#'     framing for non-parametric tests ("X tends to have larger/smaller values than Y")
#'     and mean-based framing for parametric tests.}
#'   \item{assumptions}{Data frame summarizing assumption check results.}
#'   \item{data_summary}{Data frame with comprehensive descriptive statistics per group.}
#'   \item{outcome}{Name of the outcome variable.}
#'   \item{group_var}{Name of the group variable.}
#'   \item{paired}{The paired setting used.}
#'   \item{test_type}{The test_type strategy used.}
#'   \item{test_alpha}{The significance level used.}
#'   \item{subgroup_analysis}{Named list of subgroup results, nested as
#'     \code{subgroup_var -> level -> run_diff result}. \code{NULL} if no subgroup
#'     specified.}
#'   \item{warnings}{Character vector of any warnings generated.}
#'   \item{parameters}{List of all parameters used in the analysis.}
#'   \item{raw_data}{Data frame containing the variables used.}
#'
#'   When \code{outcome} is a character vector of length > 1, returns a named list where
#'   each element is a \code{"run_diff"} object for the corresponding outcome variable.
#'   If \code{summary_table = TRUE}, the list also contains:
#'   \itemize{
#'     \item \code{$summary_table} — overall tidy summary across all outcomes.
#'     \item \code{$summary_table_<sg>} (one per subgroup variable \code{sg}) — a named
#'       list of per-level summary tables, e.g.
#'       \code{$summary_table_Gender$Male}, \code{$summary_table_Gender$Female}.
#'   }
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
#' @references Bartlett, M. S. (1937). \emph{Properties of sufficiency and statistical tests}. Proceedings of the Royal Society of London Series A 160, 268-282. doi:10.1098/rspa.1937.0109.
#' @references Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988) \emph{The New S Language}. Wadsworth & Brooks/Cole.
#' @references Benjamini, Y., and Hochberg, Y. (1995). \emph{Controlling the false discovery rate: a practical and powerful approach to multiple testing}. Journal of the Royal Statistical Society Series B, 57, 289-300. doi:10.1111/j.2517-6161.1995.tb02031.x.
#' @references Benjamini, Y., and Yekutieli, D. (2001). \emph{The control of the false discovery rate in multiple testing under dependency}. Annals of Statistics, 29, 1165-1188. doi:10.1214/aos/1013699998.
#' @references Ben-Shachar M, Ludecke D, Makowski D (2020). \emph{effectsize: Estimation of Effect Size Indices and Standardized Parameters}. Journal of Open Source Software, 5(56), 2815. doi: 10.21105/joss.02815
#' @references Chambers, J. M., Freeny, A and Heiberger, R. M. (1992) \emph{Analysis of variance; designed experiments}. Chapter 5 of Statistical Models in S eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' @references Chambers, J. M. and Hastie, T. J. (1992) \emph{Statistical Models in S}. Wadsworth & Brooks/Cole.
#' @references Cliff, N. (1993). \emph{Dominance statistics: Ordinal analyses to answer ordinal questions}. Psychological bulletin, 114(3), 494.
#' @references Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences (2nd Ed.)}. New York: Routledge.
#' @references Cureton, E. E. (1956). \emph{Rank-biserial correlation}. Psychometrika, 21(3), 287-290.
#' @references Dallal, G.E. and Wilkinson, L. (1986): \emph{An analytic approximation to the distribution of Lilliefors' test for normality}. The American Statistician, 40, 294-296.
#' @references David F. Bauer (1972). \emph{Constructing confidence sets using rank statistics}. Journal of the American Statistical Association 67, 687-690. doi:10.1080/01621459.1972.10481279.
#' @references Delacre, M., Lakens, D., Ley, C., Liu, L., & Leys, C. (2021, May 7). \emph{Why Hedges' g*s based on the non-pooled standard deviation should be reported with Welch's t-test}. doi:10.31234/osf.io/tu6mp
#' @references Dunn, O. J. (1964) \emph{Multiple comparisons using rank sums Technometrics}, 6(3):241-252.
#' @references Eric Langford (2006) \emph{Quartiles in Elementary Statistics}, Journal of Statistics Education 14 3; doi:10.1080/10691898.2006.11910589
#' @references Fox, J. (2016) \emph{Applied Regression Analysis and Generalized Linear Models}, Third Edition. Sage.
#' @references Fox J, Weisberg S (2019). \emph{_An R Companion to Applied Regression_}, Third edition. Sage, Thousand Oaks CA. <https://www.john-fox.ca/Companion/>.
#' @references Glass, G. V. (1965). \emph{A ranking variable analogue of biserial correlation: Implications for short-cut item analysis}. Journal of Educational Measurement, 2(1), 91-95.
#' @references Gross J, Ligges U (2015). \emph{_nortest: Tests for Normality_}. doi:10.32614/CRAN.package.nortest <https://doi.org/10.32614/CRAN.package.nortest>, R package version 1.0-4, <https://CRAN.R-project.org/package=nortest>.
#' @references Hedges, L. V. & Olkin, I. (1985). \emph{Statistical methods for meta-analysis}. Orlando, FL: Academic Press.
#' @references Hochberg, Y. (1988). \emph{A sharper Bonferroni procedure for multiple tests of significance}. Biometrika, 75, 800-803. doi:10.2307/2336325
#' @references Holm, S. (1979). \emph{A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics}, 6, 65-70. \url{https://www.jstor.org/stable/4615733}.
#' @references Hommel, G. (1988). \emph{A stagewise rejective multiple test procedure based on a modified Bonferroni test}. Biometrika, 75, 383-386. doi:10.2307/2336190.
#' @references Hunter, J. E., & Schmidt, F. L. (2004). \emph{Methods of meta-analysis: Correcting error and bias in research findings}. Sage.
#' @references Hyndman, R. J. and Fan, Y. (1996) \emph{Sample quantiles in statistical packages}, American Statistician 50, 361-365. doi:10.2307/2684934.
#' @references Kassambara A (2023). \emph{_rstatix: Pipe-Friendly Framework for Basic Statistical Tests_}. doi:10.32614/CRAN.package.rstatix <https://doi.org/10.32614/CRAN.package.rstatix>, R package version 0.7.2, <https://CRAN.R-project.org/package=rstatix>.
#' @references Kelley, T. (1935) \emph{An unbiased correlation ratio measure}. Proceedings of the National Academy of Sciences. 21(9). 554-559.
#' @references Kerby, D. S. (2014). \emph{The simple difference formula: An approach to teaching nonparametric correlation}. Comprehensive Psychology, 3, 11-IT.
#' @references King, B. M., & Minium, E. W. (2008). \emph{Statistical reasoning in the behavioral sciences}. John Wiley & Sons Inc.
#' @references Komsta L, Novomestky F (2022). \emph{_moments: Moments, Cumulants, Skewness, Kurtosis and Related Tests_}. doi:10.32614/CRAN.package.moments <https://doi.org/10.32614/CRAN.package.moments>, R package version 0.14.1, <https://CRAN.R-project.org/package=moments>.
#' @references Makkonen, L. and Pajari, M. (2014) \emph{Defining Sample Quantiles by the True Rank Probability}, Journal of Probability and Statistics; Hindawi Publ.Corp. doi:10.1155/2014/326579
#' @references Myles Hollander and Douglas A. Wolfe (1973). \emph{Nonparametric Statistical Methods}. New York: John Wiley & Sons. Pages 27-33 (one-sample), 68-75 (two-sample). Or second edition (1999).
#' @references Myles Hollander and Douglas A. Wolfe (1973), \emph{Nonparametric Statistical Methods}. New York: John Wiley & Sons. Pages 115-120.
#' @references Muller K, Wickham H (2025). \emph{_tibble: Simple Data Frames_}. doi:10.32614/CRAN.package.tibble <https://doi.org/10.32614/CRAN.package.tibble>, R package version 3.3.0, <https://CRAN.R-project.org/package=tibble>.
#' @references Olejnik, S., & Algina, J. (2003). \emph{Generalized eta and omega squared statistics: measures of effect size for some common research designs}. Psychological methods, 8(4), 434.
#' @references Patrick Royston (1982). \emph{Algorithm AS 181: The W test for Normality}. Applied Statistics, 31, 176-180. doi:10.2307/2347986.
#' @references Patrick Royston (1982). \emph{An extension of Shapiro and Wilk's W test for normality to large samples}. Applied Statistics, 31, 115-124. doi:10.2307/2347973.
#' @references Patrick Royston (1995). \emph{Remark AS R94: A remark on Algorithm AS 181: The W test for normality}. Applied Statistics, 44, 547-551. doi:10.2307/2986146.
#' @references Ruxton, G.D., and Beauchamp, G. (2008) \emph{'Time for some a priori thinking about post hoc testing'}, Behavioral Ecology, 19(3), pp. 690-693. doi: 10.1093/beheco/arn020.
#' @references Sangseok Lee, Dong Kyu Lee. \emph{What is the proper way to apply the multiple comparison test?}. Korean J Anesthesiol. 2018;71(5):353-360.
#' @references Sarkar, S. (1998). \emph{Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture}. Annals of Statistics, 26, 494-504. doi:10.1214/aos/1013699998.
#' @references Sarkar, S., and Chang, C. K. (1997). \emph{The Simes method for multiple hypothesis testing with positively dependent test statistics}. Journal of the American Statistical Association, 92, 1601-1608. doi:10.2307/2965431.
#' @references Shaffer, J. P. (1995). \emph{Multiple hypothesis testing. Annual Review of Psychology}, 46, 561-584. doi:10.1146/annurev.ps.46.020195.003021.
#' @references Steiger, J. H. (2004). \emph{Beyond the F test: Effect size confidence intervals and tests of close fit in the analysis of variance and contrast analysis}. Psychological Methods, 9, 164-182.
#' @references Stephens, M.A. (1974): \emph{EDF statistics for goodness of fit and some comparisons}. Journal of the American Statistical Association, 69, 730-737.
#' @references Thode Jr., H.C. (2002): \emph{Testing for Normality}. Marcel Dekker, New York.
#' @references Tomczak, M., & Tomczak, E. (2014). \emph{The need to report effect size estimates revisited}. An overview of some recommended measures of effect size.
#' @references Wickham H, Averick M, Bryan J, Chang W, McGowan LD, Francois R, Grolemund G, Hayes A, Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Muller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019). \emph{"Welcome to the tidyverse."} _Journal of Open Source Software_, *4*(43), 1686. doi:10.21105/joss.01686 <https://doi.org/10.21105/joss.01686>.
#' @references Wicklin, R. (2017) \emph{Sample quantiles: A comparison of 9 definitions}; SAS Blog. https://blogs.sas.com/content/iml/2017/05/24/definitions-sample-quantiles.htm
#' @references Wright, S. P. (1992). \emph{Adjusted P-values for simultaneous inference}. Biometrics, 48, 1005-1013. doi:10.2307/2532694.
#' @references R Core Team (2025). \emph{_R: A Language and Environment for Statistical Computing_}. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/>.
#'
#' @author John Lennon L. Calorio
#'
#' @examples
#' \dontrun{
#' # --- Basic Usage with 3+ Groups ---
#' data(iris)
#' res_iris <- run_diff(x = iris, outcome = "Sepal.Length", group = "Species")
#' print(res_iris)
#' summary(res_iris)
#'
#' # --- Two-Sample Independent Comparison ---
#' data(mtcars)
#' mtcars$am <- factor(mtcars$am, labels = c("automatic", "manual"))
#' res_mtcars <- run_diff(x = mtcars, outcome = "mpg", group = "am")
#' print(res_mtcars)
#'
#' # --- Paired-Sample Comparison ---
#' data(sleep)
#' res_sleep <- run_diff(x = sleep, outcome = "extra", group = "group", paired = TRUE)
#' summary(res_sleep)
#'
#' # --- Multi-outcome with summary table (e.g., metabolomics screening) ---
#' data(iris)
#' outcomes <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
#' res_multi <- run_diff(
#'   x             = iris,
#'   outcome       = outcomes,
#'   group         = "Species",
#'   summary_table = TRUE,
#'   verbose       = FALSE
#' )
#' # Access per-outcome results as usual
#' print(res_multi$Sepal.Length)
#' # Filter features with significant omnibus p-value
#' sig_features <- res_multi$summary_table[res_multi$summary_table$significant, ]
#' print(sig_features)
#'
#' # --- Multi-outcome + subgroup summary tables ---
#' # Suppose df has columns: metabolites..., "Disease", "Gender"
#' res_sg <- run_diff(
#'   x             = df,
#'   outcome       = metabolite_cols,
#'   group         = "Disease",
#'   subgroup      = "Gender",
#'   summary_table = TRUE,
#'   verbose       = FALSE
#' )
#' # Overall summary table
#' res_sg$summary_table
#' # Per-subgroup summary tables
#' res_sg$summary_table_Gender$Male
#' res_sg$summary_table_Gender$Female
#'
#' # --- Forcing a Non-Parametric Test ---
#' data(iris)
#' res_iris_np <- run_diff(
#'   x         = iris,
#'   outcome   = "Sepal.Length",
#'   group     = "Species",
#'   test_type = "nonparametric"
#' )
#' print(res_iris_np)
#'
#' # --- Forcing a Parametric Test & Custom P-Value Adjustment ---
#' data(ToothGrowth)
#' ToothGrowth$dose <- as.factor(ToothGrowth$dose)
#' res_tooth <- run_diff(
#'   x               = ToothGrowth,
#'   outcome         = "len",
#'   group           = "dose",
#'   test_type       = "parametric",
#'   p_adjust_method = "bonferroni"
#' )
#' print(res_tooth$posthoc_result)
#'
#' # --- Custom Group Ordering ---
#' data(PlantGrowth)
#' res_plant <- run_diff(
#'   x           = PlantGrowth,
#'   outcome     = "weight",
#'   group       = "group",
#'   group_order = c("trt1", "ctrl", "trt2"),
#'   verbose     = FALSE
#' )
#' print(res_plant$data_summary)
#' 
#' # --- One-way Repeated Measures ANOVA ---
#' # Does HeartRate change across time points (Pre, During, Post)?
#' res_rm1 <- run_diff(
#'   x          = df_long,
#'   outcome    = "HeartRate",
#'   within     = "Time",
#'   subject_id = "ID",
#'   verbose    = FALSE
#' )
#' print(res_rm1)
#' summary(res_rm1)
#'
#' # --- Two-way Repeated Measures ANOVA ---
#' # Does an outcome change across two within-subject factors?
#' # Requires a dataset where subjects are measured under all combinations
#' # of two within-subject factors (e.g., Time × Condition).
#' # Create a minimal example with two within factors:
#' set.seed(1)
#' df_2rm <- data.frame(
#'   id        = rep(paste0("S", 1:10), each = 6),
#'   time      = rep(rep(c("T1", "T2", "T3"), each = 2), 10),
#'   condition = rep(c("A", "B"), times = 30),
#'   score     = rnorm(60, mean = 50, sd = 10)
#' )
#' res_rm2 <- run_diff(
#'   x          = df_2rm,
#'   outcome    = "score",
#'   within     = c("time", "condition"),
#'   subject_id = "id",
#'   verbose    = FALSE
#' )
#' print(res_rm2)
#' res_rm2$anova_table
#'
#' # --- Two-way Mixed ANOVA ---
#' # Does HeartRate differ by Group (between) and Time (within),
#' # and is there a Group x Time interaction?
#' res_mixed <- run_diff(
#'   x          = df_long,
#'   outcome    = "HeartRate",
#'   group      = "Group",
#'   within     = "Time",
#'   subject_id = "ID",
#'   verbose    = FALSE
#' )
#' print(res_mixed)
#' summary(res_mixed)
#' # Full ANOVA table with all effects and GG-corrected p-values if needed:
#' res_mixed$anova_table
#' # Significant post-hoc pairs:
#' res_mixed$posthoc_result[res_mixed$posthoc_result$significant %in% TRUE, ]
#' }
#'
#' @export
run_diff <- function(
    x,
    outcome,
    group                 = NULL,
    within                = NULL,
    subject_id            = NULL,
    paired                = FALSE,
    test_type             = c("auto", "parametric", "nonparametric"),
    normality_method      = c("auto", "shapiro", "lilliefors"),
    alpha_normality       = 0.05,
    alpha_variance        = 0.05,
    alpha_sphericity      = 0.05,
    type_sosquares        = 2,
    rm_effect_size        = "pes",
    correction            = "auto",
    test_alpha            = 0.05,
    min_n_threshold       = 3,
    calculate_effect_size = TRUE,
    perform_posthoc       = TRUE,
    p_adjust_method       = "BH",
    group_order           = NULL,
    subgroup              = NULL,
    summary_table         = FALSE,
    verbose               = TRUE,
    num_cores             = 1
) {

  # ---------------------------------------------------------------------------
  # Multi-outcome dispatch
  # ---------------------------------------------------------------------------
  if (length(outcome) > 1) {

    if (!is.logical(summary_table) || length(summary_table) != 1L)
      stop("'summary_table' must be a single logical value (TRUE or FALSE).")

    results <- lapply(outcome, function(oc) {
      tryCatch(
        run_diff(
          x                     = x,
          outcome               = oc,
          group                 = group,
          within                = within,
          subject_id            = subject_id,
          paired                = paired,
          test_type             = test_type,
          normality_method      = normality_method,
          alpha_normality       = alpha_normality,
          alpha_variance        = alpha_variance,
          alpha_sphericity      = alpha_sphericity,
          type_sosquares        = type_sosquares,
          rm_effect_size        = rm_effect_size,
          correction            = correction,
          test_alpha            = test_alpha,
          min_n_threshold       = min_n_threshold,
          calculate_effect_size = calculate_effect_size,
          perform_posthoc       = perform_posthoc,
          p_adjust_method       = p_adjust_method,
          group_order           = group_order,
          subgroup              = subgroup,
          summary_table         = FALSE,   # single-outcome call; no recursion
          verbose               = verbose,
          num_cores             = num_cores
        ),
        error = function(e) {
          warning(paste0("run_diff failed for outcome '", oc, "': ", e$message),
                  call. = FALSE)
          NULL
        }
      )
    })
    names(results) <- outcome

    # -------------------------------------------------------------------------
    # summary_table block
    # -------------------------------------------------------------------------
    if (summary_table) {

      # 1. Overall summary table (across all outcomes, no subgroup filter)
      results[["summary_table"]] <- .build_summary_table(results, outcome, test_alpha)

      # 2. Per-subgroup summary tables
      #    For each subgroup variable sg and each level lvl, extract the
      #    per-outcome subgroup result and build a separate summary table.
      #    Result is stored as results[[ "summary_table_<sg>" ]][[ lvl ]].
      if (!is.null(subgroup)) {
        for (sg in subgroup) {

          # Determine the levels that were analysed (same logic used inside
          # the single-outcome subgroup block)
          sg_levels <- if (is.factor(x[[sg]])) {
            levels(droplevels(as.factor(x[[sg]])))
          } else {
            sort(unique(na.omit(as.character(x[[sg]]))))
          }

          sg_tables <- list()   # will become results[[ "summary_table_<sg>" ]]

          for (lvl in sg_levels) {

            # For each outcome, pull the run_diff result that lives at
            #   results[[ outcome ]]$subgroup_analysis[[ sg ]][[ lvl ]]
            # and build a flat list keyed by outcome name so
            # .build_summary_table() can consume it directly.
            level_results <- stats::setNames(
              lapply(outcome, function(oc) {
                res <- results[[oc]]
                if (is.null(res)) return(NULL)
                res$subgroup_analysis[[sg]][[lvl]]   # NULL if missing
              }),
              outcome
            )

            sg_tables[[lvl]] <- .build_summary_table(level_results, outcome, test_alpha)
          }

          # Attach to the top-level results list under the name
          # "summary_table_<sg>", e.g. "summary_table_Gender"
          results[[ paste0("summary_table_", sg) ]] <- sg_tables
        }
      }
    }

    return(results)
  }

  # ===========================================================================
  # Single-outcome path  (unchanged from your original — paste your full
  # existing single-outcome code here as-is, starting from step 1)
  # ===========================================================================

  # --- 1. Input Validation and Preparation ---
  test_type        <- match.arg(test_type)
  normality_method <- match.arg(normality_method)
  warnings_list    <- character(0)

  if (!is.logical(summary_table) || length(summary_table) != 1L)
    stop("'summary_table' must be a single logical value (TRUE or FALSE).")

  # --- Subgroup Validation ---
  if (!is.null(subgroup)) {
    if (!is.character(subgroup))
      stop("'subgroup' must be a character vector of column names.")
    missing_sg <- setdiff(subgroup, names(x))
    if (length(missing_sg) > 0)
      stop(paste("Subgroup column(s) not found in data:", paste(missing_sg, collapse = ", ")))
    non_cat <- subgroup[sapply(subgroup, function(s) is.numeric(x[[s]]) && !is.factor(x[[s]]))]
    if (length(non_cat) > 0)
      stop(paste("Subgroup column(s) must be categorical (factor or character):",
                 paste(non_cat, collapse = ", ")))
  }

  # --- RM-specific parameter validation ---
  if (!type_sosquares %in% c(1L, 2L, 3L))
    stop("'type_sosquares' must be 1, 2, or 3.")
  type_sosquares <- as.integer(type_sosquares)

  valid_rm_es <- c("ges", "pes", "both")
  if (!all(rm_effect_size %in% c("ges", "pes")))
    stop("'rm_effect_size' must be \"ges\", \"pes\", or c(\"ges\", \"pes\") for both.")
  # Normalise: sort so ges always comes first when both requested
  rm_effect_size <- intersect(c("ges", "pes"), rm_effect_size)

  correction <- match.arg(correction, c("auto", "GG", "HF", "none"))

  # --- within / subject_id / mixed-ANOVA Validation ---
  is_rm_design    <- !is.null(within)
  is_mixed_design <- !is.null(within) && !is.null(group)

  if (is_rm_design) {
    if (!requireNamespace("rstatix", quietly = TRUE))
      stop("Package 'rstatix' is required for repeated-measures and mixed ANOVA. ",
           "Please install it with: install.packages('rstatix')")
    if (is.null(subject_id))
      stop("'subject_id' must be specified for repeated-measures or mixed-ANOVA designs.")
    if (!is.character(within) || length(within) < 1)
      stop("'within' must be a character vector of within-subject factor column name(s).")
    missing_w <- setdiff(within, names(x))
    if (length(missing_w) > 0)
      stop(paste("Within-subject column(s) not found in data:",
                 paste(missing_w, collapse = ", ")))
    if (!subject_id %in% names(x))
      stop(paste0("subject_id column '", subject_id, "' not found in data."))
  }

  avail_cores                         <- parallel::detectCores()
  if (is.na(avail_cores)) avail_cores <- 1

  if (identical(num_cores, "max")) {
    num_cores <- max(1L, avail_cores - 2L)
    if (verbose)
      message(sprintf("Parallel processing: Using 'max' setting (%d cores active, 2 reserved).",
                      num_cores))
  } else if (is.numeric(num_cores)) {
    if (num_cores < 1 || num_cores > avail_cores || num_cores %% 1 != 0)
      stop(sprintf("'num_cores' must be an integer between 1 and %d, or the string 'max'.",
                   avail_cores))
    num_cores <- as.integer(num_cores)
    if (num_cores > 1 && verbose)
      message(sprintf("Parallel processing: Using %d cores.", num_cores))
  } else {
    stop("Invalid 'num_cores' argument. Must be an integer or 'max'.")
  }

  # if (!is.data.frame(x))
  #   stop("'x' must be a data frame.")
  # if (!all(c(outcome, group) %in% names(x)))
  #   stop(paste0("Outcome variable '", outcome, "' or group variable '", group,
  #               "' not found in data."))

  # y   <- x[[outcome]]
  # grp <- as.factor(x[[group]])

  # if (!is.numeric(y))
  #   stop("Outcome variable must be numeric.")

  # if (!is.null(group_order)) {
  #   if (!all(group_order %in% levels(grp)))
  #     stop("group_order contains levels not present in the data.")
  #   grp <- factor(grp, levels = group_order)
  # }

  # complete_cases <- complete.cases(y, grp)
  # if (sum(!complete_cases) > 0) {
  #   msg           <- sprintf("Removing %d rows with missing values.", sum(!complete_cases))
  #   warnings_list <- c(warnings_list, msg)
  #   if (verbose) message(msg)
  # }
  # y   <- y[complete_cases]
  # grp <- droplevels(grp[complete_cases])

  # if (nlevels(grp) < 2)
  #   stop("Grouping factor must have at least 2 levels after removing missing values.")

  # groups   <- levels(grp)
  # n_groups <- length(groups)
  # group_ns <- table(grp)
  # min_n    <- min(group_ns)

  # if (verbose) {
  #   message("\n=== Auto Compare Analysis ===")
  #   message(sprintf("Outcome: %s | Group: %s | Paired: %s", outcome, group, paired))
  #   message(sprintf("Number of groups: %d", n_groups))
  #   message(sprintf("Sample sizes: %s", paste(names(group_ns), "=", group_ns, collapse = ", ")))
  #   message(sprintf("Test Type: %s", test_type))
  # }

  if (!is.data.frame(x))
    stop("'x' must be a data frame.")

  # For RM designs, group may be NULL (pure within) or a between factor.
  # Validate outcome always; validate group only when it is specified.
  if (!outcome %in% names(x))
    stop(paste0("Outcome variable '", outcome, "' not found in data."))
  if (!is.null(group) && !group %in% names(x))
    stop(paste0("Group variable '", group, "' not found in data."))

  # Classical path needs y / grp. For RM designs these are still extracted
  # so assumption checks and verbose output work, but grp may be a dummy
  # if group is NULL (pure RM). The RM block returns early before using grp
  # for testing, so a dummy factor is harmless.
  y <- x[[outcome]]
  if (!is.numeric(y))
    stop("Outcome variable must be numeric.")

  if (!is.null(group)) {
    grp <- as.factor(x[[group]])
    if (!is.null(group_order)) {
      if (!all(group_order %in% levels(grp)))
        stop("group_order contains levels not present in the data.")
      grp <- factor(grp, levels = group_order)
    }
    complete_cases <- complete.cases(y, grp)
    if (sum(!complete_cases) > 0) {
      msg           <- sprintf("Removing %d rows with missing values.", sum(!complete_cases))
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    }
    y   <- y[complete_cases]
    grp <- droplevels(grp[complete_cases])
    if (nlevels(grp) < 2)
      stop("Grouping factor must have at least 2 levels after removing missing values.")
    groups   <- levels(grp)
    n_groups <- length(groups)
    group_ns <- table(grp)
    min_n    <- min(group_ns)
  } else {
    # Pure RM design: create a dummy single-level grp so downstream code
    # that references grp/groups/n_groups doesn't break before the early return.
    grp      <- factor(rep("all", nrow(x)))
    groups   <- levels(grp)
    n_groups <- 1L
    group_ns <- table(grp)
    min_n    <- as.integer(nrow(x))
  }

  if (verbose && !is_rm_design) {
    message("\n=== Auto Compare Analysis ===")
    message(sprintf("Outcome: %s | Group: %s | Paired: %s", outcome, group, paired))
    message(sprintf("Number of groups: %d", n_groups))
    message(sprintf("Sample sizes: %s", paste(names(group_ns), "=", group_ns, collapse = ", ")))
    message(sprintf("Test Type: %s", test_type))
  }

  # --- 2. Initialize Assumption Checks and Determine Test Family ---
  assumptions <- data.frame(
    check = character(), result = character(),
    p_value = numeric(), decision = character(),
    stringsAsFactors = FALSE
  )

  if (!is_rm_design) {
    assumptions <- rbind(assumptions, data.frame(
      check    = "Minimum sample size",
      result   = sprintf("Min n = %d", min_n),
      p_value  = NA_real_,
      decision = ifelse(min_n >= min_n_threshold, "PASS", "WARNING"),
      stringsAsFactors = FALSE
    ))
    if (min_n < min_n_threshold) {
      msg           <- sprintf("Warning: Minimum sample size (%d) is below threshold (%d)",
                               min_n, min_n_threshold)
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    }
  }

  if (min_n < min_n_threshold) {
    msg           <- sprintf("Warning: Minimum sample size (%d) is below threshold (%d)",
                             min_n, min_n_threshold)
    warnings_list <- c(warnings_list, msg)
    if (verbose) message(msg)
  }

  # if (n_groups == 2) {
  #   test_family <- ifelse(paired, "paired_two_sample", "independent_two_sample")
  # } else if (n_groups > 2) {
  #   test_family <- ifelse(paired, "repeated_measures_anova", "independent_anova")
  # }

  # if (paired && n_groups == 2) {
  #   if (group_ns[1] != group_ns[2])
  #     stop("For paired tests, groups must have equal sample sizes after removing missing values.")
  # }
  if (is_rm_design) {
    has_between <- !is.null(group) &&
                   group %in% names(x) &&
                   length(unique(na.omit(as.character(x[[group]])))) > 1
    if (length(within) == 1 && !has_between) {
      test_family <- "one_way_rm_anova"
    } else if (length(within) >= 2 && !has_between) {
      test_family <- "two_way_rm_anova"
    } else if (length(within) >= 1 && has_between) {
      test_family <- "mixed_anova"
    } else {
      test_family <- "one_way_rm_anova"
    }
  } else if (n_groups == 2) {
    test_family <- ifelse(paired, "paired_two_sample", "independent_two_sample")
  } else if (n_groups > 2) {
    test_family <- ifelse(paired, "repeated_measures_anova", "independent_anova")
  }

  # --- 3. Helper Function: Check Normality ---
  # check_normality <- function(data_vector, group_label = NULL, method = normality_method) {
  #   n <- length(data_vector)
  #   if (n < 3)
  #     return(list(violated = TRUE, p_value = NA_real_, method = "Insufficient data"))
  #   use_method <- if (method == "auto") { if (n < 50) "shapiro" else "lilliefors" } else method
  #   if (use_method == "shapiro") {
  #     shapiro_result <- shapiro.test(data_vector)
  #     p_val          <- shapiro_result$p.value
  #     return(list(violated = p_val < alpha_normality, p_value = p_val,
  #                 method = "Shapiro-Wilk test",
  #                 label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
  #   } else {
  #     if (!requireNamespace("nortest", quietly = TRUE)) {
  #       warning("Package 'nortest' not available. Falling back to Shapiro-Wilk test.")
  #       shapiro_result <- shapiro.test(data_vector)
  #       p_val          <- shapiro_result$p.value
  #       return(list(violated = p_val < alpha_normality, p_value = p_val,
  #                   method = "Shapiro-Wilk test (fallback)",
  #                   label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
  #     }
  #     lillie_result <- nortest::lillie.test(data_vector)
  #     p_val         <- lillie_result$p.value
  #     return(list(violated = p_val < alpha_normality, p_value = p_val,
  #                 method = "Lilliefors (K-S) test",
  #                 label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
  #   }
  # }
  check_normality <- function(data_vector, group_label = NULL, method = normality_method) {
  n <- length(data_vector)
  if (n < 3)
    return(list(violated = TRUE, p_value = NA_real_, method = "Insufficient data",
                label = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))

  if (length(unique(data_vector)) < 2) {
    msg <- if (!is.null(group_label))
      sprintf("All values in group '%s' are identical; normality assumption violated.", group_label)
    else
      "All values are identical; normality assumption violated."
    warnings_list <<- c(warnings_list, msg)
    return(list(violated = TRUE, p_value = NA_real_,
                method = "Skipped (zero variance)",
                label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
  }

  use_method <- if (method == "auto") { if (n < 50) "shapiro" else "lilliefors" } else method
  if (use_method == "shapiro") {
    shapiro_result <- shapiro.test(data_vector)
    p_val          <- shapiro_result$p.value
    return(list(violated = p_val < alpha_normality, p_value = p_val,
                method = "Shapiro-Wilk test",
                label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
  } else {
    if (!requireNamespace("nortest", quietly = TRUE)) {
      warning("Package 'nortest' not available. Falling back to Shapiro-Wilk test.")
      shapiro_result <- shapiro.test(data_vector)
      p_val          <- shapiro_result$p.value
      return(list(violated = p_val < alpha_normality, p_value = p_val,
                  method = "Shapiro-Wilk test (fallback)",
                  label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
    }
    lillie_result <- nortest::lillie.test(data_vector)
    p_val         <- lillie_result$p.value
    return(list(violated = p_val < alpha_normality, p_value = p_val,
                method = "Lilliefors (K-S) test",
                label  = if (!is.null(group_label)) sprintf("Normality: %s", group_label) else "Normality"))
  }
}

  # --- 4. Assumption Checking ---
  normality_violated <- FALSE
  variance_violated  <- FALSE

  if (test_type == "auto" && !is_rm_design) {
    if (n_groups == 2) {
      run_norm_check <- function(g) check_normality(y[grp == g], g, normality_method)
      norm_results   <- list()
      if (num_cores > 1) {
        if (.Platform$OS.type == "windows") {
          cl <- parallel::makeCluster(num_cores)
          on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
          parallel::clusterExport(cl,
                                  varlist = c("y", "grp", "check_normality",
                                              "normality_method", "alpha_normality"),
                                  envir   = environment())
          parallel::clusterEvalQ(cl, { if (requireNamespace("nortest", quietly = TRUE)) library(nortest) })
          norm_results <- parallel::parLapply(cl, groups, run_norm_check)
        } else {
          norm_results <- parallel::mclapply(groups, run_norm_check, mc.cores = num_cores)
        }
      } else {
        norm_results <- lapply(groups, run_norm_check)
      }
      for (norm_check in norm_results) {
        if (inherits(norm_check, "try-error")) { if (verbose) warning("A parallel worker failed."); next }
        assumptions <- rbind(assumptions, data.frame(
          check = norm_check$label, result = norm_check$method,
          p_value = norm_check$p_value,
          decision = ifelse(norm_check$violated, "VIOLATED", "OK"),
          stringsAsFactors = FALSE))
        if (norm_check$violated) normality_violated <- TRUE
      }
    } else {
      fit        <- lm(y ~ grp)
      norm_check <- check_normality(fit$residuals, NULL, normality_method)
      norm_check$label <- "Normality of residuals"
      assumptions <- rbind(assumptions, data.frame(
        check = norm_check$label, result = norm_check$method,
        p_value = norm_check$p_value,
        decision = ifelse(norm_check$violated, "VIOLATED", "OK"),
        stringsAsFactors = FALSE))
      if (norm_check$violated) normality_violated <- TRUE
    }
    if (!paired) {
      variance_check_result <- tryCatch({
        if (requireNamespace("car", quietly = TRUE)) {
          car::leveneTest(y, grp)$`Pr(>F)`[1]
        } else {
          if (verbose) message("Package 'car' not found. Using Bartlett's test.")
          bartlett.test(x = y, g = grp)$p.value
        }
      }, error = function(e) {
        msg <- "Error in variance test. Assuming homogeneity."
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message(msg)
        1
      })
      p_val          <- variance_check_result
      decision       <- ifelse(p_val < alpha_variance, "VIOLATED", "OK")
      test_used_name <- if (requireNamespace("car", quietly = TRUE)) "Levene's test" else "Bartlett's test"
      assumptions <- rbind(assumptions, data.frame(
        check = "Homogeneity of variance", result = test_used_name,
        p_value = p_val, decision = decision, stringsAsFactors = FALSE))
      if (decision == "VIOLATED") variance_violated <- TRUE
    }
  }

  # --- 5. Select and Perform Test ---
  test_used       <- NULL
  test_result     <- NULL
  parametric_used <- TRUE
  fitted_model    <- NULL

  if (!is_rm_design) {
    if (test_type == "auto") {
      if (test_family == "independent_two_sample") {
        if (normality_violated) {
          test_used       <- "Mann-Whitney U test (Wilcoxon rank-sum)"
          # test_result     <- wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]])
          test_result <- withCallingHandlers(
              stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]]),
              warning = function(w) {
                if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
                  warnings_list <<- c(warnings_list, conditionMessage(w))
                invokeRestart("muffleWarning")
              }
            )
          parametric_used <- FALSE
        } else if (variance_violated) {
          test_used   <- "Welch's t-test"
          test_result <- t.test(x = y[grp == groups[1]], y = y[grp == groups[2]], var.equal = FALSE)
        } else {
          test_used   <- "Student's t-test"
          test_result <- t.test(x = y[grp == groups[1]], y = y[grp == groups[2]], var.equal = TRUE)
        }
      } else if (test_family == "paired_two_sample") {
        g1 <- y[grp == groups[1]]; g2 <- y[grp == groups[2]]
        if (normality_violated) {
          test_used <- "Wilcoxon signed-rank test"; test_result <- wilcox.test(g1, g2, paired = TRUE)
          parametric_used <- FALSE
        } else {
          test_used <- "Paired t-test"; test_result <- t.test(g1, g2, paired = TRUE)
        }
      } else if (test_family == "independent_anova") {
        if (normality_violated) {
          test_used       <- "Kruskal-Wallis test"
          test_result     <- kruskal.test(x = y, g = grp)
          parametric_used <- FALSE
        } else if (variance_violated) {
          test_used   <- "Welch's ANOVA"
          test_result <- oneway.test(formula(y ~ grp), var.equal = FALSE)
        } else {
          test_used    <- "One-way ANOVA"
          fitted_model <- aov(y ~ grp)
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
        test_result <- t.test(x = y[grp == groups[1]], y = y[grp == groups[2]], var.equal = TRUE)
      } else if (test_family == "paired_two_sample") {
        test_used   <- "Paired t-test (forced)"
        test_result <- t.test(y[grp == groups[1]], y[grp == groups[2]], paired = TRUE)
      } else if (test_family == "independent_anova") {
        test_used    <- "One-way ANOVA (forced)"
        fitted_model <- aov(y ~ grp)
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
        # test_result <- wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]])
        test_result <- withCallingHandlers(
          # stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]]),
          stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]]),
          warning = function(w) {
            if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
              warnings_list <<- c(warnings_list, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )
      } else if (test_family == "paired_two_sample") {
        test_used   <- "Wilcoxon signed-rank test (forced)"
        test_result <- wilcox.test(y[grp == groups[1]], y[grp == groups[2]], paired = TRUE)
      } else if (test_family == "independent_anova") {
        test_used   <- "Kruskal-Wallis test (forced)"
        test_result <- kruskal.test(x = y, g = grp)
      }
    }
  }

  # if (test_family == "repeated_measures_anova")
    # stop("Repeated measures ANOVA is not yet implemented.")
  # ===========================================================================
  # RM / Mixed ANOVA path — executed before the classical effect-size /
  # post-hoc blocks; returns early with a fully structured result list.
  # ===========================================================================
  if (test_family %in% c("one_way_rm_anova", "two_way_rm_anova", "mixed_anova")) {

    # ---- 0.  Build analysis data frame ----------------------------------------
    # Keep only the columns we need; drop rows with any NA in key columns.
    key_cols <- c(subject_id, within, outcome)
    if (test_family == "mixed_anova") key_cols <- c(key_cols, group)
    key_cols <- unique(key_cols)

    df_rm <- x[, key_cols, drop = FALSE]

    # Coerce within-factor(s) and between-factor to factor
    for (wf in within) df_rm[[wf]] <- as.factor(df_rm[[wf]])
    if (test_family == "mixed_anova") df_rm[[group]] <- as.factor(df_rm[[group]])
    df_rm[[subject_id]] <- as.factor(df_rm[[subject_id]])

    # Remove incomplete cases
    cc_rm <- stats::complete.cases(df_rm)
    if (sum(!cc_rm) > 0) {
      msg <- sprintf("Removing %d rows with missing values (RM/Mixed ANOVA).", sum(!cc_rm))
      warnings_list <- c(warnings_list, msg)
      if (verbose) message(msg)
    }
    df_rm <- df_rm[cc_rm, , drop = FALSE]

    # ---- 0b.  Safe syntactic name for outcome ---------------------------------
    # Column names with spaces/hyphens/special characters break formula parsing
    # inside rstatix::anova_test(). Use a guaranteed-safe internal name and
    # restore the original name in outputs.
    safe_outcome <- ".outcome_internal_"
    df_rm[[safe_outcome]] <- df_rm[[outcome]]
    # All subsequent anova_test / pairwise_t_test calls use safe_outcome.
    # We restore the original outcome name in data_summary at step 6.

    # ---- 1.  Assumption: Normality of residuals --------------------------------
    assumptions_rm <- data.frame(
      check = character(), result = character(),
      p_value = numeric(), decision = character(),
      stringsAsFactors = FALSE
    )

    # Fit a simple lm on the outcome to get residuals for normality check
    rhs <- if (test_family == "mixed_anova") {
      paste(c(group, within), collapse = " * ")
    } else {
      paste(within, collapse = " * ")
    }
    lm_formula <- stats::as.formula(paste(safe_outcome, "~", rhs))
    lm_fit     <- tryCatch(stats::lm(lm_formula, data = df_rm), error = function(e) NULL)

    if (!is.null(lm_fit)) {
      nc <- check_normality(lm_fit$residuals, NULL, normality_method)
      nc$label <- "Normality of residuals"
      assumptions_rm <- rbind(assumptions_rm, data.frame(
        check   = nc$label, result  = nc$method,
        p_value = nc$p_value,
        decision = ifelse(nc$violated, "VIOLATED", "OK"),
        stringsAsFactors = FALSE
      ))
      if (nc$violated) {
        msg <- "Normality of residuals violated; interpret RM/mixed ANOVA results with caution."
        warnings_list <- c(warnings_list, msg)
        if (verbose) message("Warning: ", msg)
      }
    }

    # ---- 2.  Run rstatix::anova_test ------------------------------------------
    # Build the formula for anova_test.
    # For mixed ANOVA:   outcome ~ between + Error(subject/within1*within2)
    # For RM ANOVA:      outcome ~ Error(subject/within1*within2)
    # rstatix::anova_test() accepts a 'wid' and 'within' / 'between' approach
    # which is cleaner and handles Mauchly automatically.

    if (!requireNamespace("rstatix", quietly = TRUE))
      stop("Package 'rstatix' is required for RM/mixed ANOVA.")

    anova_args <- list(
      data        = df_rm,
      dv          = safe_outcome,
      wid         = subject_id,
      within      = within,
      type        = type_sosquares,
      effect.size = rm_effect_size[1],   # primary pass; see below for "both"
      detailed    = TRUE
    )
    if (test_family == "mixed_anova") {
      anova_args[["between"]] <- group
    }

    aov_res <- tryCatch(
      do.call(rstatix::anova_test, anova_args),
      error = function(e) {
        stop("rstatix::anova_test failed: ", e$message)
      }
    )

    # Extract ANOVA table applying the user's sphericity correction choice
    aov_tbl <- as.data.frame(rstatix::get_anova_table(aov_res, correction = correction))

    # ---- Handle "both" effect sizes -----------------------------------------
    # rstatix::anova_test() can only compute one effect size at a time when the
    # internal branch is a strict if/else (pes vs ges). When the user requests
    # both, we run a second pass for the other metric and merge the column in.
    if (length(rm_effect_size) == 2) {
      second_es <- rm_effect_size[2]   # the one not yet in aov_tbl
      anova_args2              <- anova_args
      anova_args2$effect.size  <- second_es
      aov_res2 <- tryCatch(
        do.call(rstatix::anova_test, anova_args2),
        error = function(e) {
          msg <- paste0("Second effect size ('", second_es, "') pass failed: ", e$message)
          warnings_list <<- c(warnings_list, msg)
          if (verbose) message("Warning: ", msg)
          NULL
        }
      )
      if (!is.null(aov_res2)) {
        aov_tbl2 <- as.data.frame(rstatix::get_anova_table(aov_res2, correction = correction))
        # Merge second effect size column if not already present
        if (second_es %in% names(aov_tbl2) && !second_es %in% names(aov_tbl)) {
          # Match rows by Effect name to be safe
          idx <- match(aov_tbl$Effect, aov_tbl2$Effect)
          aov_tbl[[second_es]] <- aov_tbl2[[second_es]][idx]
        }
      }
    }

    # ---- 3.  Sphericity assumption (Mauchly's test) ----------------------------
    # aov_res$`Mauchly's Test for Sphericity` is present when a within factor
    # has > 2 levels.
    if (!is.null(aov_res$`Mauchly's Test for Sphericity`)) {
      mauchly_tbl <- as.data.frame(aov_res$`Mauchly's Test for Sphericity`)
      for (i in seq_len(nrow(mauchly_tbl))) {
        sph_p    <- mauchly_tbl$p[i]
        sph_eff  <- as.character(mauchly_tbl$Effect[i])
        sph_dec  <- ifelse(!is.na(sph_p) & sph_p < alpha_sphericity,
                           "VIOLATED", "OK")
        assumptions_rm <- rbind(assumptions_rm, data.frame(
          check    = paste0("Sphericity (Mauchly): ", sph_eff),
          result   = sprintf("W = %.4f", mauchly_tbl$W[i]),
          p_value  = sph_p,
          decision = sph_dec,
          stringsAsFactors = FALSE
        ))
        if (sph_dec == "VIOLATED") {
          msg <- sprintf(
            "Sphericity violated for '%s' (Mauchly p = %.4f). GG correction applied.",
            sph_eff, sph_p)
          warnings_list <- c(warnings_list, msg)
          if (verbose) message("Warning: ", msg)
        }
      }
      # Also record GG/HF corrections if present
      if (!is.null(aov_res$`Sphericity Corrections`)) {
        sph_cor <- as.data.frame(aov_res$`Sphericity Corrections`)
        for (i in seq_len(nrow(sph_cor))) {
          assumptions_rm <- rbind(assumptions_rm, data.frame(
            check    = paste0("GG epsilon: ", sph_cor$Effect[i]),
            result   = sprintf("eps = %.4f", sph_cor$GGe[i]),
            p_value  = sph_cor$`p[GG]`[i],
            decision = "INFO",
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    # ---- 4.  Effect size -------------------------------------------------------
    # ges (generalized eta-squared) is already in aov_tbl$ges
    # We also expose the primary/omnibus effect for the top-level effect_size slot.
    # For one-way RM: the single within effect.
    # For two-way RM / mixed: the interaction term (last row) as primary.
    primary_row <- aov_tbl[nrow(aov_tbl), ]   # interaction or last main effect

    # Primary effect size for the top-level $effect_size slot.
    # Uses rm_effect_size[1] (the first requested metric).
    primary_es_col <- rm_effect_size[1]
    es_rm <- if (primary_es_col %in% names(aov_tbl) &&
                 !is.na(primary_row[[primary_es_col]])) {
      es_val  <- as.numeric(primary_row[[primary_es_col]])
      es_label <- if (primary_es_col == "ges") "Generalized eta-squared (ges)"
                  else                         "Partial eta-squared (pes)"
      # Field (2013) thresholds apply to both ges and pes
      mag <- if (es_val >= 0.14) "large" else if (es_val >= 0.06) "medium" else "small"
      list(
        estimate  = es_val,
        ci_low    = NA_real_,
        ci_high   = NA_real_,
        magnitude = mag,
        metric    = es_label,
        interpretation = sprintf("Effect size: %s (%s = %.3f)", mag, primary_es_col, es_val)
      )
    } else NULL

    # ---- 5.  Post-hoc tests ----------------------------------------------------
    posthoc_rm <- NULL

    if (perform_posthoc && requireNamespace("rstatix", quietly = TRUE)) {

      sig_effects <- aov_tbl$Effect[!is.na(aov_tbl$p) & aov_tbl$p < test_alpha]

      if (length(sig_effects) > 0) {
        ph_list <- list()

        for (eff in sig_effects) {
          eff_factors    <- trimws(strsplit(eff, ":")[[1]])
          within_in_eff  <- intersect(eff_factors, within)
          between_in_eff <- if (!is.null(group)) intersect(eff_factors, group) else character(0)

          # Skip effects that contain neither a within nor a between factor
          if (length(within_in_eff) == 0 && length(between_in_eff) == 0) next

          ph_one <- tryCatch({
            if (length(within_in_eff) > 0) {
              # -----------------------------------------------------------
              # Within or interaction effect: use paired pairwise_t_test.
              # For interactions, group the data by the between factor first
              # so comparisons are made within each between-factor level.
              # -----------------------------------------------------------
              ph_formula <- stats::as.formula(
                paste(safe_outcome, "~", within_in_eff[1])
              )
              if (length(between_in_eff) > 0) {
                # Interaction: run pairwise test grouped by between factor
                rstatix::pairwise_t_test(
                  data            = df_rm,
                  formula         = ph_formula,
                  paired          = TRUE,
                  p.adjust.method = p_adjust_method,
                  detailed        = TRUE
                ) |> as.data.frame()
              } else {
                # Main within effect only
                rstatix::pairwise_t_test(
                  data            = df_rm,
                  formula         = ph_formula,
                  paired          = TRUE,
                  p.adjust.method = p_adjust_method,
                  detailed        = TRUE
                ) |> as.data.frame()
              }
            } else {
              # -----------------------------------------------------------
              # Pure between-subject effect: unpaired pairwise_t_test
              # -----------------------------------------------------------
              ph_formula <- stats::as.formula(
                paste(safe_outcome, "~", between_in_eff[1])
              )
              rstatix::pairwise_t_test(
                data            = df_rm,
                formula         = ph_formula,
                paired          = FALSE,
                p.adjust.method = p_adjust_method,
                detailed        = TRUE
              ) |> as.data.frame()
            }
          }, error = function(e) {
            msg <- paste0("Post-hoc failed for effect '", eff, "': ", e$message)
            warnings_list <<- c(warnings_list, msg)
            if (verbose) message(msg)
            NULL
          })

          if (!is.null(ph_one) && nrow(ph_one) > 0) {
            ph_one$effect       <- eff
            ph_one$posthoc_test <- if (length(within_in_eff) > 0)
              "Pairwise t-test (paired)" else "Pairwise t-test (independent)"
            ph_list[[eff]] <- ph_one
          }
        }

        if (length(ph_list) > 0) {
          # Bind all effects; fill missing columns with NA
          all_cols <- unique(unlist(lapply(ph_list, names)))
          ph_list  <- lapply(ph_list, function(df) {
            for (cn in setdiff(all_cols, names(df))) df[[cn]] <- NA
            df[, all_cols, drop = FALSE]
          })
          posthoc_rm           <- do.call(rbind, ph_list)
          rownames(posthoc_rm) <- NULL

          # ---- p.adj.signif and significant ----
          if ("p.adj" %in% names(posthoc_rm)) {
            posthoc_rm$significant <- posthoc_rm$p.adj < test_alpha
            posthoc_rm$p.adj.signif <- ifelse(
              posthoc_rm$p.adj < test_alpha,
              ifelse(posthoc_rm$p.adj < 0.001, "***",
              ifelse(posthoc_rm$p.adj < 0.01,  "**",
              ifelse(posthoc_rm$p.adj < 0.05,  "*",
                     sprintf("* (p<%.2f)", test_alpha)))),
              "ns"
            )
            posthoc_rm$p.adj.format <- ifelse(
              posthoc_rm$p.adj < 0.001, "p < 0.001",
              sprintf("p = %.3f", posthoc_rm$p.adj)
            )
          }
          if ("p" %in% names(posthoc_rm)) {
            posthoc_rm$p.format <- ifelse(
              is.na(posthoc_rm$p), NA_character_,
              ifelse(posthoc_rm$p < 0.001, "p < 0.001",
                     sprintf("p = %.3f", posthoc_rm$p))
            )
          }

          # ---- mean1 / mean2 from within[1] marginal means ----
          mean_tbl <- stats::setNames(
            as.numeric(tapply(df_rm[[safe_outcome]], df_rm[[within[1]]], mean, na.rm = TRUE)),
            levels(as.factor(df_rm[[within[1]]]))
          )
          posthoc_rm$mean1 <- mean_tbl[match(posthoc_rm$group1, names(mean_tbl))]
          posthoc_rm$mean2 <- mean_tbl[match(posthoc_rm$group2, names(mean_tbl))]

          # ---- interpretation ----
          posthoc_rm$interpretation <- ifelse(
            !posthoc_rm$significant %in% TRUE,
            "No significant difference",
            ifelse(
              !is.na(posthoc_rm$mean1) & !is.na(posthoc_rm$mean2) &
                posthoc_rm$mean1 < posthoc_rm$mean2,
              paste0(posthoc_rm$group1, " has lower mean than ", posthoc_rm$group2),
              paste0(posthoc_rm$group1, " has higher mean than ", posthoc_rm$group2)
            )
          )
        }
      }
    }

    # ---- 6.  Data summary (per-cell means) -------------------------------------
    # For RM: summarise by within factor(s); for mixed: by between × within
    summary_factors <- if (test_family == "mixed_anova") c(group, within) else within
    ds_rm <- tryCatch({
      agg <- stats::aggregate(
        stats::as.formula(paste(safe_outcome, "~",
                                paste(summary_factors, collapse = " + "))),
        data = df_rm, FUN = function(v)
          c(n      = length(v),
            mean   = mean(v),
            sd     = stats::sd(v),
            se     = stats::sd(v) / sqrt(length(v)),
            median = stats::median(v),
            q25    = stats::quantile(v, 0.25, names = FALSE),
            q75    = stats::quantile(v, 0.75, names = FALSE),
            min    = min(v),
            max    = max(v))
      )
      # aggregate with c() returns a matrix column; flatten
      stats_mat <- as.data.frame(agg[[safe_outcome]])
      names(stats_mat) <- c("n", "mean", "sd", "se", "median", "q25", "q75", "min", "max")
      cbind(agg[, summary_factors, drop = FALSE], stats_mat)
    }, error = function(e) {
      warning("Could not build data_summary for RM/mixed ANOVA: ", e$message)
      NULL
    })

    # ---- 7.  Verbose output ----------------------------------------------------
    if (verbose) {
      message("\n--- RM/Mixed ANOVA Results ---")
      message(sprintf("Design: %s", test_family))
      print(aov_tbl)
      if (!is.null(posthoc_rm) && nrow(posthoc_rm) > 0) {
        message(sprintf("\n--- Post-Hoc (Pairwise t-test, adj: %s) ---", p_adjust_method))
        sig_ph <- posthoc_rm[posthoc_rm$significant %in% TRUE, ]
        if (nrow(sig_ph) > 0) print(sig_ph) else message("No significant pairwise differences.")
      }
    }

    # ---- 8.  Construct test_result stub for compatibility ----------------------
    # The rest of the package (summary_table, print.run_diff) reads
    # test_result$statistic, test_result$p.value, test_result$parameter.
    # We expose the primary (interaction or last) effect's values.
    test_result_rm <- list(
      statistic = stats::setNames(as.numeric(primary_row$F), "F value"),
      parameter = stats::setNames(
        c(as.numeric(primary_row$DFn), as.numeric(primary_row$DFd)),
        c("DFn", "DFd")
      ),
      p.value = as.numeric(primary_row$p),
      method  = test_family,
      anova_table = aov_tbl   # full table accessible if needed
    )

    # ---- 9.  Return ------------------------------------------------------------
    result_rm <- list(
      test_used         = switch(test_family,
        one_way_rm_anova = "One-way Repeated Measures ANOVA",
        two_way_rm_anova = "Two-way Repeated Measures ANOVA",
        mixed_anova      = "Two-way Mixed ANOVA"
      ),
      test_result       = test_result_rm,
      effect_size       = es_rm,
      posthoc_result    = posthoc_rm,
      assumptions       = assumptions_rm,
      parametric        = TRUE,
      data_summary      = ds_rm,
      anova_table       = aov_tbl,          # full ANOVA table (all effects)
      sphericity        = aov_res$`Mauchly's Test for Sphericity`,
      sphericity_corrections = aov_res$`Sphericity Corrections`,
      outcome           = outcome,
      group_var         = group,
      within_vars       = within,
      subject_id_var    = subject_id,
      paired            = paired,
      test_type         = test_type,
      test_alpha        = test_alpha,
      subgroup_analysis = NULL,
      warnings          = warnings_list,
      parameters        = list(
        test_type             = test_type,
        normality_method      = normality_method,
        alpha_normality       = alpha_normality,
        alpha_variance        = alpha_variance,
        alpha_sphericity      = alpha_sphericity,
        type_sosquares        = type_sosquares,
        rm_effect_size        = rm_effect_size,
        correction            = correction,
        test_alpha            = test_alpha,
        min_n_threshold       = min_n_threshold,
        calculate_effect_size = calculate_effect_size,
        perform_posthoc       = perform_posthoc,
        p_adjust_method       = p_adjust_method,
        group_order           = group_order,
        subgroup              = subgroup,
        within                = within,
        subject_id            = subject_id,
        paired                = paired
      ),
      raw_data = df_rm
    )
    class(result_rm) <- "run_diff"
    return(result_rm)
  }
  # End RM/mixed ANOVA path

  # --- 6. Calculate Effect Size ---
  effect_size_result <- NULL

  if (calculate_effect_size && !is.null(test_result) &&
      requireNamespace("effectsize", quietly = TRUE)) {
    effect_size_result <- tryCatch({
      if (test_family %in% c("independent_two_sample", "paired_two_sample")) {
        is_paired_fam <- (test_family == "paired_two_sample")
        if (parametric_used) {
          es     <- effectsize::cohens_d(
                      x         = y[grp == groups[1]],
                      y         = y[grp == groups[2]],
                      pooled_sd = !variance_violated,
                      paired    = is_paired_fam,
                      verbose   = FALSE
                    )
          interp <- effectsize::interpret(es, rules = "cohen1988")
          list(estimate = es$Cohens_d, ci_low = es$CI_low, ci_high = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric = if (is_paired_fam) "Cohen's d (paired)" else "Cohen's d",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation), es$Cohens_d))
        } else {
          es     <- effectsize::rank_biserial(
                      x       = y[grp == groups[1]],
                      y       = y[grp == groups[2]],
                      paired  = is_paired_fam,
                      verbose = FALSE
                    )
          interp <- effectsize::interpret(es, rules = "funder2019")
          list(estimate = es$r_rank_biserial, ci_low = es$CI_low, ci_high = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric = if (is_paired_fam) "Rank-biserial correlation (paired)"
                        else "Rank-biserial correlation",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation),
                                        es$r_rank_biserial))
        }
      } else if (test_family == "independent_anova") {
        if (parametric_used) {
          es_model <- if (is.null(fitted_model)) {
            warnings_list <<- c(warnings_list,
              "Eta-squared calculated from standard ANOVA (Welch's ANOVA used for test).")
            aov(y ~ grp)
          } else fitted_model
          es     <- effectsize::eta_squared(es_model, verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate = es$Eta2, ci_low = es$CI_low, ci_high = es$CI_high,
               magnitude = as.character(interp$Interpretation), metric = "Eta-squared",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation), es$Eta2))
        } else {
          es     <- effectsize::rank_epsilon_squared(x = y, groups = grp, verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate  = es$rank_epsilon_squared,
               ci_low    = if (!is.null(es$CI_low))  es$CI_low  else NA_real_,
               ci_high   = if (!is.null(es$CI_high)) es$CI_high else NA_real_,
               magnitude = as.character(interp$Interpretation),
               metric    = "Rank epsilon-squared",
               interpretation = sprintf("Effect size: %s (%.3f)",
                                        as.character(interp$Interpretation),
                                        es$rank_epsilon_squared))
        }
      }
    }, error = function(e) {
      msg <- paste("Could not calculate effect size:", e$message)
      warnings_list <<- c(warnings_list, msg)
      if (verbose) message(msg)
      NULL
    })
  }

  # --- 7. Perform Post-Hoc Tests ---
  posthoc_result <- NULL

  if (perform_posthoc && n_groups > 2 &&
    !is.null(test_result) && test_result$p.value < test_alpha &&
    requireNamespace("rstatix", quietly = TRUE)) {
      posthoc_result <- tryCatch({
        test_data <- tibble::tibble(y = y, grp = grp)
        ph        <- NULL

        if (test_used %in% c("One-way ANOVA", "One-way ANOVA (forced)")) {
          ph              <- as.data.frame(rstatix::tukey_hsd(test_data, y ~ grp))
          # Drop rstatix's own hardcoded p.adj.signif immediately — we recompute
          # it below using test_alpha so it respects the user's chosen threshold.
          ph$p.adj.signif <- NULL
          ph$posthoc_test <- "Tukey HSD"
          ph$statistic    <- NA_real_   # tukey_hsd() genuinely has no per-pair statistic
          ph$p            <- NA_real_   # tukey_hsd() genuinely has no raw p
          # Use match() not [] subsetting to safely extract n per group from table()
          n_lookup        <- as.numeric(group_ns)
          names(n_lookup) <- names(group_ns)
          ph$n1           <- n_lookup[match(ph$group1, names(n_lookup))]
          ph$n2           <- n_lookup[match(ph$group2, names(n_lookup))]

        } else if (test_used == "Welch's ANOVA") {
          ph              <- as.data.frame(rstatix::games_howell_test(test_data, y ~ grp))
          # Drop rstatix's own hardcoded p.adj.signif immediately — same reason.
          ph$p.adj.signif <- NULL
          ph$posthoc_test <- "Games-Howell"
          # games_howell_test() returns statistic, n1, n2 natively — do not overwrite.
          # games_howell_test() does NOT return a raw p (only p.adj) — p stays absent.

        } else if (test_used %in% c("Kruskal-Wallis test", "Kruskal-Wallis test (forced)")) {
          ph              <- as.data.frame(rstatix::dunn_test(test_data, y ~ grp,
                                                              p.adjust.method = p_adjust_method))
          # Drop rstatix's own hardcoded p.adj.signif immediately — same reason.
          ph$p.adj.signif <- NULL
          ph$posthoc_test <- "Dunn's test"
          # dunn_test() returns statistic, p, and p.adj natively — do not overwrite.
        }

        # ---- p.format (raw p, only present for Dunn's) ----
        if ("p" %in% names(ph) && !all(is.na(ph$p))) {
          ph$p.format <- ifelse(ph$p < 0.001, "p < 0.001", sprintf("p = %.3f", ph$p))
        } else {
          if (!"p" %in% names(ph)) ph$p <- NA_real_
          ph$p.format <- NA_character_
        }

        # ---- p.adj.signif and significant — both respect test_alpha ----
        if ("p.adj" %in% names(ph)) {
          ph$significant  <- ph$p.adj < test_alpha

          ph$p.adj.signif <- ifelse(
            ph$p.adj < test_alpha,
            ifelse(ph$p.adj < 0.001, "***",
            ifelse(ph$p.adj < 0.01,  "**",
            ifelse(ph$p.adj < 0.05,  "*",
                  sprintf("* (p<%.2f)", test_alpha)))),  # significant but above 0.05
            "ns"
          )

          ph$p.adj.format <- ifelse(
            ph$p.adj < 0.001, "p < 0.001", sprintf("p = %.3f", ph$p.adj)
          )
        } else {
          ph$p.adj.signif <- NA_character_
          ph$p.adj.format <- NA_character_
          ph$significant  <- NA
        }

        # ---- mean lookup via match() to avoid merge() side-effects ----
        summary_stats <- data.frame(
          group = groups,
          mean  = as.numeric(tapply(y, grp, mean)),
          stringsAsFactors = FALSE
        )
        ph$mean1 <- summary_stats$mean[match(ph$group1, summary_stats$group)]
        ph$mean2 <- summary_stats$mean[match(ph$group2, summary_stats$group)]

        # ---- interpretation: driven by ph$significant (respects test_alpha) ----
        if ("significant" %in% names(ph)) {
          if (parametric_used) {
            ph$interpretation <- ifelse(
              !ph$significant %in% TRUE,              # vectorized; handles NA safely
              "No significant difference between groups",
              ifelse(ph$mean1 < ph$mean2,
                    paste0(ph$group1, " has lower mean than ",  ph$group2),
                    paste0(ph$group1, " has higher mean than ", ph$group2))
            )
          } else {
            ranks      <- rank(y)
            mean_ranks <- tapply(ranks, grp, mean)
            ph$interpretation <- mapply(
              function(g1, g2, sig) {
                if (!isTRUE(sig)) return("No significant difference between groups")
                if (mean_ranks[g1] < mean_ranks[g2])
                  paste0(g1, " tends to have smaller values than ", g2)
                else
                  paste0(g1, " tends to have larger values than ", g2)
              },
              ph$group1, ph$group2, ph$significant,
              SIMPLIFY = TRUE, USE.NAMES = FALSE
            )
          }
        } else {
          ph$interpretation <- NA_character_
        }

        # ---- assemble final columns ----
        final_cols <- c("group1", "group2", "n1", "n2", "statistic", "p", "p.adj",
                        "posthoc_test", "p.format", "p.adj.format", "p.adj.signif",
                        "significant", "mean1", "mean2", "interpretation")
        for (col in final_cols) if (!col %in% names(ph)) ph[[col]] <- NA
        ph[, final_cols]

      }, error = function(e) {
        msg <- paste("Could not perform post-hoc test:", e$message)
        warnings_list <<- c(warnings_list, msg)
        if (verbose) message(msg)
        NULL
      })
  }

  # --- 8. Comprehensive Data Summary ---
  present_groups <- levels(droplevels(grp))
  present_ns     <- table(droplevels(grp))
  data_summary <- data.frame(
    group  = present_groups,
    n      = as.numeric(present_ns),
    mean   = as.numeric(tapply(y, droplevels(grp), mean)),
    sd     = as.numeric(tapply(y, droplevels(grp), sd)),
    se     = as.numeric(tapply(y, droplevels(grp), function(v) sd(v) / sqrt(length(v)))),
    median = as.numeric(tapply(y, droplevels(grp), median)),
    q25    = as.numeric(tapply(y, droplevels(grp), function(v) quantile(v, 0.25))),
    q75    = as.numeric(tapply(y, droplevels(grp), function(v) quantile(v, 0.75))),
    min    = as.numeric(tapply(y, droplevels(grp), min)),
    max    = as.numeric(tapply(y, droplevels(grp), max)),
    row.names = NULL
  )
  if (requireNamespace("moments", quietly = TRUE)) {
    data_summary$skewness <- as.numeric(tapply(y, droplevels(grp), moments::skewness))
    data_summary$kurtosis <- as.numeric(tapply(y, droplevels(grp), moments::kurtosis))
  }

  # --- 9. Print Results ---
  if (verbose) {
    message("\n--- Test Results ---")
    message(sprintf("Test Selected: %s", test_used))
    message(sprintf("Test statistic: %.4f", test_result$statistic))
    if (test_result$p.value < 0.001) message("P-value: < 0.001") else
      message(sprintf("P-value: %.4f", test_result$p.value))
    if (!is.null(effect_size_result)) message(sprintf("%s", effect_size_result$interpretation))
    if (test_result$p.value < test_alpha)
      message(sprintf("Interpretation: Significant difference detected (p < %.2f) ***", test_alpha))
    else
      message(sprintf("Interpretation: No significant difference detected (p >= %.2f)", test_alpha))
    if (!is.null(posthoc_result)) {
      message(sprintf("\n--- Post-Hoc Comparisons (%s) ---", unique(posthoc_result$posthoc_test)))
      sig_comparisons   <- sum(posthoc_result$significant, na.rm = TRUE)
      total_comparisons <- nrow(posthoc_result)
      message(sprintf("Significant pairwise differences: %d out of %d",
                      sig_comparisons, total_comparisons))
    }
  }

  # --- 9b. Subgroup Analysis ---
  subgroup_analysis <- NULL

  if (!is.null(subgroup)) {
    subgroup_analysis <- list()
    for (sg_var in subgroup) {
      sg_levels <- if (is.factor(x[[sg_var]])) {
        levels(droplevels(as.factor(x[[sg_var]])))
      } else {
        sort(unique(na.omit(as.character(x[[sg_var]]))))
      }
      sg_var_results <- list()
      for (lvl in sg_levels) {
        x_sub   <- x[as.character(x[[sg_var]]) == lvl, , drop = FALSE]
        if (nrow(x_sub) == 0) {
          warning(paste0("Subgroup '", sg_var, "' level '", lvl, "' has no observations. Skipping."))
          next
        }
        grp_sub <- droplevels(as.factor(x_sub[[group]]))
        if (nlevels(grp_sub) < 2) {
          warning(paste0("Subgroup '", sg_var, "' = '", lvl,
                         "': fewer than 2 group levels after subsetting. Skipping."))
          next
        }
        sg_result <- tryCatch(
          run_diff(
            x = x_sub, outcome = outcome, group = group, 
            within = within, subject_id = subject_id,
            paired = paired,
            test_type = test_type, normality_method = normality_method,
            alpha_normality = alpha_normality, alpha_variance = alpha_variance,
            alpha_sphericity = alpha_sphericity,
            test_alpha = test_alpha, min_n_threshold = min_n_threshold,
            calculate_effect_size = calculate_effect_size,
            perform_posthoc = perform_posthoc, p_adjust_method = p_adjust_method,
            group_order = group_order, verbose = verbose, num_cores = num_cores,
            summary_table = FALSE, subgroup = NULL
          ),
          error = function(e) {
            warning(paste0("run_diff failed for subgroup '", sg_var, "' = '", lvl,
                           "': ", e$message), call. = FALSE)
            NULL
          }
        )
        if (!is.null(sg_result)) sg_var_results[[lvl]] <- sg_result
      }
      subgroup_analysis[[sg_var]] <- sg_var_results
    }
  }

  # --- 10. Return Object ---
  result <- list(
    test_used         = test_used,
    test_result       = test_result,
    effect_size       = effect_size_result,
    posthoc_result    = posthoc_result,
    assumptions       = assumptions,
    parametric        = parametric_used,
    data_summary      = data_summary,
    outcome           = outcome,
    group_var         = group,
    paired            = paired,
    test_type         = test_type,
    test_alpha        = test_alpha,
    subgroup_analysis = subgroup_analysis,
    warnings          = warnings_list,
    parameters        = list(
      test_type             = test_type,
      normality_method      = normality_method,
      alpha_normality       = alpha_normality,
      alpha_variance        = alpha_variance,
      alpha_sphericity      = alpha_sphericity,
      type_sosquares        = type_sosquares,
      rm_effect_size        = rm_effect_size,
      correction            = correction,
      test_alpha            = test_alpha,
      min_n_threshold       = min_n_threshold,
      calculate_effect_size = calculate_effect_size,
      perform_posthoc       = perform_posthoc,
      p_adjust_method       = p_adjust_method,
      group_order           = group_order,
      subgroup              = subgroup,
      within                = within,
      subject_id            = subject_id,
      paired                = paired
    ),
    raw_data = tibble::tibble(outcome = y, group = grp)
  )

  return(result)
  }


# =============================================================================
# S3 methods  (unchanged)
# =============================================================================

#' Print method for run_diff objects
#' @param x An object from \code{run_diff}.
#' @param ... Further arguments passed to or from other methods.
#' @export
print.run_diff <- function(x, ...) {
  cat("\n=== Automatic Statistical Comparison ===\n\n")
  cat(sprintf("Test Used: %s\n", x$test_used))
  cat(sprintf("Parametric: %s\n\n", x$parametric))
  cat("Test Results:\n")
  cat(sprintf("  Statistic: %.4f\n", x$test_result$statistic[[1]]))
  if (!is.null(x$test_result$p.value)) {
    if (x$test_result$p.value < 0.001) {
      cat("  P-value: < 0.001 ***\n")
    } else {
      cat(sprintf("  P-value: %.4f", x$test_result$p.value))
      if      (x$test_result$p.value < 0.01) cat(" **")
      else if (x$test_result$p.value < 0.05) cat(" *")
      cat("\n")
    }
  }
  if (!is.null(x$effect_size))
    cat(sprintf("  %s\n", x$effect_size$interpretation))
  cat("\n")
  # ---- full ANOVA table for RM/mixed designs ----
  if (!is.null(x$anova_table)) {
    cat("ANOVA Table (all effects):\n")
    print(x$anova_table, row.names = FALSE)
  } else if (!is.null(x$data_summary)) {
    cat("Data Summary:\n")
    # Safely select only columns that actually exist
    show_cols <- c("group", "n", "mean", "sd", "median")
    show_cols <- show_cols[show_cols %in% names(x$data_summary)]
    if (length(show_cols) >= 2) {
      print(x$data_summary[, show_cols, drop = FALSE], row.names = FALSE)
    } else {
      print(x$data_summary, row.names = FALSE)
    }
  }
  if (!is.null(x$posthoc_result)) {
    cat(sprintf("\n\nPost-Hoc Tests (%s):\n",
                paste(unique(x$posthoc_result$posthoc_test), collapse = ", ")))
    sig_only <- x$posthoc_result[x$posthoc_result$significant %in% TRUE, , drop = FALSE]
    # Select columns that exist
    ph_show <- c("group1", "group2", "effect", "p.adj", "p.adj.signif")
    ph_show <- ph_show[ph_show %in% names(sig_only)]
    if (nrow(sig_only) > 0) {
      cat("Significant comparisons:\n")
      print(sig_only[, ph_show, drop = FALSE], row.names = FALSE)
    } else {
      cat("No significant pairwise differences found.\n")
    }
  }
  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) cat(sprintf("  - %s\n", w))
  }
}


#' Summary method for run_diff objects
#' @param object An object from \code{run_diff}.
#' @param ... Further arguments passed to or from other methods.
#' @export
summary.run_diff <- function(object, ...) {
  cat("\n=== Detailed Statistical Comparison Summary ===\n\n")
  cat(sprintf("Outcome Variable:  %s\n", object$outcome))
  cat(sprintf("Grouping Variable: %s\n", object$group_var))
  cat(sprintf("Paired Analysis:   %s\n", object$paired))
  cat(sprintf("Test Strategy:     %s (alpha = %.2f)\n", object$test_type, object$test_alpha))
  cat(sprintf("Test Selected:     %s\n\n", object$test_used))
  if (!is.null(object$within_vars))
    cat(sprintf("Within-Subject Factor(s): %s\n",
                paste(object$within_vars, collapse = ", ")))
  if (!is.null(object$subject_id_var))
    cat(sprintf("Subject ID Variable:      %s\n", object$subject_id_var))
  cat("--- Descriptive Statistics ---\n")
  print(object$data_summary, row.names = FALSE)
  cat("\n--- Assumption Checks ---\n")
  print(object$assumptions, row.names = FALSE)
  # Full ANOVA table for RM/mixed designs
  if (!is.null(object$anova_table)) {
    cat("\n--- Full ANOVA Table (all effects) ---\n")
    print(object$anova_table, row.names = FALSE)
  }
  if (!is.null(object$sphericity)) {
    cat("\n--- Mauchly's Test for Sphericity ---\n")
    print(as.data.frame(object$sphericity), row.names = FALSE)
  }
  if (!is.null(object$sphericity_corrections)) {
    cat("\n--- Sphericity Corrections (GG / HF) ---\n")
    print(as.data.frame(object$sphericity_corrections), row.names = FALSE)
  }
  cat("\n--- Statistical Test Results ---\n")
  if (is.null(object$anova_table)) {
    # Classical design: print the htest object normally
    print(object$test_result)
  } else {
    # RM/mixed: anova_table already printed above; just show primary effect summary
    tr <- object$test_result
    cat(sprintf("Primary effect: F(%s, %s) = %.4f, p = %s\n",
                tr$parameter["DFn"], tr$parameter["DFd"],
                tr$statistic[[1]],
                if (tr$p.value < 0.001) "< 0.001" else sprintf("%.4f", tr$p.value)))
  }
  if (!is.null(object$effect_size)) {
    cat("\n--- Effect Size ---\n")
    cat(sprintf("Metric:    %s\n", object$effect_size$metric))
    cat(sprintf("Estimate:  %.3f", object$effect_size$estimate))
    if (!is.na(object$effect_size$ci_low))
      cat(sprintf(" [95%% CI: %.3f, %.3f]", object$effect_size$ci_low, object$effect_size$ci_high))
    cat("\n")
    cat(sprintf("Magnitude: %s\n", object$effect_size$magnitude))
  }
  if (!is.null(object$posthoc_result)) {
    cat(sprintf("\n--- Post-Hoc Pairwise Comparisons (%s) ---\n",
                unique(object$posthoc_result$posthoc_test)))
    cat(sprintf("P-value adjustment method: %s\n\n", object$parameters$p_adjust_method))
    cols_to_show <- c("group1", "group2", "n1", "n2", "statistic", "p",
                      "p.adj", "p.adj.signif", "interpretation")
    cols_to_show <- cols_to_show[cols_to_show %in% names(object$posthoc_result)]
    print(object$posthoc_result[, cols_to_show], row.names = FALSE)
  }
  cat("\n--- Interpretation ---\n")
  if (object$test_result$p.value < object$test_alpha) {
    cat(sprintf("Result: SIGNIFICANT difference detected between groups (p < %.2f).\n",
                object$test_alpha))
    if (!is.null(object$effect_size))
      cat(sprintf("The effect size is %s (%s = %.3f).\n",
                  tolower(object$effect_size$magnitude),
                  object$effect_size$metric,
                  object$effect_size$estimate))
  } else {
    cat(sprintf("Result: NO significant difference detected between groups (p >= %.2f).\n",
                object$test_alpha))
  }
  if (length(object$warnings) > 0) {
    cat("\n--- Warnings ---\n")
    for (w in object$warnings) cat(sprintf("  - %s\n", w))
  }
  invisible(object)
}


# =============================================================================
#  .run_diff_single for plot_diff function
# =============================================================================

.run_diff_single <- function(
    x,
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
    verbose               = FALSE,
    num_cores             = 1L
) {
  test_type        <- match.arg(test_type)
  normality_method <- match.arg(normality_method)
  warnings_list    <- character(0)

  # ---- num_cores ----
  avail_cores <- parallel::detectCores()
  if (is.na(avail_cores)) avail_cores <- 1L
  if (identical(num_cores, "max")) {
    num_cores <- max(1L, avail_cores - 2L)
  } else {
    num_cores <- as.integer(num_cores)
  }

  # ---- basic checks ----
  if (!is.data.frame(x))
    stop("'x' must be a data frame.")
  if (!outcome %in% names(x))
    stop(sprintf("Outcome '%s' not found in data.", outcome))
  if (!group %in% names(x))
    stop(sprintf("Group '%s' not found in data.", group))

  y   <- x[[outcome]]
  grp <- as.factor(x[[group]])
  if (!is.numeric(y))
    stop("Outcome variable must be numeric.")

  if (!is.null(group_order)) {
    if (!all(group_order %in% levels(grp)))
      stop("'group_order' contains levels not present in the data.")
    grp <- factor(grp, levels = group_order)
  }

  cc <- stats::complete.cases(y, grp)
  if (sum(!cc) > 0) {
    msg           <- sprintf("Removing %d rows with missing values.", sum(!cc))
    warnings_list <- c(warnings_list, msg)
  }
  y   <- y[cc]
  grp <- droplevels(grp[cc])

  if (nlevels(grp) < 2)
    stop("Grouping factor must have at least 2 levels after removing missing values.")

  groups   <- levels(grp)
  n_groups <- length(groups)
  group_ns <- table(grp)
  min_n    <- min(group_ns)

  if (min_n < min_n_threshold)
    warnings_list <- c(warnings_list,
                       sprintf("Min sample size (%d) below threshold (%d).",
                               min_n, min_n_threshold))

  # ---- test family ----
  if (n_groups == 2) {
    test_family <- if (paired) "paired_two_sample" else "independent_two_sample"
  } else {
    test_family <- if (paired) "repeated_measures_anova" else "independent_anova"
  }
  if (test_family == "repeated_measures_anova")
    stop("Repeated measures ANOVA is not yet implemented.")
  if (paired && n_groups == 2 && group_ns[1] != group_ns[2])
    stop("Paired tests require equal group sizes after removing missing values.")

  # ---- normality helper ----
  # .check_norm <- function(data_vector, group_label = NULL) {
  #   n <- length(data_vector)
  #   if (n < 3)
  #     return(list(violated = TRUE, p_value = NA_real_,
  #                 method = "Insufficient data", label = group_label))
  #   use_m <- if (normality_method == "auto") {
  #     if (n < 50) "shapiro" else "lilliefors"
  #   } else normality_method

  #   if (use_m == "shapiro") {
  #     res <- stats::shapiro.test(data_vector)
  #     return(list(violated = res$p.value < alpha_normality,
  #                 p_value  = res$p.value,
  #                 method   = "Shapiro-Wilk",
  #                 label    = group_label))
  #   } else {
  #     if (!requireNamespace("nortest", quietly = TRUE)) {
  #       res <- stats::shapiro.test(data_vector)
  #       warnings_list <<- c(warnings_list,
  #                           "nortest unavailable; fell back to Shapiro-Wilk.")
  #       return(list(violated = res$p.value < alpha_normality,
  #                   p_value  = res$p.value,
  #                   method   = "Shapiro-Wilk (fallback)",
  #                   label    = group_label))
  #     }
  #     res <- nortest::lillie.test(data_vector)
  #     return(list(violated = res$p.value < alpha_normality,
  #                 p_value  = res$p.value,
  #                 method   = "Lilliefors (K-S)",
  #                 label    = group_label))
  #   }
  # }
  .check_norm <- function(data_vector, group_label = NULL) {
  n <- length(data_vector)
  if (n < 3)
    return(list(violated = TRUE, p_value = NA_real_,
                method = "Insufficient data", label = group_label))

  if (length(unique(data_vector)) < 2) {
    msg <- if (!is.null(group_label))
      sprintf("All values in group '%s' are identical; normality assumption violated.", group_label)
    else
      "All values are identical; normality assumption violated."
    warnings_list <<- c(warnings_list, msg)
    return(list(violated = TRUE, p_value = NA_real_,
                method = "Skipped (zero variance)", label = group_label))
  }

  use_m <- if (normality_method == "auto") {
    if (n < 50) "shapiro" else "lilliefors"
  } else normality_method

  if (use_m == "shapiro") {
    res <- stats::shapiro.test(data_vector)
    return(list(violated = res$p.value < alpha_normality,
                p_value  = res$p.value,
                method   = "Shapiro-Wilk",
                label    = group_label))
  } else {
    if (!requireNamespace("nortest", quietly = TRUE)) {
      res <- stats::shapiro.test(data_vector)
      warnings_list <<- c(warnings_list,
                          "nortest unavailable; fell back to Shapiro-Wilk.")
      return(list(violated = res$p.value < alpha_normality,
                  p_value  = res$p.value,
                  method   = "Shapiro-Wilk (fallback)",
                  label    = group_label))
    }
    res <- nortest::lillie.test(data_vector)
    return(list(violated = res$p.value < alpha_normality,
                p_value  = res$p.value,
                method   = "Lilliefors (K-S)",
                label    = group_label))
  }
}

  # ---- assumption checks ----
  normality_violated <- FALSE
  variance_violated  <- FALSE
  norm_results_raw   <- list()
  variance_p         <- NA_real_
  variance_method    <- NA_character_

  if (test_type == "auto") {
    if (n_groups == 2) {
      run_nc <- function(g) .check_norm(y[grp == g], g)
      norm_results_raw <- if (num_cores > 1 && .Platform$OS.type != "windows") {
        parallel::mclapply(groups, run_nc, mc.cores = num_cores)
      } else if (num_cores > 1) {
        cl <- parallel::makeCluster(num_cores)
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        parallel::clusterExport(
          cl,
          varlist = c("y", "grp", ".check_norm", "normality_method", "alpha_normality"),
          envir   = environment()
        )
        parallel::clusterEvalQ(cl, {
          if (requireNamespace("nortest", quietly = TRUE)) library(nortest)
        })
        parallel::parLapply(cl, groups, run_nc)
      } else {
        lapply(groups, run_nc)
      }
    } else {
      fit              <- stats::lm(y ~ grp)
      nr               <- .check_norm(fit$residuals, "residuals")
      norm_results_raw <- list(nr)
    }

    for (nr in norm_results_raw) {
      if (inherits(nr, "try-error")) next
      if (isTRUE(nr$violated)) normality_violated <- TRUE
    }

    if (!paired) {
      variance_p <- tryCatch({
        if (requireNamespace("car", quietly = TRUE)) {
          variance_method <- "Levene's test"
          car::leveneTest(y, grp)$`Pr(>F)`[1]
        } else {
          variance_method <- "Bartlett's test"
          stats::bartlett.test(x = y, g = grp)$p.value
        }
      }, error = function(e) {
        warnings_list <<- c(warnings_list,
                            "Variance test failed; assumed homogeneous.")
        1
      })
      if (!is.na(variance_p) && variance_p < alpha_variance)
        variance_violated <- TRUE
    }
  }

  # ---- select & run test ----
  test_used       <- NULL
  test_result     <- NULL
  parametric_used <- TRUE
  fitted_model    <- NULL

  if (test_type == "auto") {
    if (test_family == "independent_two_sample") {
      if (normality_violated) {
        test_used       <- "Mann-Whitney U test (Wilcoxon rank-sum)"
        test_result     <- stats::wilcox.test(x = y[grp == groups[1]],
                                              y = y[grp == groups[2]])
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's t-test"
        test_result <- stats::t.test(x = y[grp == groups[1]],
                                     y = y[grp == groups[2]], var.equal = FALSE)
      } else {
        test_used   <- "Student's t-test"
        test_result <- stats::t.test(x = y[grp == groups[1]],
                                     y = y[grp == groups[2]], var.equal = TRUE)
      }
    } else if (test_family == "paired_two_sample") {
      g1 <- y[grp == groups[1]]
      g2 <- y[grp == groups[2]]
      if (normality_violated) {
        test_used       <- "Wilcoxon signed-rank test"
        test_result     <- stats::wilcox.test(g1, g2, paired = TRUE)
        parametric_used <- FALSE
      } else {
        test_used   <- "Paired t-test"
        test_result <- stats::t.test(g1, g2, paired = TRUE)
      }
    } else {
      # independent_anova
      if (normality_violated) {
        test_used       <- "Kruskal-Wallis test"
        test_result     <- stats::kruskal.test(x = y, g = grp)
        parametric_used <- FALSE
      } else if (variance_violated) {
        test_used   <- "Welch's ANOVA"
        test_result <- stats::oneway.test(stats::formula(y ~ grp), var.equal = FALSE)
      } else {
        test_used    <- "One-way ANOVA"
        fitted_model <- stats::aov(y ~ grp)
        sa           <- summary(fitted_model)[[1]]
        test_result  <- list(
          statistic = c("F value" = sa$`F value`[1]),
          parameter = c(df1 = sa$Df[1], df2 = sa$Df[2]),
          p.value   = sa$`Pr(>F)`[1],
          method    = test_used
        )
      }
    }
  } else if (test_type == "parametric") {
    if (test_family == "independent_two_sample") {
      test_used   <- "Student's t-test (forced)"
      test_result <- stats::t.test(x = y[grp == groups[1]],
                                   y = y[grp == groups[2]], var.equal = TRUE)
    } else if (test_family == "paired_two_sample") {
      test_used   <- "Paired t-test (forced)"
      test_result <- stats::t.test(y[grp == groups[1]], y[grp == groups[2]],
                                   paired = TRUE)
    } else {
      test_used    <- "One-way ANOVA (forced)"
      fitted_model <- stats::aov(y ~ grp)
      sa           <- summary(fitted_model)[[1]]
      test_result  <- list(
        statistic = c("F value" = sa$`F value`[1]),
        parameter = c(df1 = sa$Df[1], df2 = sa$Df[2]),
        p.value   = sa$`Pr(>F)`[1],
        method    = test_used
      )
    }
  } else {
    # nonparametric
    parametric_used <- FALSE
    if (test_family == "independent_two_sample") {
      test_used   <- "Mann-Whitney U test (forced)"
      test_result <- stats::wilcox.test(x = y[grp == groups[1]],
                                        y = y[grp == groups[2]])
    } else if (test_family == "paired_two_sample") {
      test_used   <- "Wilcoxon signed-rank test (forced)"
      test_result <- withCallingHandlers(
        stats::wilcox.test(x = y[grp == groups[1]], y = y[grp == groups[2]], paired = TRUE),
        warning = function(w) {
          if (grepl("ties", conditionMessage(w), ignore.case = TRUE))
            warnings_list <<- c(warnings_list, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    } else {
      test_used   <- "Kruskal-Wallis test (forced)"
      test_result <- stats::kruskal.test(x = y, g = grp)
    }
  }

  # ---- effect size ----
  es_result <- NULL
  if (calculate_effect_size && !is.null(test_result) &&
      requireNamespace("effectsize", quietly = TRUE)) {
    es_result <- tryCatch({
      if (test_family %in% c("independent_two_sample", "paired_two_sample")) {
        is_paired_fam <- (test_family == "paired_two_sample")
        if (parametric_used) {
          es     <- effectsize::cohens_d(
                      x         = y[grp == groups[1]],
                      y         = y[grp == groups[2]],
                      #pooled_sd = !variance_violated,
                      paired    = TRUE,
                      verbose   = FALSE
                    )
          interp <- effectsize::interpret(es, rules = "cohen1988")
          list(estimate  = es$Cohens_d,
               ci_low    = es$CI_low,
               ci_high   = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric    = if (is_paired_fam) "Cohen's d (paired)" else "Cohen's d")
        } else {
          es     <- effectsize::rank_biserial(
                      x       = y[grp == groups[1]],
                      y       = y[grp == groups[2]],
                      paired  = TRUE,
                      verbose = FALSE
                    )
          interp <- effectsize::interpret(es, rules = "funder2019")
          list(estimate  = es$r_rank_biserial,
               ci_low    = es$CI_low,
               ci_high   = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric    = if (is_paired_fam) "Rank-biserial correlation (paired)"
                           else "Rank-biserial correlation")
        }
      } else {
        # independent_anova
        if (parametric_used) {
          es_model <- if (is.null(fitted_model)) {
            warnings_list <<- c(warnings_list,
              "Eta-squared calculated from standard ANOVA (Welch's ANOVA used for test).")
            stats::aov(y ~ grp)
          } else fitted_model
          es     <- effectsize::eta_squared(es_model, verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate  = es$Eta2,
               ci_low    = es$CI_low,
               ci_high   = es$CI_high,
               magnitude = as.character(interp$Interpretation),
               metric    = "Eta-squared")
        } else {
          es     <- effectsize::rank_epsilon_squared(x = y, groups = grp,
                                                     verbose = FALSE)
          interp <- effectsize::interpret(es, rules = "field2013")
          list(estimate  = es$rank_epsilon_squared,
               ci_low    = if (!is.null(es$CI_low))  es$CI_low  else NA_real_,
               ci_high   = if (!is.null(es$CI_high)) es$CI_high else NA_real_,
               magnitude = as.character(interp$Interpretation),
               metric    = "Rank epsilon-squared")
        }
      }
    }, error = function(e) {
      warnings_list <<- c(warnings_list, paste("Effect size failed:", e$message))
      NULL
    })
  }

  # ---- post-hoc ----
  posthoc_result <- NULL
  if (perform_posthoc && n_groups > 2 &&
      !is.null(test_result) && test_result$p.value < test_alpha &&
      requireNamespace("rstatix", quietly = TRUE)) {
    posthoc_result <- tryCatch({
      td <- tibble::tibble(y = y, grp = grp)
      ph <- NULL

      if (test_used %in% c("One-way ANOVA", "One-way ANOVA (forced)")) {
        ph              <- rstatix::tukey_hsd(td, y ~ grp)
        ph$posthoc_test <- "Tukey HSD"
        ph$statistic    <- NA_real_
        ph$p            <- NA_real_
        nl              <- stats::setNames(as.numeric(group_ns), names(group_ns))
        ph$n1           <- nl[ph$group1]
        ph$n2           <- nl[ph$group2]
      } else if (test_used == "Welch's ANOVA") {
        ph              <- rstatix::games_howell_test(td, y ~ grp)
        ph$posthoc_test <- "Games-Howell"
        ph$p            <- NA_real_
      } else if (test_used %in% c("Kruskal-Wallis test",
                                   "Kruskal-Wallis test (forced)")) {
        ph              <- rstatix::dunn_test(td, y ~ grp,
                                              p.adjust.method = p_adjust_method)
        ph$posthoc_test <- "Dunn's test"
      }

      if (!is.null(ph)) {
        if ("p.adj" %in% names(ph)) {
          ph$p.adj.signif <- ifelse(ph$p.adj < 0.001, "***",
                             ifelse(ph$p.adj < 0.01,  "**",
                             ifelse(ph$p.adj < 0.05,  "*", "ns")))
          ph$significant  <- ph$p.adj < test_alpha
        }

        # Stochastic dominance interpretation via mean ranks
        rks        <- rank(y)
        mean_ranks <- tapply(rks, grp, mean)
        ph$interpretation <- mapply(function(g1, g2, sig) {
          if (is.na(sig) || sig == "ns") return("No significant difference")
          if (mean_ranks[g1] < mean_ranks[g2])
            paste0(g1, " tends to have smaller values than ", g2)
          else
            paste0(g1, " tends to have larger values than ", g2)
        },
        ph$group1, ph$group2,
        if ("p.adj.signif" %in% names(ph)) ph$p.adj.signif else rep("ns", nrow(ph)),
        SIMPLIFY = TRUE, USE.NAMES = FALSE)
      }
      ph
    }, error = function(e) {
      warnings_list <<- c(warnings_list, paste("Post-hoc failed:", e$message))
      NULL
    })
  }

  # ---- data summary ----
  dg <- droplevels(grp)
  data_summary <- data.frame(
    group     = levels(dg),
    n         = as.numeric(table(dg)),
    mean      = as.numeric(tapply(y, dg, mean)),
    sd        = as.numeric(tapply(y, dg, stats::sd)),
    se        = as.numeric(tapply(y, dg,
                  function(v) stats::sd(v) / sqrt(length(v)))),
    median    = as.numeric(tapply(y, dg, stats::median)),
    q25       = as.numeric(tapply(y, dg,
                  function(v) stats::quantile(v, 0.25))),
    q75       = as.numeric(tapply(y, dg,
                  function(v) stats::quantile(v, 0.75))),
    min       = as.numeric(tapply(y, dg, min)),
    max       = as.numeric(tapply(y, dg, max)),
    row.names = NULL
  )

  # ---- return rich list (mirrors old run_diff structure for plot_diff) ----
  list(
    test_used          = test_used,
    test_result        = test_result,
    effect_size        = es_result,
    posthoc_result     = posthoc_result,
    norm_results       = norm_results_raw,
    variance_p         = variance_p,
    variance_method    = variance_method,
    normality_violated = normality_violated,
    variance_violated  = variance_violated,
    parametric         = parametric_used,
    data_summary       = data_summary,
    outcome            = outcome,
    group_var          = group,
    groups             = groups,
    n_groups           = n_groups,
    group_ns           = group_ns,
    paired             = paired,
    test_type          = test_type,
    test_alpha         = test_alpha,
    warnings           = warnings_list,
    raw_data           = tibble::tibble(outcome = y, group = grp)
  )
}



# =============================================================================
#  Helper: build a summary_table data frame from a named list of run_diff
#          objects (one element per outcome).  Used both for the main table
#          and for each subgroup × level slice.
# =============================================================================
.build_summary_table <- function(results_list, outcome_names, test_alpha) {

  # ------------------------------------------------------------------
  # Pass 1: collect all group names that appear across any outcome so
  # we can build a consistent, union set of average_* columns.
  # ------------------------------------------------------------------
  all_groups <- character(0)
  for (oc in outcome_names) {
    res <- results_list[[oc]]
    if (!is.null(res) && !is.null(res$data_summary)) {
      all_groups <- union(all_groups, as.character(res$data_summary$group))
    }
  }
  avg_col_names <- if (length(all_groups) > 0)
    paste0("average_", all_groups) else character(0)

  # ------------------------------------------------------------------
  # Pass 2: build one row per outcome
  # ------------------------------------------------------------------
  tbl_rows <- lapply(outcome_names, function(oc) {
    res <- results_list[[oc]]

    # ---- skeleton row for failed / NULL outcomes ----
    if (is.null(res)) {
      base <- data.frame(
        outcome               = oc,
        test_used             = NA_character_,
        statistic             = NA_real_,
        df                    = NA_character_,
        p_value               = NA_real_,
        effect_size_metric    = NA_character_,
        effect_size_estimate  = NA_real_,
        effect_size_ci_low    = NA_real_,
        effect_size_ci_high   = NA_real_,
        effect_size_magnitude = NA_character_,
        significant           = NA,
        n_significant_posthoc = NA_integer_,
        posthoc_pairs         = NA_integer_,
        stringsAsFactors      = FALSE
      )
      # fill every average_* column with NA
      for (cn in avg_col_names) base[[cn]] <- NA_real_
      return(base)
    }

    tr <- res$test_result
    es <- res$effect_size
    ph <- res$posthoc_result
    ds <- res$data_summary   # data frame: columns group, n, mean, …

    df_val <- if (!is.null(tr$parameter) && length(tr$parameter) > 0)
      paste(names(tr$parameter), round(tr$parameter, 4), sep = " = ", collapse = ", ")
    else NA_character_

    n_sig_ph <- if (!is.null(ph) && nrow(ph) > 0)
      sum(ph$significant, na.rm = TRUE) else NA_integer_
    total_ph <- if (!is.null(ph) && nrow(ph) > 0) nrow(ph) else NA_integer_

    # ---- per-group averages ----
    avg_vals <- stats::setNames(
      vector("list", length(avg_col_names)),
      avg_col_names
    )
    for (cn in avg_col_names) {
      grp_label <- sub("^average_", "", cn)
      if (!is.null(ds) && grp_label %in% ds$group) {
        avg_vals[[cn]] <- ds$mean[ds$group == grp_label]
      } else {
        avg_vals[[cn]] <- NA_real_
      }
    }

    # ---- assemble ----
    row <- data.frame(
      outcome               = oc,
      test_used             = res$test_used,
      statistic             = if (!is.null(tr$statistic))
                                as.numeric(tr$statistic[[1]]) else NA_real_,
      df                    = df_val,
      p_value               = if (!is.null(tr$p.value)) tr$p.value else NA_real_,
      stringsAsFactors      = FALSE
    )

    # insert average_* columns right after p_value
    for (cn in avg_col_names) row[[cn]] <- as.numeric(avg_vals[[cn]])

    # append effect-size + significance columns
    row$effect_size_metric    <- if (!is.null(es)) es$metric    else NA_character_
    row$effect_size_estimate  <- if (!is.null(es)) es$estimate  else NA_real_
    row$effect_size_ci_low    <- if (!is.null(es)) es$ci_low    else NA_real_
    row$effect_size_ci_high   <- if (!is.null(es)) es$ci_high   else NA_real_
    row$effect_size_magnitude <- if (!is.null(es)) es$magnitude else NA_character_
    row$significant           <- if (!is.null(tr$p.value))
                                   tr$p.value < test_alpha else NA
    row$n_significant_posthoc <- n_sig_ph
    row$posthoc_pairs         <- total_ph

    row
  })

  tbl <- do.call(rbind, tbl_rows)
  rownames(tbl) <- NULL
  tbl
}