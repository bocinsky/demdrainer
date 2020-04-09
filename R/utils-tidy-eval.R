#' Tidy eval helpers
#'
#' These functions provide tidy eval-compatible ways to capture
#' symbols (`sym()`, `syms()`, `ensym()`), expressions (`expr()`,
#' `exprs()`, `enexpr()`), and quosures (`quo()`, `quos()`, `enquo()`).
#' To learn more about tidy eval and how to use these tools, read
#' <http://rlang.tidyverse.org/articles/tidy-evaluation.html>
#'
#' @name tidyeval
#' @keywords internal
#' @aliases          quo quos enquo sym syms ensym expr exprs enexpr quo_name
#' @importFrom rlang quo quos enquo sym syms ensym expr exprs enexpr quo_name
#' @export           quo quos enquo sym syms ensym expr exprs enexpr quo_name
#' @importFrom rlang UQ UQS .data :=
NULL

# Make CRAN check not complain about "." and package data
if (getRversion() >= "2.15.1") utils::globalVariables(c("mukey",
                                                        "muname",
                                                        "FCode",
                                                        "AreaSqKm",
                                                        "GNIS_Name"))

#' Polygon of bounding box around McPhee Reservoir, Colorado
#'
#' A dataset containing the bounding box around McPhee Reservoir, Colorado.
#'
#' @format A geometry set with one feature.
"mcphee"
