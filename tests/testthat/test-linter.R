#if (requireNamespace("lintr", quietly = TRUE)) {
#  context("lints")
#  test_that("Package Style", {
#    lintr::expect_lint_free(linters = with_defaults(object_name_linter = NULL,
#                                                    line_length_linter = NULL,
#                                                    spaces_left_parentheses_linter = NULL))
#  })
#}
