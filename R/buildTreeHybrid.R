#' @title Generate decision tree from single-cell clustering output.
#' @description Uses decision tree procudure to generate a set of rules for each
#'  cell cluster defined by a single-cell clustering.  Splits are determined by
#'  one of two metrics at each split: a one-off metric to determine rules for
#'  identifying clusters by a single feature, and a balanced metric to determine
#'  rules for identifying sets of similar clusters.
#' @param features A L(features) by N(samples) numeric matrix.
#' @param class A vector of K label assignemnts.
#' @param oneoffMetric A character string. What one-off metric to run, either
#'  `modified F1` or  `pairwise AUC`.
#' @param threshold A numeric value. The threshold for the oneoff metric to use
#'  between 0 and 1, 0.95 by default. Smaller values will result is more one-off
#'  splits.
#' @param reuseFeatures Logical.  Whether or not a feature can be used more than
#'  once on the same cluster. Default is TRUE.
#' @param altSplit Logical. Whether or not to force a marker for clusters that
#'  are solely defined by the absence of markers. Defulault is TRUE
#' @param consecutiveOneoff Logical. Whether or not to allow one-off splits at
#'  consecutive brances. Default it TRUE
#' @return A named list with five elements.
#' \itemize{
#'   \item rules - A named list with one `data.frame` for every label. Each
#' `data.frame` has five columns and gives the set of rules for disinguishing
#'  each label.
#'   \itemize{
#'    \item feature - Feature identifier.
#'    \item direction - Relationship to feature value, -1 if less than, 1 if
#'  greater than.
#'    \item value - The feature value which defines the decision boundary
#'    \item stat - The performance value returned by the splitting metric for
#'  this split.
#'    \item statUsed - Which performance metric was used. "IG" if information
#'  gain and "OO" if one-off.
#'    \item level - The level of the tree at which is rule was defined. 1 is the
#'  level of the first split of the tree.
#'   }
#'  \item dendro - A dendrogram object of the decision tree output
#'  \item summaryMatrix - A K(labels) by L(features) matrix representation of
#'  the decision rules. Columns denote features and rows denote labels. Non-0
#'  values denote instances where a feature was used on a given label. Positive
#'  and negative values denote whether the values of the label for that feature
#'  were greater than or less than the decision threshold, respectively. The
#'  magnitude of Non-0 values denote the level at which the feature was used,
#'  where the first split has a magnitude of 1. Note, if reuse_features = TRUE,
#'  only the final usage of a feature for a given label is shown.
#'  \item prediction - A character vector of label of predictions of the
#'  training data using the final model. "MISSING" if label prediction was
#'  ambiguous.
#'  \item performance - A named list denoting the training performance of the
#'  model.
#'  \itemize{
#'   \item accuracy - (number correct/number of samples) for the whole set of
#'  samples.
#'   \item balAcc - mean sensitivity across all labels
#'   \item meanPrecision - mean precision across all labels
#'   \item correct - the number of correct predictions of each label
#'   \item sizes - the number of actual counts of each label
#'   \item sensitivity - the sensitivity of the prediciton of each label.
#'   \item precision - the precision of the prediciton of each label.
#'  }
#' }
#' @examples
#' library(M3DExampleData)
#' counts <- M3DExampleData::Mmus_example_list$data
#' # Subset 500 genes for fast clustering
#' counts <- as.matrix(counts[1501:2000, ])
#' # Cluster genes ans samples each into 10 modules
#' cm <- celda_CG(counts = counts, L = 10, K = 5, verbose = FALSE)
#' # Get features matrix and cluster assignments
#' factorized <- factorizeMatrix(counts, cm)
#' features <- factorized$proportions$cell
#' class <- clusters(cm)$z
#' # Generate Decision Tree
#' DecTree <- buildTreeHybrid(features,
#'     class,
#'     oneoffMetric = "modified F1",
#'     threshold = 1,
#'     consecutiveOneoff = FALSE)
#'
#' # Plot dendrogram
#' plotDendro(DecTree)
#' @export
buildTreeHybrid <- function(features,
    class,
    oneoffMetric = c("modified F1", "pairwise AUC"),
    threshold = 0.95,
    reuseFeatures = FALSE,
    altSplit = TRUE,
    consecutiveOneoff = TRUE) {

    if (ncol(features) != length(class)) {
        stop("Number of columns of features must equal length of class")
    }
    if (any(is.na(class))) {
        stop("NA class values")
    }
    if (any(is.na(features))) {
        stop("NA feature values")
    }

    # Match the oneoffMetric argument
    oneoffMetric <- match.arg(oneoffMetric)

    # Transpose features
    features <- t(features)

    # Set class to factor
    class <- as.factor(class)

    # Generate list of tree levels
    tree <- .generateTreeList(
        features,
        class,
        oneoffMetric,
        threshold,
        reuseFeatures,
        consecutiveOneoff)

    # Add alternative node for the solely down-regulated leaf
    if (altSplit) {
        tree <- .addAlternativeSplit(tree, features, class)
    }

    # Format tree output for plotting and generate summary statistics
    DTsummary <- .summarizeTree(tree, features, class)

    return(DTsummary)
}
