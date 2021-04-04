#include <RcppEigen.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// Function to sum the counts for all cells that belong to a population together
//Performs the same function as .colSumByGroup but for sparse matrices
// [[Rcpp::export]]
Rcpp::NumericMatrix colSumByGroupSparse(
    const Eigen::MappedSparseMatrix<double> &counts,
    const IntegerVector &group)
{
  // Perform error checking
  if (counts.cols() != group.size())
  {
    stop("Length of 'group' must be equal to the number of columns in 'counts'.");
  }
  if (min(group) < 1 || max(group) > counts.cols())
  {
    stop("The entries in 'group' need to be between 1 and the number of columns in 'counts'.");
  }
//Need to check that group is a factor, maybe get levels

  NumericMatrix x(counts.rows(), max(group)); 

  for (int j = 0; j < counts.cols(); ++j)
  {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_)
    {
      x(i_.index(), group[j] - 1) += i_.value();
    }
  }

  return x;
}


// Function to sum the counts for all features that belong to a module together
//Performs the same function as .rowSumByGroup but for sparse matrices
// [[Rcpp::export]]
Rcpp::NumericMatrix rowSumByGroupSparse(
    const Eigen::MappedSparseMatrix<double> &counts,
    const IntegerVector &group)
{
  // Perform error checking
  if (counts.rows() != group.size())
  {
    stop("Length of 'group' must be equal to the number of rows in 'counts'.");
  }
  if (min(group) < 1 || max(group) > counts.rows())
  {
    stop("The entries in 'group' need to be between 1 and the number of rows in 'counts'.");
  }
  //Need to check that group is a factor, maybe get levels
  
  NumericMatrix x(max(group), counts.cols()); 
  
  for (int j = 0; j < counts.cols(); ++j)
  {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_)
    {
      x(group[i_.index()] - 1, j) += i_.value();
    }
  }
  
  return x;
}







// Function to change the counts for a cells previously assigned to one group
// that have now been assigned to another group. Performs the same function
// as .colSumByGroupChange but for sparse matrices
// [[Rcpp::export]]
Rcpp::NumericMatrix colSumByGroupChangeSparse(
    const Eigen::MappedSparseMatrix<double> &counts,
    const NumericMatrix &px,
    const IntegerVector &group,
    const IntegerVector &pgroup)
{
  // Perform error checking
  if (counts.cols() != group.size())
  {
    stop("Length of 'group' must be equal to the number of columns in 'counts'.");
  }
  if(group.size() != pgroup.size())
  {
    stop("Length of 'group' must equal 'pgroup'.");  
  }
  if (min(group) < 1 || max(group) > px.cols())
  {
    stop("The entries in 'group' need to be between 1 and the number of columns in 'ps'.");
  }
  if (min(pgroup) < 1 || max(pgroup) > px.cols())
  {
    stop("The entries in 'pgroup' need to be between 1 and the number of columns in 'ps'.");
  }
  if(px.rows() != counts.rows()) 
  {
    stop("'px' and 'counts' must have the same number of rows.");
  }

  NumericMatrix x = px; 
  
  for (int j = 0; j < counts.cols(); ++j)
  {
    if(group[j] != pgroup[j]) {
      for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_)
      {
        x(i_.index(), group[j] - 1) += i_.value();
        x(i_.index(), pgroup[j] - 1) -= i_.value();
      }
    }
  }
  
  return x;
}





// Function to change the counts for a features previously assigned to one group
// that have now been assigned to another group. Performs the same function
// as .rowSumByGroupChange but for sparse matrices
// [[Rcpp::export]]
Rcpp::NumericMatrix rowSumByGroupChangeSparse(
    const Eigen::MappedSparseMatrix<double> &counts,
    const NumericMatrix &px,
    const IntegerVector &group,
    const IntegerVector &pgroup)
{
  // Perform error checking
  if (counts.rows() != group.size())
  {
    stop("Length of 'group' must be equal to the number of rows in 'counts'.");
  }
  if(group.size() != pgroup.size())
  {
    stop("Length of 'group' must equal 'pgroup'.");
  }
  if (min(group) < 1 || max(group) > px.rows())
  {
    stop("The entries in 'group' need to be between 1 and the number of rows in 'px'.");
  }
  if (min(pgroup) < 1 || max(pgroup) > px.rows())
  {
    stop("The entries in 'pgroup' need to be between 1 and the number of rows in 'px'.");
  }
  if(px.cols() != counts.cols()) 
  {
    stop("'px' and 'counts' must have the same number of rows.");
  }
  
  NumericMatrix x = px;
  
  for (int j = 0; j < counts.cols(); ++j)
  {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_)
    {
      if(group[i_.index()] != pgroup[i_.index()]) {
        x(group[i_.index()] - 1, j) += i_.value();
        x(pgroup[i_.index()] - 1, j) -= i_.value();
      }
    }
  }
  
  return x;
}


