#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

SEXP _perplexityG(SEXP R_x, SEXP R_phi, SEXP R_psi, SEXP R_group)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  int nl = nlevels(R_group);
  
  // If the grouping variable is not a factor, throw an error
  if (!isFactor(R_group)) {
    error("The grouping argument must be a factor");
  }
  // If the length of the grouping variable and matrix do not match, throw an error
  if (LENGTH(R_group) != nr) {
    error("The length of the grouping argument must match the number of rows in the matrix.");
  }
  if (ncols(R_phi) != nc) {
    error("The R_phi and R_x must have the same number of colums.");
  }  
  if (nrows(R_phi) != nl) {
    error("R_phi must have the same number of rows as the number of levels in R_group.");
  }  
  if (nrows(R_psi) != nr) {
    error("The R_psi and R_x must have the same number of rows.");
  }  
  if (ncols(R_psi) != nl) {
    error("R_phi must have the same number of columns as the number of levels in R_group.");
  }  
  
  // Create pointers
  int *group = INTEGER(R_group);
  double *phi = REAL(R_phi);
  double *psi = REAL(R_psi);
  int *x = INTEGER(R_x);  
  
  // Make sure values are not NA and within the range of the number of rows
  for (i = 0; i < nr; i++) {
    if(group[i] == NA_INTEGER || group[i] < 0 || group[i] > nr) {
      error("Labels in group and pgroup must not be NA and must less than or equal to the number of rows in the matrix.");
    }
  }  
  
  // Allocate a variable for the return matrix
  
  
  double ans = 0;
  // Multiply the probabilties, log transform, and multiply against the counts to derive log(p(x))
  for (j = 0; j < nc; j++) {
    for (i = 0; i < nr; i++) {
      ans += x[j * nr + i] * log(phi[j * nl + (group[i]-1)] * psi[nr * (group[i]-1) + i]);
    }
  }
  
  SEXP R_ans = PROTECT(allocVector(REALSXP, 1));
  REAL(R_ans)[0] = ans;
  
  UNPROTECT(1);
  return(R_ans);
}  

