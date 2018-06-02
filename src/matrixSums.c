#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>

SEXP rowSumByGroup(SEXP R_x, SEXP R_group)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  int *x = INTEGER(R_x);

  // If the grouping variable is not a factor, throw an error
  if (!isFactor(R_group)) {
    error("The grouping argument must be a factor");
  }

  int *group = INTEGER(R_group);
  int nl = nlevels(R_group);
  // If the sizes of the grouping variable and matrix do not match, throw an error
  if (LENGTH(R_group) != nr) {
    error("The length of the grouping argument must match the number of rows in the matrix");
  }

  // Allocate a variable for the return matrix
  SEXP R_ans;
  PROTECT(R_ans = allocMatrix(INTSXP, nl, nc));
  // Set a pointer to the return matrix and initialize the memory
  int *ans = INTEGER(R_ans);
  Memzero(ans, nl * nc);

  // Sum the totals for each element of the 'group' variable
  // Note: columns are iterated over before rows because the compiler appears to store expressions like
  //       'j * nr' in a temporary variable (as they do not change within inner loop);
  //       swapping the order of the outer and inner loops slows down the code ~10X
  for (j = 0; j < nc; j++) {
    for (i = 0; i < nr; i++) {
      ans[j * nl + group[i] - 1] += x[j * nr + i];
    }
  }

  UNPROTECT(1);
  return(R_ans);
}

SEXP colSumByGroup(SEXP R_x, SEXP R_group)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  int *x = INTEGER(R_x);

  // If the grouping variable is not a factor, throw an error
  if (!isFactor(R_group)) {
    error("The grouping argument must be a factor");
  }

  int *group = INTEGER(R_group);
  int nl = nlevels(R_group);
  // If the sizes of the grouping variable and matrix do not match, throw an error
  if (LENGTH(R_group) != nc) {
    error("The length of the grouping argument must match the number of columns in the matrix");
  }
  
  // Allocate a variable for the return matrix
  SEXP R_ans;
  PROTECT(R_ans = allocMatrix(INTSXP, nr, nl));
  // Set a pointer to the return matrix and initialize the memory
  int *ans = INTEGER(R_ans);
  Memzero(ans, nr * nl);

  // Sum the totals for each element of the 'group' variable
  // Note: columns are iterated over before rows because the compiler appears to store expressions like
  //       'j * nr' in a temporary variable (as they do not change within inner loop);
  //       swapping the order of the outer and inner loops slows down the code ~10X
  for (j = 0; j < nc; j++) {
    for (i = 0; i < nr; i++) {
      ans[(group[j] - 1) * nr + i] += x[j * nr + i];
    }
  }

  UNPROTECT(1);
  return(R_ans);
}







SEXP rowSumByGroupChange(SEXP R_x, SEXP R_px, SEXP R_group, SEXP R_pgroup)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  int *x = INTEGER(R_x);
  int *px = INTEGER(R_px);    
  int *group = INTEGER(R_group);
  int *pgroup = INTEGER(R_pgroup);
  
  // If the grouping variable is not a factor, throw an error
  if (!isFactor(R_group) | !isFactor(R_pgroup)) {
    error("The grouping arguments must be factors");
  }
  int nl = nlevels(R_group);
  int nlp = nlevels(R_pgroup);
  if(nl != nlp || nl != nrows(R_px)) {
    error("group and pgroup must have the same number of levels equal to row number of px");
  }
  
  if(nc != ncols(R_px)) {
    error("x and the previously summed matrix, px, must have the same number of columns.");
  }
  if(length(R_group) != length(R_pgroup) || length(R_group) != nr) {
    error("group label and previous group label must be the same length as the number of rows in x.");
  }
  
  // Create vector containing indices where group is different than pgroup
  // First, calculate the number of differences in group labels
  int ndiff = 0;
  for(i = 0; i < nr; i++) {
    
    if(pgroup[i] != group[i]) {
      ndiff++;
    }
  }
  // Second, add the index of the differences to a vector
  int group_diff_ix[ndiff-1];
  j = 0;
  for(i = 0; i < nr; i++) {
    if(group[i] != pgroup[i]) {
      group_diff_ix[j] = i;
      j++;
    }
  }
  
  // Sum the totals for each element of the 'group' variable
  // Note: columns are iterated over before rows because the compiler appears to store expressions like
  //       'j * nr' in a temporary variable (as they do not change within inner loop);
  //       swapping the order of the outer and inner loops slows down the code ~10X
  int row_ix;
  for (j = 0; j < nc; j++) {
    for (i = 0; i < ndiff; i++) {
      row_ix = group_diff_ix[i];
      px[j * nl + (pgroup[row_ix] - 1)] -= x[j * nr + row_ix];
      px[j * nl + (group[row_ix] - 1)] += x[j * nr + row_ix];      
    }
  }

  return(R_px);
}



SEXP colSumByGroupChange(SEXP R_x, SEXP R_px, SEXP R_group, SEXP R_pgroup)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  int *x = INTEGER(R_x);
  int *px = INTEGER(R_px);  
  int *group = INTEGER(R_group);
  int *pgroup = INTEGER(R_pgroup);
  
  // If the grouping variable is not a factor, throw an error
  if (!isFactor(R_group) | !isFactor(R_pgroup)) {
    error("The grouping arguments must be factors");
  }
  int nl = nlevels(R_group);
  int nlp = nlevels(R_pgroup);
  if(nl != nlp || nl != ncols(R_px)) {
    error("group and pgroup must have the same number of levels equal to column number of px");
  }

  if(nr != nrows(R_px)) {
    error("x and the previously summed matrix, pxc must have the same number of rows");
  }
  if(length(R_group) != length(R_pgroup) || length(R_group) != nc) {
    error("group label and previous group label must be the same length as the number of columns in x.");
  }
  
  // Sum the totals for each element of the 'group' variable,
  // But only where the group label is different than the previous label
  // Note: columns are iterated over before rows because the compiler appears to store expressions like
  //       'j * nr' in a temporary variable (as they do not change within inner loop);
  //       swapping the order of the outer and inner loops slows down the code ~10X
  for (j = 0; j < nc; j++) {
    if(group[j] != pgroup[j]) {
      for (i = 0; i < nr; i++) {
        px[group[j] * nr + i - 1] += x[j * nr + i];
        px[pgroup[j] * nr + i - 1] -= x[j * nr + i];
      }
    }
  }

  return(R_px);
}

