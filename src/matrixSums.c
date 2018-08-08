#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>



SEXP _rowSumByGroup(SEXP R_x, SEXP R_group)
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

SEXP _colSumByGroup(SEXP R_x, SEXP R_group)
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







SEXP _rowSumByGroupChange(SEXP R_x, SEXP R_px, SEXP R_group, SEXP R_pgroup)
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
  
  int g_ix;
  int pg_ix;
  for (i = 0; i < nr; i++) {
    if(pgroup[i] != group[i]) {
      for (j = 0; j < nc; j++) {
        pg_ix = (j * nl) + (pgroup[i] - 1);
        g_ix = (j * nl) + (group[i] - 1);
        px[pg_ix] -= x[j * nr + i];
        px[g_ix] += x[j * nr + i];      
      }
    }
  }

  return(R_px);
}




SEXP _colSumByGroupChange(SEXP R_x, SEXP R_px, SEXP R_group, SEXP R_pgroup)
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
        px[(group[j]-1) * nr + i] += x[j * nr + i];
        px[(pgroup[j]-1) * nr + i] -= x[j * nr + i];
      }
    }
  }

  return(R_px);
}




SEXP _mvAdd(SEXP R_m, SEXP R_v, SEXP R_by_row) {
  int i, j; 
  int nr = nrows(R_m);
  int nc = ncols(R_m);
  double *x = REAL(R_m); 
  double *v = REAL(R_v);
  int by_row = asLogical(R_by_row);
	
  // If the length of the vector and the length of the corresponding dimension of the matrix do not match, throw an error
  if  (!by_row) {
  	if (LENGTH(R_v) != nc ) {
    	error("The length of the vector must match the number of columns in the matrix");
    	}
  }
  if ( by_row) {
  	if (LENGTH(R_v) != nr ) {
  		error("The length of the vector must match the number of rows in the matrix");
  	}
  }

  // Allocate a variable for the return matrix
  SEXP R_ans;
  PROTECT(R_ans = allocMatrix(REALSXP, nr, nc));
  // Set a pointer to the return matrix and initialize the memory
  double *ans = REAL(R_ans);
  Memzero(ans, nr * nc);

  if (by_row) { 
  	// element-wise addition of the i-th row of the matrix to i-th element of the vector
  	for (j = 0; j < nc; j++) {
    	for (i = 0; i < nr; i++) {
      		ans[j * nr + i] = x[j * nr + i] + v[i];
  		}
  	}
  }
  if (!by_row) {
  	// element-wise addition of the j-th column of the matrix to j-th element of the vector
  	for (j = 0; j < nc; j++) {
    	for (i = 0; i < nr; i++) {
      		ans[j * nr + i] = x[j * nr + i] + v[j];
    	}
  	}
  }
  UNPROTECT(1);
  return(R_ans);
}





SEXP _mvMult(SEXP R_m, SEXP R_v, SEXP R_by_row)
{
  int i, j;
  int nr = nrows(R_m);
  int nc = ncols(R_m);
  double *x = REAL(R_m);
  double *v = REAL(R_v);
  int by_row = asLogical(R_by_row); 

  // If the length of the vector and the length of the corresponding dimension of the matrix do not match, throw an error
  if  (!by_row) {
  	if (LENGTH(R_v) != nc ) {
    	error("The length of the vector must match the number of columns in the matrix");
    	}
  }
  if ( by_row) {
  	if (LENGTH(R_v) != nr ) {
  		error("The length of the vector must match the number of rows in the matrix");
  	}
  }

  // Allocate a variable for the return matrix
  SEXP R_ans;
  PROTECT(R_ans = allocMatrix(REALSXP, nr, nc));
  // Set a pointer to the return matrix and initialize the memory
  double *ans = REAL(R_ans);
  Memzero(ans, nr * nc);

  if (by_row) { 
  	// element-wise multiplication of the i-th row of the matrix to i-th element of the vector
  	for (j = 0; j < nc; j++) {
    	for (i = 0; i < nr; i++) {
      		ans[j * nr + i] = x[j * nr + i] * v[i];
  		}
  	}
  }
  if (!by_row) {
  	// element-wise multiplication of the j-th column of the matrix to j-th element of the vector
  	for (j = 0; j < nc; j++) {
    	for (i = 0; i < nr; i++) {
      		ans[j * nr + i] = x[j * nr + i] * v[j];
    	}
  	}
  }
  UNPROTECT(1);
  return(R_ans);
}





SEXP _rowSumByGroup_numeric(SEXP R_x, SEXP R_group)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  double *x = REAL(R_x);

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
  PROTECT(R_ans = allocMatrix(REALSXP, nl, nc));
  // Set a pointer to the return matrix and initialize the memory
  double *ans = REAL(R_ans);
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



SEXP _colSumByGroup_numeric(SEXP R_x, SEXP R_group)
{
  int i, j;
  int nr = nrows(R_x);
  int nc = ncols(R_x);
  double *x = REAL(R_x);

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
  PROTECT(R_ans = allocMatrix(REALSXP, nr, nl));
  // Set a pointer to the return matrix and initialize the memory
  double *ans = REAL(R_ans);
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



