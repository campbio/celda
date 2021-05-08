#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp)]]

//Contains a version that is more simple to implement and understand as it matches the equations more directly 
// However it is many times slower. It is useful to use as a sanity check when modifying the faster function. 
// [[Rcpp::export]]
NumericVector cG_calcGibbsProbY_Simple(const IntegerMatrix counts,
	IntegerVector nGbyTS,
	IntegerMatrix nTSbyC,
	IntegerVector nbyTS,
	IntegerVector nbyG, 
	const IntegerVector y,  
	const int L,  
	const int index, 
	const double gamma, 
	const double beta, 
	const double delta) {
  int index0 = index - 1;
  int current_y = y[index0] - 1;
  
  NumericVector probs(L);

  nGbyTS[current_y] -= 1;
  nbyTS[current_y] -= nbyG[index0];
  nTSbyC(current_y,_) = nTSbyC(current_y,_) - counts(index0,_);

  for (int i = 0; i < L; i++) {
    nGbyTS[i] += 1;
    nbyTS[i] += nbyG[index0];    
    nTSbyC(i,_) = nTSbyC(i,_) + counts(index0,_);

    probs[i] += sum(lgamma(nGbyTS + gamma));
    probs[i] += sum(lgamma(nTSbyC + beta));	  
    probs[i] += sum(lgamma(nGbyTS * delta));
    probs[i] -= sum(lgamma(nbyTS + (nGbyTS * delta)));
    
    nGbyTS[i] -= 1;
    nbyTS[i] -= nbyG[index0];        
    nTSbyC(i,_) = nTSbyC(i,_) - counts(index0,_);
  }
  
  nGbyTS[current_y] += 1;
  nbyTS[current_y] += nbyG[index0];  
  nTSbyC(current_y,_) = nTSbyC(current_y,_) + counts(index0,_);
  
  return(probs);
}



// [[Rcpp::export]]
NumericVector cG_CalcGibbsProbY_ori(const int index,
	const IntegerMatrix& counts,
	const IntegerMatrix& nTSbyC,
	const IntegerVector& nbyTS,
	const IntegerVector& nGbyTS,	
	const IntegerVector& nbyG,
	const IntegerVector& y, 
	const int L, 
	const int nG,
	const NumericVector& lg_beta,
	const NumericVector& lg_gamma,
	const NumericVector& lg_delta,
	const double delta) {

  int index0 = index - 1;
  int current_y = y[index0] - 1;
  int i, j, k;
  
  NumericVector probs(L);
  NumericVector nTSbyC_prob1(L);
  NumericVector nTSbyC_prob2(L);  
  
  // Calculate probabilities related to the "n.TS.by.C" part of equation one time up front
  // The first vector represents when the current feature is added to that module
  // The second vector represents when the current feature is NOT added to that module
  for (int col = 0; col < counts.ncol(); col++) {
    j = col * L + current_y; // Index for the current module in the n.TS.by.C matrix
    k = col * nG + index0; // Index for the current feature in counts matrix
	for (int row = 0; row < L; row++) {
	  if (row == current_y) {
		nTSbyC_prob1[row] += lg_beta[nTSbyC[j] - counts[k]];
		nTSbyC_prob2[row] += lg_beta[nTSbyC[j]];
	  } else {
		nTSbyC_prob1[row] += lg_beta[nTSbyC[col * L + row]];
		nTSbyC_prob2[row] += lg_beta[nTSbyC[col * L + row] + counts[k]];
	  }
	}
  }

  // Calculate the probabilities for each module
  // If statements determine whether to add or subtract counts from each probability
  for (i = 0; i < L; i++) {
	for(j = 0; j < L; j++) {
	  if((i == j) & (i != current_y)) {
		probs[i] += lg_gamma[nGbyTS[j] + 1];
		probs[i] += nTSbyC_prob2[j];
		probs[i] += lg_delta[nGbyTS[j] + 1];
		probs[i] -= lgamma(nbyTS[j] + nbyG[index0] + ((nGbyTS[j] + 1) * delta));
	  } else if ((j == current_y) & (i != current_y)) {
		probs[i] += lg_gamma[nGbyTS[j] - 1];;
		probs[i] += nTSbyC_prob1[j];		  
		probs[i] += lg_delta[nGbyTS[j] - 1];
		probs[i] -= lgamma(nbyTS[j] - nbyG[index0] + ((nGbyTS[j] - 1) * delta));
	  } else {
		probs[i] += lg_gamma[nGbyTS[j]];;
		probs[i] += lg_delta[nGbyTS[j]];
		probs[i] -= lgamma(nbyTS[j] + (nGbyTS[j] * delta));
		
		if(j == current_y) {
		  probs[i] += nTSbyC_prob2[j];
		} else {
		  probs[i] += nTSbyC_prob1[j];          
		}
	  } 
	}
  }

  return(probs);
}


// [[Rcpp::export]]
NumericVector cG_CalcGibbsProbY_fastRow(const int index,
	const IntegerMatrix& counts,
	const IntegerMatrix& nTSbyC,
	const IntegerVector& nbyTS,
	const IntegerVector& nGbyTS,	
	const IntegerVector& nbyG,
	const IntegerVector& y, 
	const int L, 
	const int nG,
	const NumericVector& lg_beta,
	const NumericVector& lg_gamma,
	const NumericVector& lg_delta,
	const double delta) {

  int index0 = index - 1;
  int current_y = y[index0] - 1;
  int i;
  
  NumericVector probs(L);
  
  // Calculate probabilities related to the "n.TS.by.C" part of equation one time up front
  // The first case of if statement represents when the current feature is already added to that module
  // The second case represents when the current feature is NOT YET added to that module
  for (i = 0; i < L; i++) {
    if (i == current_y ) {
      for (int col = 0; col < counts.ncol(); col++) {
        probs[i] += lg_beta[nTSbyC(i, col)];
	probs[i] -= lg_beta[nTSbyC(i, col) - counts(index0, col)];
      }
    } else {
      for (int col = 0; col < counts.ncol(); col++) {
        probs[i] += lg_beta[nTSbyC(i, col) + counts(index0, col)];
	probs[i] -= lg_beta[nTSbyC(i, col)];
      }
    }
  }

  // Calculate the probabilities for each module
  // If statements determine whether to add or subtract counts from each probability
  for (i = 0; i < L; i++) {
    if (i == current_y) {
      probs[i] += lg_gamma[nGbyTS[i]];
      probs[i] -= lg_gamma[nGbyTS[i] - 1];
      probs[i] += lg_delta[nGbyTS[i]];
      probs[i] -= lg_delta[nGbyTS[i] - 1];
      probs[i] += lgamma(nbyTS[i] - nbyG[index0] + (nGbyTS[i] - 1) * delta);
      probs[i] -= lgamma(nbyTS[i] + nGbyTS[i] * delta);
    } else {
      probs[i] += lg_gamma[nGbyTS[i] + 1];
      probs[i] -= lg_gamma[nGbyTS[i]];
      probs[i] += lg_delta[nGbyTS[i] + 1];
      probs[i] -= lg_delta[nGbyTS[i]];
      probs[i] += lgamma(nbyTS[i] + nGbyTS[i] * delta);
      probs[i] -= lgamma(nbyTS[i] + nbyG[index0] + (nGbyTS[i] + 1) * delta);
    }
  }

  return(probs);
}


// [[Rcpp::export]]
NumericVector cG_CalcGibbsProbY(const int index,
	const NumericVector& counts,
	const IntegerMatrix& nTSbyC,
	const IntegerVector& nbyTS,
	const IntegerVector& nGbyTS,	
	const IntegerVector& nbyG,
	const IntegerVector& y, 
	const int L, 
	const int nG,
	const NumericVector& lg_beta,
	const NumericVector& lg_gamma,
	const NumericVector& lg_delta,
	const double delta) {

  int index0 = index - 1;
  int current_y = y[index0] - 1;
  //int current_y = y;
  int i;
  int j,k;
  
  NumericVector probs(L);
  
  // Calculate probabilities related to the "n.TS.by.C" part of equation one time up front
  // The first case of if statement represents when the current feature is already added to that module
  // The second case represents when the current feature is NOT YET added to that module
  for (int col = 0; col < counts.length(); col++) {
    //k = col * nG + index0; // Index for the current feature in counts matrix, used when the whole matrix was passed to this function rather than just the row
    k = col; // Index for the current feature in counts matrix
    for (i = 0; i < L; i++) {
      j = col * L + i; // Index for the current module in the n.TS.by.C matrix
      if (i == current_y) {
        probs[i] += lg_beta[nTSbyC[j]];
	probs[i] -= lg_beta[nTSbyC[j] - counts[k]];
      } else {
        probs[i] += lg_beta[nTSbyC[j] + counts[k]];
	probs[i] -= lg_beta[nTSbyC[j]];
      }
    }
  }


  // Calculate the probabilities for each module
  // If statements determine whether to add or subtract counts from each probability
  for (i = 0; i < L; i++) {
    if (i == current_y) {
      probs[i] += lg_gamma[nGbyTS[i]];
      probs[i] -= lg_gamma[nGbyTS[i] - 1];
      probs[i] += lg_delta[nGbyTS[i]];
      probs[i] -= lg_delta[nGbyTS[i] - 1];
      probs[i] += lgamma(nbyTS[i] - nbyG[index0] + (nGbyTS[i] - 1) * delta);
      probs[i] -= lgamma(nbyTS[i] + nGbyTS[i] * delta);
    } else {
      probs[i] += lg_gamma[nGbyTS[i] + 1];
      probs[i] -= lg_gamma[nGbyTS[i]];
      probs[i] += lg_delta[nGbyTS[i] + 1];
      probs[i] -= lg_delta[nGbyTS[i]];
      probs[i] += lgamma(nbyTS[i] + nGbyTS[i] * delta);
      probs[i] -= lgamma(nbyTS[i] + nbyG[index0] + (nGbyTS[i] + 1) * delta);
    }
  }

  return(probs);
}
