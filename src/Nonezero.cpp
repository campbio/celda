#include <Rcpp.h>
using namespace Rcpp ;

//' get row and column indices of none zero elements in the matrix
////' 
////' @param R_counts A matrix
////' @return An integer matrix where each row is a row, column indices pair 
//// [[Rcpp::export]]
SEXP nonzero(NumericMatrix R_counts) {

    IntegerVector row(1);
		IntegerVector col(1);
		NumericVector val(1);

    int nR = R_counts.nrow();
    int nC = R_counts.ncol();
		double x;

    for (int c = 0; c < nC; c++) {
        for (int r = 0; r < nR; r++) {
            if (x = R_counts[c * nR + r] != 0) {
                row.push_back(r + 1);
                col.push_back(c + 1);
                val.push_back(x);
						}
				}
		}

    List res;
		res["row"] = row;
		res["col"] = col;
		res["val"] = val;

		return(res);
}

