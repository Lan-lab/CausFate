#include <Rcpp.h>
#include <algorithm>
#include <random>

//' @useDynLib causfate

// [[Rcpp::export]]
Rcpp::DataFrame dataFrameShuffle(Rcpp::DataFrame df)
{
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(df.begin(),df.end(),g);
    return df;
}
