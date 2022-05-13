#include <iostream>
#include <fstream>
#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>
#include <Rcpp.h>

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS


//*******************************************************************//
//              spatially informed deconvolution:CARD                        //
//*******************************************************************//
//' SpatialDeconv function based on Conditional Autoregressive model
//' @param XinputIn The input of normalized spatial data
//' @param UIn The input of cell type specific basis matrix B
//' @param WIn The constructed W weight matrix from Gaussian kernel
//' @param phiIn The phi value
//' @param max_iterIn Maximum iterations
//' @param epsilonIn epsilon for convergence 
//' @param initV Initial matrix of cell type compositions V
//' @param initb Initial vector of cell type specific intercept
//' @param initSigma_e2 Initial value of residual variance
//' @param initLambda Initial vector of cell type sepcific scalar. 
//'
//' @return A list
//'
//' @export
// [[Rcpp::export]]
SEXP CARDref(SEXP XinputIn, SEXP UIn, SEXP WIn, SEXP phiIn, SEXP max_iterIn, SEXP epsilonIn, SEXP initV, SEXP initb, SEXP initSigma_e2, SEXP initLambda)
{    
    try {
        // read in the data
        arma::mat Xinput = as<mat>(XinputIn);
        arma::mat U = as<mat>(UIn);
        arma::mat W = as<mat>(WIn);
        double phi = as<double>(phiIn);
        int max_iter = Rcpp::as<int>(max_iterIn);
        double epsilon = as<double>(epsilonIn);
        arma::mat V = as<mat>(initV);
        arma::vec b = as<vec>(initb);
        double sigma_e2 = as<double>(initSigma_e2);
        arma::vec lambda = as<vec>(initLambda);
        // initialize some useful items
        int nSample = (int)Xinput.n_cols; // number of spatial sample points
        int mGene = (int)Xinput.n_rows; // number of genes in spatial deconvolution
        int k = (int)U.n_cols; // number of cell type
        arma::mat L = zeros<mat>(nSample,nSample);
        arma::mat D = zeros<mat>(nSample,nSample);
        arma::mat V_old = zeros<mat>(nSample,k);
        arma::mat UtU = zeros<mat>(k,k);
        arma::mat VtV = zeros<mat>(k,k);
        arma::vec colsum_W = zeros<vec>(nSample);
        arma::mat UtX = zeros<mat>(k,nSample);
        arma::mat XtU = zeros<mat>(nSample,k);
        arma::mat UtXV = zeros<mat>(k,k);
        arma::mat temp = zeros<mat>(k,k);
        arma::mat part1 = zeros<mat>(nSample,k);
        arma::mat part2 = zeros<mat>(nSample,k);
        arma::vec updateV_k = zeros<vec>(k);
        arma::vec updateV_den_k = zeros<vec>(k);
        arma::vec vecOne = ones<vec>( nSample);
        arma::vec diag_UtU = zeros<vec>(k);
        bool logicalLogL = FALSE;
        double obj = 0;
        double obj_old = 0;
        double normNMF = 0;
        double logX = 0;
        double logV = 0;
        double alpha = 1.0;
        double beta = nSample / 2.0;
        double logSigmaL2 = 0.0;
        double accu_L = 0.0;
        double trac_xxt = accu(Xinput % Xinput);
        
        // initialize values
        // constant matrix caculations for increasing speed 
        UtX = U.t() * Xinput;
        XtU = UtX.t();
        colsum_W = sum(W,1);
        D =  diagmat(colsum_W);// diagnol matrix whose entries are column
        L = D -  phi*W; // graph laplacian
        accu_L = accu(L);
        UtXV = UtX * V;
        VtV = V.t() * V;
        UtU = U.t() * U;
        diag_UtU = UtU.diag();
        // calculate initial objective function 
        normNMF = trac_xxt - 2.0 * trace(UtXV) + trace(UtU * VtV);
        logX = -(double)(mGene * nSample) * 0.5 * log(sigma_e2) - 0.5 * (double)(normNMF / sigma_e2);
        temp = (V.t() - b * vecOne.t()) * L * (V - vecOne * b.t());
        logV = - (double)(nSample) * 0.5 * sum(log(lambda )) - 0.5 * (sum(temp.diag() / lambda )); 
        logSigmaL2 = -(alpha + 1.0) * sum(log(lambda)) - sum(beta / lambda);
        obj_old = logX + logV + logSigmaL2;
        V_old = V;
        // iteration starts
        for(int i = 1; i <= max_iter; ++i) {
            logV = 0.0;  
            b = sum(V.t() * L, 1) / accu_L;
            lambda = (temp.diag() / 2.0 + beta ) / (double(nSample) / 2.0 + alpha + 1.0);  
            part1 = sigma_e2 * (D * V + phi * colsum_W * b.t());
            part2 = sigma_e2 * (phi * W * V + colsum_W * b.t());
            for(int nCT = 0; nCT < k; ++nCT){
                updateV_den_k = lambda(nCT) * (V.col(nCT) * diag_UtU(nCT) + (V * UtU.col(nCT) - V.col(nCT) * diag_UtU(nCT))) +  part1.col(nCT);
                updateV_k = (lambda(nCT) * XtU.col(nCT) + part2.col(nCT)) / updateV_den_k;
                V.col(nCT) %= updateV_k;
            }
            UtXV = UtX * V;
            VtV = V.t() * V;
            normNMF = trac_xxt - 2.0 * trace(UtXV) + trace(UtU * VtV);
            sigma_e2 = normNMF / (double)(mGene * nSample);
            temp = (V.t() - b * vecOne.t()) * L * (V - vecOne * b.t());
            logX = -(double)(nSample * mGene) * 0.5 * log(sigma_e2) - 0.5 * (double)(normNMF / sigma_e2);
            logV = - (double)(nSample) * 0.5 * sum(log(lambda))- 0.5 * (sum(temp.diag() / lambda )); 
            logSigmaL2 = -(alpha + 1.0) * sum(log(lambda)) - sum(beta / lambda);
            obj = logX + logV + logSigmaL2;
            logicalLogL = (obj > obj_old) && (abs(obj - obj_old) * 2.0 / abs(obj + obj_old) < epsilon);
            if(isnan(obj) || (sqrt(accu((V - V_old) % (V - V_old)) / double(nSample * k))  < epsilon) || logicalLogL){
               if(i > 5){ // run at least 5 iterations 
               break;
           }
       }else{
        obj_old = obj;
        V_old = V;
       }
       }
       return List::create(Named("V") = V,
                           Named("sigma_e2") = sigma_e2,
                           Named("lambda") = lambda,
                           Named("b") = b,
                           Named("Obj") = obj);
        }//end try 
        catch (std::exception &ex)
        {
            forward_exception_to_r(ex);
        }
        catch (...)
        {
            ::Rf_error("C++ exception (unknown reason)...");
        }
        return R_NilValue;
} // end funcs



