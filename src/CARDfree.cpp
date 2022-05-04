#include <iostream>
#include <fstream>
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
//              Extension of spatially informed deconvolution in a reference-free version:CARDfree                        //
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
SEXP CARDfree(SEXP XinputIn, SEXP UIn, SEXP WIn, SEXP phiIn, SEXP max_iterIn, SEXP epsilonIn, SEXP initV, SEXP initb, SEXP initSigma_e2, SEXP initLambda)
{    
    try {
        // read in the data
        arma::fmat Xinput = as<fmat>(XinputIn);
        arma::fmat U = as<fmat>(UIn);
        arma::fmat W = as<fmat>(WIn);
        float phi = as<float>(phiIn);
        int max_iter = Rcpp::as<int>(max_iterIn);
        float epsilon = as<float>(epsilonIn);
        arma::fmat V = as<fmat>(initV);
        arma::fvec b = as<fvec>(initb);
        float sigma_e2 = as<float>(initSigma_e2);
        arma::fvec lambda = as<fvec>(initLambda);
        // initialize some useful items
        int nSample = (int)Xinput.n_cols; // number of spatial sample points
        int mGene = (int)Xinput.n_rows; // number of genes in spatial deconvolution
        int k = (int)U.n_cols; // number of cell type
        arma::fmat L = zeros<fmat>(nSample,nSample);
        arma::fmat D = zeros<fmat>(nSample,nSample);
        arma::fmat V_old = zeros<fmat>(nSample,k);
        arma::fmat UtU = zeros<fmat>(k,k);
        arma::fmat VtV = zeros<fmat>(k,k);
        arma::fvec colsum_W = zeros<fvec>(nSample);
        arma::fmat UtX = zeros<fmat>(k,nSample);
        arma::fmat XtU = zeros<fmat>(nSample,k);
        arma::fmat UtXV = zeros<fmat>(k,k);
        arma::fmat temp = zeros<fmat>(k,k);
        arma::fmat part1 = zeros<fmat>(nSample,k);
        arma::fmat part2 = zeros<fmat>(nSample,k);
        arma::fvec updateV_k = zeros<fvec>(k);
        arma::fvec updateV_den_k = zeros<fvec>(k);
        arma::fvec vecOne = ones<fvec>( nSample);
        arma::fvec diag_UtU = zeros<fvec>(k);
        bool logicalLogL = FALSE;
        float obj = 0;
        float obj_old = 0;
        float normNMF = 0;
        float logX = 0;
        float logV = 0;
        float alpha = 1.0;
        float beta = nSample / 2.0;
        float logSigmaL2 = 0.0;
        float accu_L = 0.0;
        float trac_xxt = accu(Xinput % Xinput);
        
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
        logX = -(float)(mGene * nSample) * 0.5 * log(sigma_e2) - 0.5 * (float)(normNMF / sigma_e2);
        temp = (V.t() - b * vecOne.t()) * L * (V - vecOne * b.t());
        logV = - (float)(nSample) * 0.5 * sum(log(lambda )) - 0.5 * (sum(temp.diag() / lambda )); 
        logSigmaL2 = -(alpha + 1.0) * sum(log(lambda)) - sum(beta / lambda);
        obj_old = logX + logV + logSigmaL2;
        V_old = V;
        // iteration starts
        for(int i = 1; i <= max_iter; ++i) {
            logV = 0.0;  
            b = sum(V.t() * L, 1) / accu_L;
            lambda = (temp.diag() / 2.0 + beta ) / (float(nSample) / 2.0 + alpha + 1.0);  
            part1 = sigma_e2 * (D * V + phi * colsum_W * b.t());
            part2 = sigma_e2 * (phi * W * V + colsum_W * b.t());
            for(int nCT = 0; nCT < k; ++nCT){
                updateV_den_k = lambda(nCT) * (V.col(nCT) * diag_UtU(nCT) + (V * UtU.col(nCT) - V.col(nCT) * diag_UtU(nCT))) +  part1.col(nCT);
                updateV_k = (lambda(nCT) * XtU.col(nCT) + part2.col(nCT)) / updateV_den_k;
                V.col(nCT) %= updateV_k;
            }
            VtV = V.t() * V;
            //update U
            U %= (Xinput * V) / (U * VtV);
            UtX = U.t() * Xinput;
            XtU = UtX.t();
            UtU = U.t() * U;
            UtXV = UtX * V;

            normNMF = trac_xxt - 2.0 * trace(UtXV) + trace(UtU * VtV);
            sigma_e2 = normNMF / (float)(mGene * nSample);
            temp = (V.t() - b * vecOne.t()) * L * (V - vecOne * b.t());
            logX = -(float)(nSample * mGene) * 0.5 * log(sigma_e2) - 0.5 * (float)(normNMF / sigma_e2);
            logV = - (float)(nSample) * 0.5 * sum(log(lambda))- 0.5 * (sum(temp.diag() / lambda )); 
            logSigmaL2 = -(alpha + 1.0) * sum(log(lambda)) - sum(beta / lambda);
            obj = logX + logV + logSigmaL2;
            logicalLogL = (obj > obj_old) && (abs(obj - obj_old) * 2.0 / abs(obj + obj_old) < epsilon);
            if(isnan(obj) || (sqrt(accu((V - V_old) % (V - V_old)) / float(nSample * k))  < epsilon) || logicalLogL){
               if(i > 5){ // run at least 5 iterations 
               break;
           }
       }else{
        obj_old = obj;
        V_old = V;
       }
       }
       return List::create(Named("V") = V,
                           Named("B") = U,
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



