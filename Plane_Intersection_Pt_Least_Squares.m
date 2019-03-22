%
% HELPER SCRIPT TO CALCULATE THE FORMULAS FOR LEAST SQUARE FITTING A POINT
% TO SEVERAL PLANES (FIND THE POINT WHICH IS CLOSEST TO ALL PLANES).
%
% This script uses an Least Squares estimator using the (X'*X)^-1 *X'*y
% 23 January 2015
%

clear all;


% X'*X is a 3x3 matrix, each element is a sum of normal products.
syms sum_nx_nx sum_nx_ny sum_nx_nz sum_ny_ny sum_ny_nz sum_nz_nz

syms sum_nx_d sum_ny_d sum_nz_d



XTX = [sum_nx_nx sum_nx_ny sum_nx_nz; 
    sum_nx_ny sum_ny_ny sum_ny_nz;
    sum_nx_nz sum_ny_nz sum_nz_nz];

XTy = [-sum_nx_d; -sum_ny_d; -sum_nz_d];

XTXInv = inv(XTX);

pOpt = XTXInv*XTy;

simple(pOpt)