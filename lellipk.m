
%
% lellipf( k, errtol)
%
% Inputs:
%
%   phi     Input angle vector size 1 or 1xN.
%   k       Input parameter vector size 1 or 1xN.
%   errtol  Error tolerance for Carlson's algorithms.
%
% Matlab function to compute Legendre's (incomplete) elliptic integral 
% K(k).  Uses a vectorized implementation of Carlson's Duplication Algorithms 
% for symmetric elliptic integrals as found in "Computing Elliptic 
% Integrals by Duplication," by B. C. Carlson, Numer. Math. 33, 1-16 (1979)
% and also found in ACM TOMS Algorithm 577.  Section 4 in the paper cited
% here describes how to convert between the symmetric elliptic integrals
% and Legendre's elliptic integrals.
%
% Returns NaN's for any argument values outside input range.
%

function f = lellipk( k, errtol)

kvec = k;

onesvec = ones(1, length(kvec) );
y = onesvec - kvec.*kvec;
f = rf( 0.0,  y, onesvec, errtol);
