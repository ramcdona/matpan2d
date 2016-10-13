
%
% lellipke( k, errtol)
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

function [K, varargout ] = lellipke( m, varargin )

if nargin < 2
    errtol = 1e-10;
else
    errtol = varargin{1};
end

onesvec = ones( size(m) );
zerosvec = zeros( size(m) );

y = 1.0 - m;

K = rf( zerosvec,  y, onesvec, errtol);

if nargout > 1
    varargout{1} = K - m .* rd( zerosvec, y, onesvec, errtol)/3.0;
end
