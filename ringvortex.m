function [ u, v ] = ringvortex( xn, rn, xm, rm )
% Velocity induced by unit-strength ring vortex around the x-axis.
%
% Vortex Element Methods for Fluid Dynamic Analysis of Engineering Systems
% R. I. Lewis, 1991, Cambridge University Press, p. 150
%
% n Subscript for the ring vortex
% m Subscript for the point of interest

% Non-dimensional coordinates
x = (xm - xn)./rn;
r = rm./rn;

A = x.^2 + ( r + 1.0 ).^2;
B = x.^2 + ( r - 1.0 ).^2;

% Matlab's elliptic integral routine is done in terms of M, where M=k^2
ksq = 4.0 * r ./ A;
[K,E] = ellipke( ksq );
% External routine, slightly slower than built-in.
% [K,E] = lellipke( ksq );

u =  - 1.0 ./ ( 2.0 * pi * rn .* sqrt( A ) ) .* ( K - ( 1.0 + 2.0 * ( r - 1.0 ) ./ B ) .* E );
v = (x./r) ./ ( 2.0 * pi * rn .* sqrt( A ) ) .* ( K - ( 1.0 + 2.0 * r ./ B ) .* E );

end
