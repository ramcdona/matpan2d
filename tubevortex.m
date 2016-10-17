function [ u, v ] = tubevortex( xn, rn, xm, rm )
% Velocity induced by unit-strength semi-infinite (downstream) tube vortex
% around the x-axis.
%
% Vortex Element Methods for Fluid Dynamic Analysis of Engineering Systems
% R. I. Lewis, 1991, Cambridge University Press, p. 169
%
% n Subscript for the ring vortex
% m Subscript for the point of interest

% Non-dimensional coordinates
x = (xm - xn)./rn;
r = rm./rn;

reqone = ( r == 1 );

A = zeros( size(r) );
A( r < 1 ) = pi;
A( reqone ) = pi/2.0;

C = x.^2 + ( r + 1.0 ).^2;

% Matlab's elliptic integral routine is done in terms of M, where M=k^2
ksq = 4.0 * r ./ C;
n = ( 4.0 * r ) ./ ( r + 1.0 ).^2;

% Built-in Matlab routines
% [K,E] = ellipke( ksq );
% PI = ellipticPi( n, ksq );

% External routine (much faster for Pi)
[K,E,PI]= lellipkepi( ksq, n );

v = 2.0 ./ ( pi .* ksq .* sqrt( C ) ) .* ( E - K .* (1.0 - ksq ./ 2 ) );

u = zeros( size (v) );
u( ~reqone ) = 1.0 ./ ( 2.0 * pi ) .* ( A(~reqone) + ( x(~reqone) ./ sqrt( C(~reqone) ) ) .* ( K(~reqone) - PI(~reqone) .* ( r(~reqone) - 1.0 ) ./ ( r(~reqone) + 1.0 ) ) );
u( reqone ) = 0.25 + x(reqone) .* K(reqone) ./ ( 2.0 * pi .* sqrt( x(reqone).^2 + 4.0 ) );

u( isnan(u) ) = 0.25;
v( isnan(v) ) = -(1+pi/2);
v( v < -(1+pi/2) ) = -(1+pi/2);

end
