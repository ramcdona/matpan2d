function [ u, v ] = ringsubvortex2( xm, rm, xs, rs, xe, re, varargin )
% Velocity induced by unit-strength ring vortex around the x-axis.
%

n = 16;

xepts = linspace( xs, xe, n + 1 );
repts = linspace( rs, re, n + 1 );

xn = ( xepts( 1:end - 1 ) + xepts( 2:end ) ) * 0.5;
rn = ( repts( 1:end - 1 ) + repts( 2:end ) ) * 0.5;

mask = abs(xn-xm) < 1e-5 & abs(rn-rm) < 1e-5;

[u,v] = ringvortex( xn(~mask), rn(~mask), xm, rm );

u = sum(u) / n;
v = sum(v) / n;

end
