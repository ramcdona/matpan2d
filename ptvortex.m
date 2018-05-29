function [ u, v ] = ptvortex( xn, yn, xm, ym )
% Velocity induced by unit-strength point vortex.
%
% Vortex Element Methods for Fluid Dynamic Analysis of Engineering Systems
% R. I. Lewis, 1991, Cambridge University Press, p. 150
%
% n Subscript for the point vortex
% m Subscript for the point of interest

dx = xm - xn;
dy = ym - yn;

rsq = dx.^2 + dy.^2;

reqzero = ( rsq == 0 );

u = ( 1.0 ./ ( 2.0 * pi ) ) .* ( dy ./ rsq );
v = ( -1.0 ./ ( 2.0 * pi ) ) .* ( dx ./ rsq );


v( reqzero ) = 0;

end
