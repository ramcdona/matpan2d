function [ u, v ] = rayvortex( xs, ys, xv, yv, x, y )
% Velocity induced by unit-strength semi-infinite constrant strength vortex panel.
%
% xs, ys Ray vortex start
% xv, yv Ray vortex direction vector
% x, y   Point of interest

big = 100;

% Make sure direction vector is unit vector.
lenv = sqrt( xv.^2 + yv.^2 );
xv = xv ./ lenv;
yv = yv ./ lenv;

xe = xs + big * xv;
ye = ys + big * yv;

dxstart = x - xs;
dystart = y - ys;
lenssq = dxstart.^2 + dystart.^2;
lens = sqrt( lenssq );
dxs = dxstart ./ lens;
dys = dystart ./ lens;

dxend = x - xe;
dyend = y - ye;
lenesq = dxend.^2 + dyend.^2;
lene = sqrt( lenesq );
dxe = dxend ./ lene;
dye = dyend ./ lene;


% dtheta = acos( dxe .* dxs + dye .* dys );
% mask = ( dxs .* dye - dxe .* dys ) < 0;
% dtheta( mask ) = -dtheta( mask ); 

dtheta = pi - acos( xv .* dxs + yv .* dys );
mask = ( dxs .* yv - xv .* dys ) < 0;
dtheta( mask ) = -dtheta( mask ); 



% Velocity in panel coordinates
xi = ( dtheta ) ./ ( 2.0 * pi );
eta = ( log( lenesq ) - log( lenssq ) ) ./ (4.0 * pi );

% Transform back to global coords.
u = xi .* xv - eta .* yv;
v = xi .* yv + eta .* xv;

end
