function [ u, v ] = panelvortex( xs, ys, xe, ye, x, y )
% Velocity induced by unit-strength constrant strength vortex panel.
%
% xs, ys Panel vortex start
% xe, ye Panel vortex end
% x, y   Point of interest

dxpan = xe - xs;
dypan = ye - ys;
lenpan = sqrt( dxpan.^2 + dypan.^2 );

dxp = dxpan ./ lenpan;
dyp = dypan ./ lenpan;

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


dtheta = acos( dxe .* dxs + dye .* dys );

mask = ( dxs .* dye - dxe .* dys ) < 0;
dtheta( mask ) = -dtheta( mask ); 


% Velocity in panel coordinates
xi = ( dtheta ) ./ ( 2.0 * pi );
eta = ( log( lenesq ) - log( lenssq ) ) ./ (4.0 * pi );

% Transform back to global coords.
u = xi .* dxp - eta .* dyp;
v = xi .* dyp + eta .* dxp;

end
