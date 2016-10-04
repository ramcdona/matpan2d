function [xep, rep, rad, Vex] = setupellipsoid( ncirc, rada, radb, xcen )

% Construct an ellipse
thetacirc = linspace( pi, 0, ncirc );

% Panel endpoints
xep = rada * cos( thetacirc ) + xcen;
rep = radb * sin( thetacirc );

% Center points
xcp = 0.5 * ( xep(2:end) + xep(1:end-1) );
rcp = 0.5 * ( rep(2:end) + rep(1:end-1) );

% Radius vector
rad = rada.^2*radb.^2.*( xcp.^2./rada.^4 + rcp.^2./radb.^4 ).^1.5;

[~, u, v, w] = ellipsoidflow( xcp, rcp, zeros(1,ncirc-1), [1 0 0], rada, radb, radb );
Vex = sqrt(u.^2 + v.^2 + w.^2);

end
