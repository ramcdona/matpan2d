function [xep, rep, rad, Vex] = setupsphere( ncirc, radius, xcen )

% Construct a circle
thetacirc = linspace( pi, 0, ncirc );

% Panel endpoints
xep = radius * cos( thetacirc ) + xcen;
rep = radius * sin( thetacirc );


rad = radius * ones( 1, ncirc - 1 );

% Center points
xcp = 0.5 * ( xep(2:end) + xep(1:end-1) ) - xcen;
rcp = 0.5 * ( rep(2:end) + rep(1:end-1) );

thetacp = atan2( rcp, xcp );

Vex = 1.5 * sin( thetacp );

end
