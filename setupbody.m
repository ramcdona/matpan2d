function [xep, rep, rad, Vex] = setupbody( ncirc, nstraight, rad, len1, len2, xnose )

thetacirc = linspace( pi, pi/2, ncirc );

xep1 = rad * cos( thetacirc ) + rad;
rep1 = rad * sin( thetacirc );

xep2 = linspace( rad, rad + len1, nstraight );
rep2 = rad * ones( 1, nstraight );

xep3 = linspace( rad + len1, rad + len1 + len2, nstraight);
rep3 = linspace( rad, 0, nstraight);

xep = [xep1 xep2(2:end) xep3(2:end)] + xnose;
rep = [rep1 rep2(2:end) rep3(2:end)];

rad = [rad * ones( 1, ncirc - 1 ) inf(1, 2 * ( nstraight - 1 ) )];

Vex = nan( 1, length( xep ) - 1 );

end
