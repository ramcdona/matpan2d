function [xep, rep, rad, Vex] = setupVKTduct( naf, chord, alpha, xoff, roff, epsilon, kappa, tau, flipaf, spacing )

%[ xep, rep ] = naca4( dig1, dig2, dig34, naf, spacing );


[xep, rep, cp_vkt] = vkt_airfoil( naf, alpha*pi/180, epsilon, kappa, tau );
xep = xep + 0.5;


% Make sure airfoil is _exactly_ closed.
xep(end) = xep(1);
rep(end) = rep(1);

xep = chord * xep;
rep = chord * rep;

if ( flipaf )
    rep = -1 * rep;
    rep = flip(rep);
    xep = flip(xep);
end

% 'Angle of attack'
A = [cos(-alpha * pi/180) -sin(-alpha * pi/180);
    sin(-alpha * pi/180)  cos(-alpha * pi/180)] * [xep rep]';
xep = A(1,:);
rep = A(2,:);

% Translate LE point.
xep = xep + xoff;
rep = rep + roff;

% Disable local curvature term
rad = inf( 1, naf - 1 );

% Vex = nan( 1, length( xep ) - 1 );
Vex = sqrt ( 1 - cp_vkt' );

end
