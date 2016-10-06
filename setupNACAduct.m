function [xep, rep, rad, Vex] = setupNACAduct( naf, chord, alpha, xoff, roff, dig1, dig2, dig34, flipaf, spacing )

[ xep, rep ] = naca4( dig1, dig2, dig34, naf, spacing );

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

Vex = nan( 1, length( xep ) - 1 );

end
