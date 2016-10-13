function [ u, v ] = ringsubvortex( xm, rm, xs, rs, xe, re, varargin )
% Velocity induced by unit-strength ring vortex around the x-axis.
%
% Binary subdivision algorithm for near-field  points
%

k = 1;       % panel length multiplier for distance check
dltol = 1e-3;  % panel length to return zero.

xn = ( xe + xs ) * 0.5;
rn = ( re + rs ) * 0.5;

dl = xrdist( xe, re, xs, rs );

if ( nargin > 6 )
    dlorig = varargin{1};
else
    dlorig = dl;
end

if ( dl/dlorig < dltol )
    u = 0.0;
    v = 0.0;
    return;
end

%de = xrdist( xm, rm, xe, re );
%ds = xrdist( xm, rm, xs, rs );
dn = xrdist( xm, rm, xn, rn );
%d = min( [de ds dn] );
d = dn;

if ( d > k * dl )  % 'Far' treatment
    [u, v] = ringvortex( xn, rn, xm, rm );
else
    
    [u1, v1] = ringsubvortex( xm, rm, xs, rs, xn, rn, dlorig );
    [u2, v2] = ringsubvortex( xm, rm, xn, rn, xe, re, dlorig );
    
    u = ( u1 + u2 ) * 0.5;
    v = ( v1 + v2 ) * 0.5;
end

end

function d = xrdist( x1, r1, x2, r2 )
d = sqrt( ( x2 - x1 ).^2 + ( r2 - r1 ).^2 );
end
