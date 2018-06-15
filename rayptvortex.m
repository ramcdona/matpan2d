function [ u, v ] = rayptvortex( xs, ys, xv, yv, x, y )
% Velocity induced by unit-strength semi-infinite constrant strength
% vortex panel approximated with many point vortices.
%
% xs, ys Ray vortex start
% xv, yv Ray vortex direction vector
% x, y   Point of interest

% Assume x,y passed in can be vectors for dot notation.  Must loop over
% point vortices explicitly.

dl = sqrt(xv.^2+yv.^2);
xv = xv ./ dl;
yv = yv ./ dl;

fardist = 20;
ds = 0.1;
n = fardist / ds;

u = zeros(size(x));
v = u;
for i=1:n
    xn = xs + xv * (i - 0.5) * ds;
    yn = ys + yv * (i - 0.5) * ds;
    [ uu, vv ] = ptvortex( xn, yn, x, y );
    u = u + uu * ds;
    v = v + vv * ds;
end

end
