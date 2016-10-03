clear all
format compact

nx = 10;
nr = 11;

xs = linspace( -2, 2, nx );
rs = linspace( 0, 2, nr );

x0 = 0;
r0 = 1;

u = nan( nr, nx );
v = u;

[x,r] = meshgrid( xs, rs );

for i=1:nx
    for j=1:nr
        [ u(j,i), v(j,i) ] = ringvortex( x0, r0, xs(i), rs(j) );
    end
end

figure(1)
quiver( x, r, u, v );