clear all
format compact

nx = 20;
nr = 33;

xs = linspace( -5, 5, nx );
rs = linspace( 0, 4, nr );
%rs(1) = 0.001;

x0 = 0;
r0 = 1;

u = nan( nr, nx );
v = u;

[x,r] = meshgrid( xs, rs );

for i=1:nx
    for j=1:nr
        [ u(j,i), v(j,i) ] = tubevortex( x0, r0, xs(i), rs(j) );
    end
end

figure(1)
quiver( x, r, u, v );
axis equal
