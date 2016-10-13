clear all
format compact

W = 1.0;
gamma = 1.5;

xgrd = unique( [-1:.1:1 -.1:.025:.1] ) + (pi-3)/1000;

rgrd = unique( [0:.1:2 0.8:.025:1.1] ) + (pi-3)/1000;

rsl = 0:.1:3 + 0.05;
xsl = -1 * ones(size(rsl)) + 0.01;


nx = length( xgrd );
nr = length( rgrd );

[xg, rg] = meshgrid( xgrd, rgrd );

xv = reshape( xg, 1, [] );
rv = reshape( rg, 1, [] );

xv = xv + rand(size(xv))*.01 - .005;
rv = rv + rand(size(rv))*.01 - .005;

rv=abs(rv);

uv = W * ones( size(xv) );
vv = zeros( size(xv) );

[ uj, vj ] = tubevortex( 0, 1, xv, rv );
uv = uv + uj * gamma;
vv = vv + vj * gamma;



tri = delaunay( xv', rv' );


verbose = 0;  % flag to report progress
maxits = 1e4; % maximum number of iterations
Ltol = 0.01;  % flowpath "out-of-triangle" tolerance
dLtol = 0.5;  % flowpath "curvature" tolerance

streamtubedn = TriStream( tri, xv, rv, uv, vv, 0.0, 1.0, verbose, maxits, Ltol, dLtol );
streamtube = TriStream( tri, xv, rv, -uv, -vv, 0.0, 1.0, verbose, maxits, Ltol, dLtol );

streamtube.x = [flip(streamtube.x) streamtubedn.x];
streamtube.y = [flip(streamtube.y) streamtubedn.y];
streamtube.t = [flip(streamtube.t) streamtubedn.t];
streamtube.done = [flip(streamtube.done) streamtubedn.done];


PlotTriStream( streamtube );



Vmagv = sqrt(uv.^2+vv.^2);
Cpv = 1 - Vmagv.^2;

figure(4)
quiver( xv, rv, uv, vv )
hold on
plot([0 0],[0 1],'k')
PlotTriStream( streamtube );
hold off
axis equal


figure(5)
trisurf( tri, xv, rv, Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
hold on
plot([0 0],[0 1],'k')
PlotTriStream( streamtube );
hold off
axis equal
view(0,90)
title('v/Vinf')


figure(6)
trisurf( tri, xv, rv, Cpv, 'LineStyle', 'none', 'FaceColor', 'interp' );
hold on
plot([0 0],[0 1],'k')
PlotTriStream( streamtube );
hold off
axis equal
view(0,90)
title('Cp')

figure(7)
FlowP = TriStream( tri, xv, rv, uv, vv, xsl, rsl, verbose, maxits, Ltol, dLtol );
PlotTriStream( FlowP );
hold on
plot([0 0],[0 1],'k')
PlotTriStream( streamtube );
hold off
axis equal

figure(8)
trisurf( tri, xv, rv, zeros(size(xv)), Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
hold on
plot([0 0],[0 1],'k')
PlotTriStream( streamtube );
hold off
hold on
PlotTriStream( FlowP, 'k' );
axis equal
view(0,90)
title('v/Vinf')