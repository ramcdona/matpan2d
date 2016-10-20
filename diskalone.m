clear all
format compact
close all

tic

% Freestream velocity
W = 1.0;

% Thrust coefficient of disk alone -- non-dimensionalized like drag
% coefficient, normalized on disk area.

CT = 5;

% Estimated cruise CT
% CT = 177.889/204.99;

% Calculate vortex tube strength and eventual jet velocity
gammainf = - W * ( sqrt( CT + 1 ) - 1 );
Wjinf = W - gammainf;

% Disk radius
rdisk = 1.0;

% Disk location
xstart = 0;
% End of contracting streamtube
xend = 6;

% Number of streamtube contraction iterations
ntstep = 10;

% Vortex ring spacing
dx = rdisk * 0.15;

% Number of vortex ring panels
npan = ( xend - xstart ) / dx;

% Set up initial geometry
xpts = linspace( xstart, xend, npan + 1 );
rpts = rdisk * ones( size( xpts ) );

xt = xpts(end);
rt = rpts(end);
gammat = gammainf;

dx = xpts(2:end) - xpts(1:end-1);
dr = rpts(2:end) - rpts(1:end-1);
ds = sqrt( dx.^2 + dr.^2 );

xcp = ( xpts(2:end) + xpts(1:end-1) ) * 0.5;
rcp = ( rpts(2:end) + rpts(1:end-1) ) * 0.5;
gammav = gammainf * ones( size( xcp ) );


disp('Relaxing streamtube')
toc

opts = odeset;
opts = odeset( opts, 'Events', @massevt );
opts = odeset( opts, 'AbsTol', 1e-5, 'RelTol', 1e-2);

for j=1:ntstep

    % Updating streamtube geometry parameters
    xcp = ( xpts(2:end) + xpts(1:end-1) ) * 0.5;
    rcp = ( rpts(2:end) + rpts(1:end-1) ) * 0.5;

    dx = xpts(2:end) - xpts(1:end-1);
    dr = rpts(2:end) - rpts(1:end-1);
    ds = sqrt( dx.^2 + dr.^2 );

    % Integrate mass flow at disk
    x0 = 0;
    mass0 = 20;
    [rdist,mdist]=ode45( @dm, [0, 1], 0, opts, W, xcp, rcp, (gammav .* ds), xt, rt, gammat, x0, mass0 );
    mass0 = mdist(end);

    % Set vortex tube radius based on jet velocity and mass flow
    rt = sqrt( mass0 / ( pi * Wjinf ) );

    % Find streamtube by continuity
    rnew = rpts;
    for i = 2:length(xpts)-2
        x0 = xpts(i);
        [rdist, mdist, re]=ode45( @dm, [0, 2], 0, opts, W, xcp, rcp, (gammav .* ds), xt, rt, gammat, x0, mass0 );
        rnew(i) = re;
    end
    % Force last two points to equal tube radius
    rnew(end-1:end) = rt;

    % Evaluate velocities at streamtube vortex locations
    % Used to adjust vortex strength along streamtube
    upts = W * ones( size(xpts) );
    vpts = zeros( size(xpts) );

    for i = 1:length(xcp)
        [ uj, vj ] = ringvortex( xcp(i), rcp(i), xpts, rpts );
        upts = upts + uj .* gammav(i) .* ds(i);
        vpts = vpts + vj .* gammav(i) .* ds(i);
    end

    for i = 1:length(xt)
        [ ut, vt ] = tubevortex( xt(i), rt(i), xpts, rpts );
        upts = upts + ut .* gammat(i);
        vpts = vpts + vt .* gammat(i);
    end

    ucp = ( upts(2:end) + upts(1:end-1) ) * 0.5;
    vcp = ( vpts(2:end) + vpts(1:end-1) ) * 0.5;

    Vpts = sqrt( upts.^2 + vpts.^2 );
    Vcp = sqrt( ucp.^2 + vcp.^2 );

    % Calculate new strength via: (vs*gamma)=const
    gammanew = ( gammainf * ( W - gammainf / 2 ) ) ./ Vcp;
    % Force last two to equal tube strength
    gammanew(end-1:end) = gammainf;

    % Track errors
    err = max(abs(rpts-rnew));
    errhist(j) = err;
    errg = max(abs(gammav-gammanew));
    errghist(j) = errg;

    % Apply updated values
    rpts = rnew;
    gammav = gammanew;
end

disp('Post-processing')
toc

figure(1)
semilogy(errhist,'x-')
hold on
semilogy(errghist,'o-')
hold off
xlabel('Streamtube and gamma iteration')
ylabel('Maximum displacement')


figure(2)
plot( xcp, -gammav, 'x-', [xt xt+2], -gammainf*[1 1],[1 1],[0 2] );
ax = axis;
ax(3) = 0; % Force y-axis to zero.
axis(ax);


if(true)

% Set up field survey points.
% Three zones -- Before disk, contracting streamtube, semi-inf streamtube
% Place in a vector for later triangulation -- does not require perfect
% grid structure.

% Before disk
xgrd = unique( [-5:.1:0-0.1] ); % xend+0.1:.1:xend+4] );
rgrd = unique( [0:.1:2 ] );
[xg, rg] = meshgrid( xgrd, rgrd );
xv = reshape( xg, 1, [] );
rv = reshape( rg, 1, [] );

% Constracting streamtube
kscale = unique( [0:.1:2 .95 .975 0.99 1.01 1.025 1.05 ] );
for i=1:length(kscale)
    xv = [ xv xpts ];
    rv = [ rv kscale(i) * rpts ];
end

% Semi-infinite streamtube
rgrdend = rt * kscale;
[xgend, rgend] = meshgrid( xend+0.1:.1:xend+4, rgrdend );
xge = reshape( xgend, 1, [] );
rge = reshape( rgend, 1, [] );
xv = [ xv xge ];
rv = [ rv rge ];

% Petrurb points randomly to improve stability of unstructured streamline
% tracing algorithm
rv = abs( rv + rand(size(rv))*.001-.0005 );

% Survey velocity
uv = W * ones( size(xv) );
vv = zeros( size(xv) );

for i = 1:length(xcp)
    [ uj, vj ] = ringvortex( xcp(i), rcp(i), xv, rv );
    uv = uv + uj .* gammav(i) .* ds(i);
    vv = vv + vj .* gammav(i) .* ds(i);
end

for i = 1:length(xt)
    [ ut, vt ] = tubevortex( xt(i), rt(i), xv, rv );
    uv = uv + ut .* gammat(i);
    vv = vv + vt .* gammat(i);
end

Vmagv = sqrt(uv.^2+vv.^2);
Cpv = 1 - Vmagv.^2;

% Triangulate survey points
tri = delaunay( xv', rv' );

% Streamline starting points
rsl = 0:.1:3 + 0.05;
xsl = -1 * ones(size(rsl)) + 0.01;

% Settings to unstructured streamline integrator
verbose = 0;  % flag to report progress
maxits = 1e5; % maximum number of iterations
Ltol = 0.01;  % flowpath "out-of-triangle" tolerance
dLtol = 0.5;  % flowpath "curvature" tolerance

% Integrate streamlines
FlowP = TriStream( tri, xv, rv, uv, vv, xsl, rsl, verbose, maxits, Ltol, dLtol );

% Generate plots
figure(4)
quiver( xv, rv, uv, vv )
hold on
plot( xpts, rpts );
plot([0 0],[0 1],'k')
hold off
axis equal


figure(5)
trisurf( tri, xv, rv, Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
hold on
plot( xpts, rpts );
plot([0 0],[0 1],'k')
hold off
axis equal
view(0,90)
title('v/Vinf')


figure(6)
trisurf( tri, xv, rv, Cpv, 'LineStyle', 'none', 'FaceColor', 'interp' );
hold on
plot( xpts, rpts );
plot([0 0],[0 1],'k')
hold off
axis equal
view(0,90)
title('Cp')


figure(7)
PlotTriStream( FlowP );
hold on
plot( xpts, rpts );
plot([0 0],[0 1],'k')
hold off
axis equal

figure(8)
trisurf( tri, xv, rv, zeros(size(xv)), Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
hold on
plot( xpts, rpts );
plot([0 0],[0 1],'k')
hold off
hold on
PlotTriStream( FlowP, 'k' );
axis equal
view(0,90)
title('v/Vinf')


mask = (rv<.01);
figure(9)
plot(xv(mask),Vmagv(mask),'o',[xgrd(1) xend+4],(W-gammainf)*[1 1],'-' )
ax=axis;
ax(3)=0;
axis(ax);
ylabel('Center-line Velocity')

disp('Done')
toc

end
