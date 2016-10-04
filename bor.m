clear all
format compact

% Body of revolution surface vorticity code from
% Vortex Element Methods for Fluid Dynamic Analysis of Engineering Systems
% R. I. Lewis, 1991, Cambridge University Press, Ch. 4.

% Freestream velocity
W = 1.0;

drawplots = true;

% 1 -- Sphere
% 2 -- BOR with aft cone
% 3 -- Ellipsoid
runcase = 3;

if ( runcase == 1 )
    % Construct a circle
    ncirc = 36;
    rad = 0.5;
    thetacirc = linspace( pi, 0, ncirc );

    % Panel endpoints
    xep = rad * cos( thetacirc );
    rep = rad * sin( thetacirc );

    % % Exact solution
    Vmag = 1.5 * sin( thetacirc );
    %
    % % [potex, uex, vex, wex] = sphereflow(xep, rep, zeros(size(xep)), [W 0 0], rad );
    % % Vmag = sqrt( uex.^2 + vex.^2 + wex.^2 );
    %
    % Cpex = 1 - Vmag.^2;

elseif (runcase == 2)

    % Construct example body
    ncirc = 10;
    nstraight = 22;
    rad = 0.16;
    thetacirc = linspace( pi, pi/2, ncirc );

    xep1 = rad * cos( thetacirc ) + rad;
    rep1 = rad * sin( thetacirc );
    rad1 = rad * ones( 1, ncirc - 1 );

    xep2 = linspace( rad, 0.77, nstraight );
    rep2 = rad * ones( 1, nstraight );

    xep3 = linspace( 0.77, 1.367128, nstraight);
    rep3 = linspace( rad, 0, nstraight);

    xep = [xep1 xep2(2:end) xep3(2:end)];
    rep = [rep1 rep2(2:end) rep3(2:end)];
else

    % Construct an ellipse
    ncirc = 101;
    rada = 5.0;
    radb = 0.5;

    thetacirc = linspace( pi, 0, ncirc );

    % Panel endpoints
    xep = rada * cos( thetacirc );
    rep = radb * sin( thetacirc );

    [pot, u, v, w] = ellipsoidflow( xep, rep, zeros(1,ncirc), [W 0 0], rada, radb, radb );
    Vmag = sqrt(u.^2 + v.^2 + w.^2);
end

npt = length( xep );
npan = npt - 1;

if ( drawplots)
    figure( 1 )
    plot( xep, rep, 'x-' )
    axis equal
end

% Process geometry
dx = xep(2:end) - xep(1:end-1);
dr = rep(2:end) - rep(1:end-1);

ds = sqrt( dx.^2 + dr.^2 );                    % Panel arclength
theta = atan2( dr, dx );                       % Panel slope angle
xcp = 0.5 * ( xep(2:end) + xep(1:end-1) );     % Panel x center point
rcp = 0.5 * ( rep(2:end) + rep(1:end-1) );     % Panel r center point

% Build up curvature term at control points.
if ( runcase == 2 )
    rad = [rad1 inf(1, npan - length(rad1) )];
elseif ( runcase == 3 )
    rad = rada.^2*radb.^2.*( xcp.^2./rada.^4 + rcp.^2./radb.^4 ).^1.5;
    % rad = inf(1, ncirc - 1 );
end

AIC = nan( npan, npan );
% i loop over panels implied by . operations
for j = 1:npan  % Loop over vortex j
    [ uj, vj ] = ringvortex( xcp(j), rcp(j), xcp, rcp );
    AIC(:,j) = ( uj .* cos( theta ) + vj .* sin( theta ) ) .* ds(j);
end

C = 4.0 * pi * rcp ./ ds;
curve = ds ./ ( 4.0 * pi * rad );  % Constant radius of curvature for sphere

% Self influence term
AICii = -0.5 - ( log( 2.0 * C ) - 0.25 ) ./ C .* cos( theta ) + curve;
% Excessively clever way to substitute the matrix diagonal.
AIC( 1 : npan + 1 : npan * npan ) = AICii;

% Setup RHS
rhs = -W * cos( theta );

% Solve (gamma is tangental velocity at control points)
gamma = AIC \ rhs';


Cp = 1 - gamma.^2;

if( drawplots )

% Post-process
figure(2)
plot( xcp, gamma, 'o', xep, rep )
ylabel('v/Vinf')

figure(3)
plot( xcp, -Cp );
ylabel('-C_p')

% Do some case-specific post-processing.
if ( runcase == 1 )
    figure(2)
    hold on
    plot( xep, Vmag )
    hold off

    % Exact solution at control points.
    Vex = 1.5 * sin( atan2( rcp, xcp ) );

    err = Vex - gamma';
    disp(['Maximum velocity error:  ' num2str(max(abs(err))) ]);
elseif ( runcase == 3 )
    figure(2)
    hold on
    plot( xep, Vmag )
    hold off


    [pot, u, v, w] = ellipsoidflow( xcp, rcp, zeros(1,npan), [W 0 0], rada, radb, radb );
    Vex = sqrt(u.^2 + v.^2 + w.^2);

    err = Vex - gamma';
    disp(['Maximum velocity error:  ' num2str(max(abs(err))) ]);

end

end
