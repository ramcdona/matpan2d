% Body of revolution surface vorticity code from
% Vortex Element Methods for Fluid Dynamic Analysis of Engineering Systems
% R. I. Lewis, 1991, Cambridge University Press, Ch. 4.

% This solver script expects some variables to be set as input before
% execution.
%
% There may be multiple input bodies/objects, each an item in a cell array.
%
% W                % Freestream velocity
% Sref             % Reference area for calculating CD
% drawplots        % Flag to draw plots and postprocess
% xgrd             % X Grid points for velocity survey
% rgrd             % R Grid points for velocity survey
% xsl              % X Points to start streamlines
% rsl              % R Points to start streamlines
% streamback       % Flag to trace streamlines backwards
% ntstep           % Number of streamtube relaxation steps
%
% names{nseg};     % Component names for output
% xepts{nseg};     % X Endpoints of panels (or initial streamtube)
% repts{nseg};     % R Endpoints of panels (or initial streamtube)
% remin{nseg};     % R Minimum for relaxed streamtube
% remax{nseg};     % R Maximum for relaxed streamtube
% kuttas{nseg};    % Flags to apply Kutta condition to each body
% props{nseg};     % Flag to treat body as actuator disk
% gammaads{nseg};  % Initial actuator disk strength vector
% Wjad{nseg};      % Far downstream disk jet velocity
% jtels{nseg};     % Index to TE lower panel
% jteus{nseg};     % Index to TE upper panel
% rads{nseg};      % Panel radius of curvature vector
% Vexs{nseg};      % Exact solution at panel endpoints for isolated body
%

nseg = length( xepts );

if ( drawplots)
    figure( 1 )
    clf
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg}, 'x-' );
    end
    hold off
    axis equal
end

npan = 0;
% Process geometry
for iseg=1:nseg
    xep = xepts{iseg};
    rep = repts{iseg};

    dx{iseg} = xep(2:end) - xep(1:end-1);
    dr{iseg} = rep(2:end) - rep(1:end-1);

    ds{iseg} = sqrt( dx{iseg}.^2 + dr{iseg}.^2 );        % Panel arclength
    theta{iseg} = atan2( dr{iseg}, dx{iseg} );           % Panel slope angle
    xcp{iseg} = 0.5 * ( xep(2:end) + xep(1:end-1) );     % Panel x center point
    rcp{iseg} = 0.5 * ( rep(2:end) + rep(1:end-1) );     % Panel r center point

    if ( ~props{iseg} ) % 'Normal' components
        npans{iseg} = length( xcp{iseg} );
        npan = npan + npans{iseg};
    else % Actuator disk doesn't contribute panels
        npans{iseg} = 0;
    end
end


AICs = cell( nseg, nseg );
for iseg=1:nseg
    for jseg = 1:nseg

        AICij = nan( npans{jseg}, npans{iseg} );

        if ( ~props{iseg} && ~props{jseg} )

            % i loop over panels implied by . operations
            for j = 1:npans{iseg}  % Loop over vortex j
                [ uj, vj ] = ringvortex( xcp{iseg}(j), rcp{iseg}(j), xcp{jseg}, rcp{jseg} );
                AICij(:,j) = ( uj .* cos( theta{jseg} ) + vj .* sin( theta{jseg} ) ) .* ds{iseg}(j);
            end

            if ( iseg == jseg )
                C = 4.0 * pi * rcp{iseg} ./ ds{iseg};
                curve = ds{iseg} ./ ( 4.0 * pi * rads{iseg} );

                % Self influence term
                AICii = -0.5 - ( log( 2.0 * C ) - 0.25 ) ./ C .* cos( theta{iseg} ) + curve;
                % Excessively clever way to substitute the matrix diagonal.
                AICij( 1 : npans{iseg} + 1 : npans{iseg} * npans{iseg} ) = AICii;

                if ( kuttas{iseg} )
                    % Correction term for back-diagonal
                    for j = 1:npans{iseg} % Loop over column j
                        AICij(j,j) = - ( sum( AICij(:,j) .* ds{iseg}' ) - AICij(j,j) * ds{iseg}(j) ) / ds{iseg}(j);
                    end
                end
            end
        end

        AICs{iseg, jseg} = AICij;
    end
end

% Assemble AIC from cell matrices.
AIC=[];
panidx=[];
segidx=[];
for iseg=1:nseg
    AICrow = [];
    for jseg = 1:nseg
        AICrow = [AICrow AICs{jseg, iseg}];
    end
    AIC = [AIC; AICrow];

    panidx = [panidx 1:npans{iseg}];
    segidx = [segidx iseg*ones(1,npans{iseg})];
end
idx = [segidx; panidx; 1:npan];
modidx = idx;

% Apply Kutta condition to AIC
for iseg=1:nseg
    if ( kuttas{iseg} )
        jtelow = find( ( idx(1,:) == iseg ) & ( idx(2,:) == jtels{iseg} ) );
        jteup = find( ( idx(1,:) == iseg ) & ( idx(2,:) == jteus{iseg} ) );

        % Apply Kutta condition to AIC
        AIC(:,jtelow) = AIC(:,jtelow) - AIC(:,jteup);
        AIC(:,jteup)=[];
        AIC(jtelow,:) = AIC(jtelow,:) - AIC(jteup,:);
        AIC(jteup,:)=[];

        % Remove index rows to match
        modidx(:,jteup)=[];
    end
end

opts = odeset;
opts = odeset( opts, 'Events', @bormassevt );
opts = odeset( opts, 'AbsTol', 1e-5, 'RelTol', 1e-2);

for itstep=1:ntstep

    rerr = 0;

    % Re-process actuator disk geometry
    for iseg=1:nseg
        if( props{iseg} )
            xep = xepts{iseg};
            rep = repts{iseg};
            rmin = remin{iseg};
            rmax = remax{iseg};

            % Limit streamtube radius
            rep( rep < rmin ) = rmin( rep < rmin );
            rep( rep > rmax ) = rmax( rep > rmax );

            % Store limited streamtube
            repts{iseg} = rep;

            dx{iseg} = xep(2:end) - xep(1:end-1);
            dr{iseg} = rep(2:end) - rep(1:end-1);

            ds{iseg} = sqrt( dx{iseg}.^2 + dr{iseg}.^2 );        % Panel arclength
            theta{iseg} = atan2( dr{iseg}, dx{iseg} );           % Panel slope angle
            xcp{iseg} = 0.5 * ( xep(2:end) + xep(1:end-1) );     % Panel x center point
            rcp{iseg} = 0.5 * ( rep(2:end) + rep(1:end-1) );     % Panel r center point
        end
    end

    % Setup & Assemble RHS
    rhs = [];
    for iseg=1:nseg
        if ( ~props{iseg} ) % 'Normal' components
            rhss{iseg} = -W * cos( theta{iseg} );

            for jseg=1:nseg
                if ( props{jseg} )

                    uad = zeros( size( xcp{iseg} ) );
                    vad = uad;
                    for j = 1:length( xcp{jseg} )  % Loop over vortex j
                        [ uj, vj ] = ringvortex( xcp{jseg}(j), rcp{jseg}(j), xcp{iseg}, rcp{iseg} );
                        uad = uad + uj * gammaad{jseg}(j) * ds{jseg}(j);
                        vad = vad + vj * gammaad{jseg}(j) * ds{jseg}(j);
                    end

                    [ uj, vj ] = tubevortex( xepts{jseg}(end), repts{jseg}(end), xcp{iseg}, rcp{iseg} );
                    uad = uad + gammaad{jseg}(end) * uj;
                    vad = vad + gammaad{jseg}(end) * vj;
                    rhss{iseg} = rhss{iseg} - ( uad .* cos( theta{iseg} ) + vad .* sin( theta{iseg} ) );
                end
            end

            rhs = [rhs  rhss{iseg}];
        end
    end

    % Apply Kutta condition to RHS
    for iseg=1:nseg
        if ( kuttas{iseg} )
            jtelow = find( ( idx(1,:) == iseg ) & ( idx(2,:) == jtels{iseg} ) );
            jteup = find( ( idx(1,:) == iseg ) & ( idx(2,:) == jteus{iseg} ) );

            % Apply Kutta condition to RHS
            rhs(:,jtelow) = rhs(:,jtelow) - rhs(:,jteup);
            rhs(:,jteup)=[];
        end
    end

    % Solve (gamma is tangental velocity at control points)
    gamma = AIC \ rhs';

    % Break gamma vector into per-body parts
    for iseg=1:nseg
        if( ~props{iseg} )
            gam = gamma( modidx(1,:) == iseg );

            % Put kutta condition back on for plotting and post-processing.
            if ( kuttas{iseg} )
                jtelow = jtels{iseg};
                jteup = jteus{iseg};

                gam = [gam(1:jteup-1); -gam(jtelow); gam(jteup+1:end)];
            end

            % Flip reversed panel velocities for plotting
            gam = gam .* sign(dx{iseg})';

            gammas{iseg} = gam';
            Cp{iseg} = 1 - gam.^2;

        else
            gammas{iseg} = gammaad{iseg};
        end
    end


    % Update streamtube shape
    for iseg=1:nseg
        if( props{iseg} )
            % Integrate mass flow at disk
            x0 = xepts{iseg}(1);
            mass0 = inf;
            [ rdist, mdist ] = ode45( @bordm, [0, repts{iseg}(1)], 0, opts, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props, x0, mass0 );
            mass0 = mdist(end);

            % Set vortex tube radius based on jet velocity and mass flow
            rt = sqrt( mass0 / ( pi * Wjad{iseg} ) );

            % Find streamtube by continuity
            rnew = repts{iseg};
            for i = 2:length(xepts{iseg})-2
                x0 = xpts(i);
                [rdist, mdist, re] = ode45( @bordm, [0, remax{iseg}(i)], 0, opts, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props, x0, mass0 );
                rnew(i) = re;
            end
            % Force last two points to equal tube radius
            rnew(end-1:end) = rt;

            rerr = max( rerr, max( abs( repts{iseg} - rnew ) ) );
            repts{iseg} = rnew;
        end
    end

end

if( drawplots )

    % Post-process
    figure(2)
    clf
    hold on
    for iseg=1:nseg
        if( ~props{iseg} )
            plot( xcp{iseg}, gammas{iseg}, 'o-' )
            plot( xcp{iseg}, Vexs{iseg} )
        end
    end
    hold off
    ylabel('v/Vinf')

    figure(3)
    clf
    hold on
    for iseg=1:nseg
        if( ~props{iseg} )
            plot( xcp{iseg}, -Cp{iseg} );
        end
    end
    hold off
    ylabel('-C_p')

end

% Calculate maximum error for cases with exact isolated solution
for iseg=1:nseg
    if ( ~isnan( Vexs{iseg}(1) ) )
        err = Vexs{iseg} - gammas{iseg};
        disp(['Maximum velocity error:  ' num2str(max(abs(err))) ]);
    end
end

% Post-process axial force.
CDi = [];
CD = 0.0;
for iseg=1:nseg
    if( ~props{iseg} )

        CDi{iseg} = - 2.0 * pi * sum( Cp{iseg}' .* sin( theta{iseg} ) .* rcp{iseg} .* ds{iseg} ) / Sref;
        CD = CD + CDi{iseg};

        disp(['CD on ' names{iseg} ':  ' num2str(CDi{iseg})])
    end
end

disp(['CD total:  ' num2str(CD)])

% Post-process velocity survey.
if( drawplots )
    nx = length( xgrd );
    nr = length( rgrd );

    [xg, rg] = meshgrid( xgrd, rgrd );

    xv = reshape( xg, 1, [] );
    rv = reshape( rg, 1, [] );

    inv = zeros( size(xv) );

    for iseg=1:nseg
        if( ~props{iseg} )
            xep = xepts{iseg};
            rep = repts{iseg};

            % Forces open polys to be closed
            % Signed minimum distance -- negative is inside.
            [dmin, x_d_min, y_d_min, is_vertex, idx_c] = p_poly_dist( xv, rv, xep, rep, true );

            % Loop to compare to nearby edge length.  Velocity survey has
            % problems within one edge length of body
            for i=1:length( inv )
                if ( idx_c(i) < length( ds{iseg} ) && idx_c(i) > 0 )
                    inv(i) = inv(i) | ( dmin(i) < ds{iseg}(idx_c(i) ) );
                end
            end
        end
    end

    [uv, vv] = borvel( xv, rv, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props );

    ug = reshape( uv, size(xg) );
    vg = reshape( vv, size(xg) );

    xv = xv( ~inv );
    rv = rv( ~inv );
    uv = uv( ~inv );
    vv = vv( ~inv );

    Cedg = [];

    for iseg=1:nseg
        if( ~props{iseg} )
            nnext = length(xv) + 1;
            mcon = length(xcp{iseg}) - 1;

            % Body on center line
            if ( abs(repts{iseg}(1)) < 1e-6 )
                xv = [xv xepts{iseg}(1)];
                rv = [rv repts{iseg}(1)];
                uv = [uv 0];
                vv = [vv 0];
                mcon = mcon + 1;
            end

            xv = [xv xcp{iseg}];
            rv = [rv rcp{iseg}];

            uv = [uv gammas{iseg} .* cos( theta{iseg} ) .* sign(dx{iseg}) ];
            vv = [vv gammas{iseg} .* sin( theta{iseg} ) .* sign(dx{iseg}) ];

            % Body on center line
            if ( abs(repts{iseg}(end)) < 1e-6 )
                xv = [xv xepts{iseg}(end)];
                rv = [rv repts{iseg}(end)];
                uv = [uv 0];
                vv = [vv 0];
                mcon = mcon + 1;
            end

            Cedg = [Cedg [ nnext:nnext+mcon; nnext+1:nnext+mcon nnext ] ];
        end
    end

    if( isempty(Cedg) )
        DT = delaunayTriangulation( xv', rv' );
        tri = DT.ConnectivityList;
        IO = true( size( tri, 1 ), 1 );
    else
        DT = delaunayTriangulation( xv', rv', Cedg' );
        tri = DT.ConnectivityList;
        IO = ~isInterior(DT);
    end

    Vmagv = sqrt(uv.^2+vv.^2);
    Cpv = 1 - Vmagv.^2;

    figure(4)
    quiver( xv, rv, uv, vv )
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg} );
    end
    hold off
    axis equal


    figure(5)
    trisurf( tri(IO,:), xv, rv, Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg} );
    end
    hold off
    axis equal
    view(0,90)
    title('v/Vinf')


    figure(6)
    trisurf( tri(IO,:), xv, rv, Cpv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg} );
    end
    hold off
    axis equal
    view(0,90)
    title('Cp')

    figure(7)
    verbose = 0;  % flag to report progress
    maxits = 1e4; % maximum number of iterations
    Ltol = 0.01;  % flowpath "out-of-triangle" tolerance
    dLtol = 0.5;  % flowpath "curvature" tolerance

    if( streamback )
        FlowP = TriStream( tri, xv, rv, -uv, -vv, xsl, rsl, verbose, maxits, Ltol, dLtol );
    else
        FlowP = TriStream( tri, xv, rv, uv, vv, xsl, rsl, verbose, maxits, Ltol, dLtol );
    end
    PlotTriStream( FlowP );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg} );
    end
    plot( xgrd, zeros(size(xgrd)),'k');
    hold off
    axis equal

    figure(8)
    trisurf( tri(IO,:), xv, rv, zeros(size(xv)), Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
    PlotTriStream( FlowP, 'k' );
    hold on;  % PlotTriStream turns hold off.
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg},'k' );
    end
    plot( xgrd, zeros(size(xgrd)),'k');
    hold off
    axis equal
    view(0,90)
    title('v/Vinf')

end
