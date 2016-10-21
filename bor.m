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
% deltaCP{nseg};   % Actuator disk pressure jump
% jtels{nseg};     % Index to TE lower panel
% jteus{nseg};     % Index to TE upper panel
% rads{nseg};      % Panel radius of curvature vector
% Vexs{nseg};      % Exact solution at panel endpoints for isolated body
%

nseg = length( xepts );

for iseg = 1:nseg

    if ( props{iseg} )  % Only check props

        gammainf = - W * ( sqrt( deltaCP{iseg} + 1 ) - 1 );
        gammas{iseg} = gammainf * ones( size( xepts{iseg} ) );

        xpts = xepts{iseg};
        rmin = remin{iseg};
        rmax = remax{iseg};

        for jseg = 1:nseg
            if( ~props{jseg} )  % Skip other props

                xtrim = xepts{jseg};
                rtrim = repts{jseg};

                for i=1:length(xpts)

                    [x,r] = polyxpoly( xpts(i) * [1 1], [rmin(i) rmax(i)], xtrim, rtrim );

                    if( ~isempty(x) )
                        if( ~kuttas{jseg} ) % Center line body
                            rmin( i ) = max( r );
                        else  % Duct
                            rmax( i ) = min( r );
                        end
                    end
                end
            end
        end

        rep = repts{iseg};
        rep( rep < rmin ) = rmin( rep < rmin );
        rep( rep > rmax ) = rmax( rep > rmax );
        repts{iseg} = rep;

        remin{iseg} = rmin;
        remax{iseg} = rmax;

    end
end

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

% Apply Kutta condition to AIC
for iseg=1:nseg
    if ( kuttas{iseg} )
        jtelow = find( ( idx(1,:) == iseg ) & ( idx(2,:) == jtels{iseg} ) );
        jteup = find( ( idx(1,:) == iseg ) & ( idx(2,:) == jteus{iseg} ) );

        % Apply Kutta condition to AIC
        AIC(jteup,:) = zeros( size( AIC(jteup,:) ) );
        AIC(jteup,jteup) = 1;
        AIC(jteup,jtelow) = 1;
    end
end

opts = odeset;
opts = odeset( opts, 'Events', @bormassevt );
opts = odeset( opts, 'AbsTol', 1e-5, 'RelTol', 1e-2);

for itstep=1:ntstep

    rerr = 0;
    gerr = 0;

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
                        uad = uad + uj * gammas{jseg}(j) * ds{jseg}(j);
                        vad = vad + vj * gammas{jseg}(j) * ds{jseg}(j);
                    end

                    [ uj, vj ] = tubevortex( xepts{jseg}(end), repts{jseg}(end), xcp{iseg}, rcp{iseg} );
                    uad = uad + gammas{jseg}(end) * uj;
                    vad = vad + gammas{jseg}(end) * vj;
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
            rhs(:,jteup) = 0;
        end
    end

    % Solve (gamma is tangental velocity at control points)
    gamma = AIC \ rhs';

    % Break gamma vector into per-body parts
    for iseg=1:nseg
        if( ~props{iseg} )
            gam = gamma( idx(1,:) == iseg );

            % Flip reversed panel velocities for plotting
            gam = gam .* sign(dx{iseg})';

            gammas{iseg} = gam';
            Cp{iseg} = 1 - gam.^2;

            for jseg=1:nseg
                if( props{jseg} )
                    xst = xepts{jseg};
                    rst = repts{jseg};

                    xst = [xst(1) xst xst(end)+10 xst(end)+10];
                    rst = [ 0 rst rst(end) 0];

                    % Forces open polys to be closed
                    % Signed minimum distance -- negative is inside.
                    [dmin, x_d_min, y_d_min, is_vertex, idx_c] = p_poly_dist( xcp{iseg}, rcp{iseg}, xst, rst, true );

                    mask = dmin <= 0.0;

                    Cp{iseg}(mask) = Cp{iseg}(mask) + deltaCP{iseg};

                end
            end
        end
    end


    rnew = repts;
    gammanew = gammas;
    % Update streamtube shape
    for iseg=1:nseg
        if( props{iseg} )
            % Integrate mass flow at disk
            x0 = xepts{iseg}(1);
            mass0 = inf;
            [ rdist, mdist ] = ode45( @bordm, [remin{iseg}(1), repts{iseg}(1)], 0, opts, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props, x0, mass0 );
            mass0 = mdist(end);

            gammainf = - W * ( sqrt( deltaCP{iseg} + 1 ) - 1 );
            Wjad = W - gammainf;

            % Set vortex tube radius based on jet velocity and mass flow
            rt = sqrt( mass0 / ( pi * Wjad ) );

            % Find streamtube by continuity
            for i = 2:length(xepts{iseg})-2
                x0 = xepts{iseg}(i);
                [rdist, mdist, re] = ode45( @bordm, [remin{iseg}(i), remax{iseg}(i)], 0, opts, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props, x0, mass0 );
                rnew{iseg}(i) = re;
            end
            % Force last two points to equal tube radius
            rnew{iseg}(end-1:end) = rt;

            [upts, vpts] = borvel( xepts{iseg}, repts{iseg}, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props );

            ucp = ( upts(2:end) + upts(1:end-1) ) * 0.5;
            vcp = ( vpts(2:end) + vpts(1:end-1) ) * 0.5;
            Vcp = sqrt( ucp.^2 + vcp.^2 );

            % Calculate vortex tube strength and eventual jet velocity
            gammainf = - W * ( sqrt( deltaCP{iseg} + 1 ) - 1 );

            % Calculate new strength via: (vs*gamma)=const
            gammanew{iseg}(1:end-1) = ( gammainf * ( W - gammainf / 2 ) ) ./ Vcp;
            % Force last two to equal tube strength
            gammanew{iseg}(end-1:end) = gammainf;

            rerr = max( rerr, max( abs( repts{iseg} - rnew{iseg} ) ) );
            gerr = max( gerr, max( abs( gammas{iseg} - gammanew{iseg} ) ) );

        end
    end

    repts = rnew;
    gammas = gammanew;


    rerrhist(itstep) = rerr;
    gerrhist(itstep) = gerr;
end

if( drawplots )

    figure(1)
    hold on
    for iseg=1:nseg
        if( props{iseg} )
            plot( xepts{iseg}, repts{iseg}, 'o-' );
        end
    end
    hold off

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

    figure(4)
    semilogy( rerrhist, 'x-' )
    hold on
    semilogy( gerrhist, 'o-' )
    legend( 'Streamtube residual', 'Gamma residual' )
    hold off

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

    rv = abs( rv + rand(size(rv))*.001-.0005 );

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
                inv(i) = inv(i) | ( dmin(i) < 0 );
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


    for iseg=1:nseg
        if( props{iseg} )
            xst = xepts{iseg};
            rst = repts{iseg};

            xst = [xst(1) xst xst(end)+10 xst(end)+10];
            rst = [ 0 rst rst(end) 0];

            % Forces open polys to be closed
            % Signed minimum distance -- negative is inside.
            [dmin, x_d_min, y_d_min, is_vertex, idx_c] = p_poly_dist( xv, rv, xst, rst, true );

            mask = dmin <= 0.0;

            Cpv(mask) = Cpv(mask) + deltaCP{iseg};

        end
    end


    figure(5)
    quiver( xv, rv, uv, vv )
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg} );
    end
    hold off
    axis equal


    figure(6)
    trisurf( tri(IO,:), xv, rv, zeros(size(xv)), Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg} );
    end
    hold off
    axis equal
    view(0,90)
    title('v/Vinf')


    figure(7)
    trisurf( tri(IO,:), xv, rv, zeros(size(xv)), Cpv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg} );
    end
    hold off
    axis equal
    view(0,90)
    title('Cp')

    figure(8)
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

    figure(9)
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


    xsurvey = [-1 0 1 2 4];

    for i=1:length(xsurvey);
        mask = ( xv == xsurvey(i) );

        r = rv(mask);
        v = vv(mask);
        u = uv(mask);

        [r, indx] = sort( r );
        v = v(indx);
        u = u(indx);

        figure(10)
        plot( v, r )
        hold on

        figure(11)
        plot( u, r );
        hold on
    end
    figure(10)
    hold off
    figure(11)
    hold off

    figure(12)
    trimesh( tri(IO,:), xv, rv, zeros(size(xv)) );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, repts{iseg},'k','LineWidth', 2.0 );
    end
    axis equal
    view(0,90)

end
