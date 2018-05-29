% Planar surface vorticity code from
% Vortex Element Methods for Fluid Dynamic Analysis of Engineering Systems
% R. I. Lewis, 1991, Cambridge University Press, Ch. 4.

% This solver script expects some variables to be set as input before
% execution.
%
% There may be multiple input bodies/objects, each an item in a cell array.
%
% W                % Freestream velocity
% Minf             % Freestream Mach number for compressibility correction
% Sref             % Reference area for calculating CD
% drawplots        % Flag to draw plots and postprocess
% xlim             % X Grid limits for velocity survey
% ylim             % Y Grid limits for velocity survey
% nsurvey          % Number of points in velocity survey grids
% nsl              % Number of streamlines to trace
% ntstep           % Number of streamtube relaxation steps
%
% names{nseg};     % Component names for output
% xepts{nseg};     % X Endpoints of panels (or initial streamtube)
% yepts{nseg};     % Y Endpoints of panels (or initial streamtube)
% kuttas{nseg};    % Flags to apply Kutta condition to each body
% props{nseg};     % Flag to treat body as actuator disk
% deltaCP{nseg};   % Actuator disk pressure jump
% jtels{nseg};     % Index to TE lower panel
% jteus{nseg};     % Index to TE upper panel
% rads{nseg};      % Panel radius of curvature vector
% Vexs{nseg};      % Exact solution at panel endpoints for isolated body
%

gammagas = 1.4;

% Centrally set linewidth for figures that use it.
lw = 2;

nseg = length( xepts );

ydisk = cell(size(xepts));
for iseg = 1:nseg

    if ( props{iseg} )  % Only check props

        gammainf = - W * ( sqrt( deltaCP{iseg} + 1 ) - 1 );
        gammas{iseg} = gammainf * ones( size( xepts{iseg} ) );
    end
end

if ( drawplots)
    figure( 1 )
    clf
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, yepts{iseg}, 'x-' );
    end
    hold off
    axis equal
end

npan = 0;
% Process geometry
for iseg=1:nseg
    xep = xepts{iseg};
    yep = yepts{iseg};

    dx{iseg} = xep(2:end) - xep(1:end-1);
    dy{iseg} = yep(2:end) - yep(1:end-1);

    ds{iseg} = sqrt( dx{iseg}.^2 + dy{iseg}.^2 );        % Panel arclength
    theta{iseg} = atan2( dy{iseg}, dx{iseg} );           % Panel slope angle
    xcp{iseg} = 0.5 * ( xep(2:end) + xep(1:end-1) );     % Panel x center point
    ycp{iseg} = 0.5 * ( yep(2:end) + yep(1:end-1) );     % Panel y center point

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
                [ uj, vj ] = ptvortex( xcp{iseg}(j), ycp{iseg}(j), xcp{jseg}, ycp{jseg} );
                AICij(:,j) = ( uj .* cos( theta{jseg} ) + vj .* sin( theta{jseg} ) ) .* ds{iseg}(j);
            end

            if ( iseg == jseg )
                C = 4.0 * pi * ycp{iseg} ./ ds{iseg};
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

    yerr = 0;
    gerr = 0;

    % Re-process actuator disk geometry
    for iseg=1:nseg
        if( props{iseg} )
            xep = xepts{iseg};
            yep = yepts{iseg};

            dx{iseg} = xep(2:end) - xep(1:end-1);
            dy{iseg} = yep(2:end) - yep(1:end-1);

            ds{iseg} = sqrt( dx{iseg}.^2 + dy{iseg}.^2 );        % Panel arclength
            theta{iseg} = atan2( dy{iseg}, dx{iseg} );           % Panel slope angle
            xcp{iseg} = 0.5 * ( xep(2:end) + xep(1:end-1) );     % Panel x center point
            ycp{iseg} = 0.5 * ( yep(2:end) + yep(1:end-1) );     % Panel r center point
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
                        [ uj, vj ] = ptvortex( xcp{jseg}(j), ycp{jseg}(j), xcp{iseg}, ycp{iseg} );
                        uad = uad + uj * gammas{jseg}(j) * ds{jseg}(j);
                        vad = vad + vj * gammas{jseg}(j) * ds{jseg}(j);
                    end

                    [ uj, vj ] = tubevortex( xepts{jseg}(end), yepts{jseg}(end), xcp{iseg}, ycp{iseg} );
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

%             % If actuator disk matches duct TE, apply gamma jump.
%             for jseg=1:nseg
%                 if( props{jseg} )
%                     % Duct TE point
%                     xte = xepts{iseg}(1);
%                     yte = yepts{iseg}(1);
%                     % Disk leading point
%                     xd = xepts{jseg}(1);
%                     yd = yepts{jseg}(1);
% 
%                     d = sqrt( ( xte - xd )^2 + ( yte - yd )^2 );
% 
%                     if ( d < 1e-3 )
%                         rhs(:,jteup) = rhs(:,jteup) + gammas{jseg}(1);
%                     end
%                 end
%             end
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
        end
    end

    xnew = xepts;
    ynew = yepts;
    gammanew = gammas;
    % Update streamtube shape
    for iseg=1:nseg
        if( props{iseg} )
            % Integrate mass flow at disk
            x0 = xepts{iseg}(1);
            y0 = yepts{iseg}(1);
            xv = xepts{iseg}(end/2 + 1) - xepts{iseg}(1);
            yv = yepts{iseg}(end/2 + 1) - yepts{iseg}(1);
            dlen = sqrt( xv^2 + yv^2 );
            xv = xv / dlen;
            yv = yv / dlen;
            
            mass0 = inf;
            
            % dm = planardm( eta, mass, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props, x0, y0, xv, yv, len, mass0 )
            [ ydist, mdist ] = ode45( @planardm, [0, 1], 0, opts, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props, x0, y0, xv, yv, dlen, mass0 );
            mass0 = mdist(end);

            gammainf = - W * ( sqrt( deltaCP{iseg} + 1 ) - 1 );
            Wjad = W - gammainf;

            % Set vortex tube radius based on jet velocity and mass flow
            rt = sqrt( mass0 / ( pi * Wjad ) );

            % Find streamtube by continuity
            for i = 2:length(xepts{iseg})-2
                x0 = xepts{iseg}(i);
                [ydist, mdist, ye] = ode45( @bordm, [yemin{iseg}(i), yemax{iseg}(i)], 0, opts, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props, x0, mass0 );
                ynew{iseg}(i) = ye;
            end
            % Force last two points to equal tube radius
            ynew{iseg}(end-1:end) = rt;

            [upts, vpts] = borvel( xepts{iseg}, yepts{iseg}, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props );

            ucp = ( upts(2:end) + upts(1:end-1) ) * 0.5;
            vcp = ( vpts(2:end) + vpts(1:end-1) ) * 0.5;
            Vcp = sqrt( ucp.^2 + vcp.^2 );

            % Calculate vortex tube strength and eventual jet velocity
            gammainf = - W * ( sqrt( deltaCP{iseg} + 1 ) - 1 );

            % Calculate new strength via: (vs*gamma)=const
            gammanew{iseg}(1:end-1) = ( gammainf * ( W - gammainf / 2 ) ) ./ Vcp;
            % Force last two to equal tube strength
            gammanew{iseg}(end-1:end) = gammainf;

            yerr = max( yerr, max( abs( yepts{iseg} - ynew{iseg} ) ) );
            gerr = max( gerr, max( abs( gammas{iseg} - gammanew{iseg} ) ) );

        end
    end

    yepts = ynew;
    gammas = gammanew;


    yerrhist(itstep) = yerr;
    gerrhist(itstep) = gerr;
end



for iseg=1:nseg
    if( ~props{iseg} )
        Cp{iseg} = 1 - gammas{iseg}.^2;

        for jseg=1:nseg
            if( props{jseg} )
                xst = xepts{jseg};
                yst = yepts{jseg};

                xst = [xst xst(end)+10 xst(end)+10];
                yst = [yst yst(end) 0];

                if ( xdisk{jseg} ~= xst(1) )
                    xst = [xdisk{jseg} xdisk{jseg} xst];
                    yst = [ 0  ydisk{jseg} yst];
                else
                    xst = [xst(1) xst];
                    yst = [ 0 yst];
                end

                % Forces open polys to be closed
                % Signed minimum distance -- negative is inside.
                [dmin, x_d_min, y_d_min, is_vertex, idx_c] = p_poly_dist( xcp{iseg}, ycp{iseg}, xst, yst, true );

                mask = dmin <= 0.0;

                Cp{iseg}(mask) = Cp{iseg}(mask) + deltaCP{jseg};
            end
        end
    end
end


if( drawplots )

    figure(1)
    hold on
    for iseg=1:nseg
        if( props{iseg} )
            plot( xepts{iseg}, yepts{iseg}, 'o-' );
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
    semilogy( yerrhist, 'x-' )
    hold on
    semilogy( gerrhist, 'o-' )
    legend( 'Streamtube residual', 'Gamma residual' )
    hold off

end

drawnow

% Post-process velocity survey.
if( drawplots )
    xv = [];
    yv = [];

    xomax = max( xlim );
    xomin = min( xlim );
    yomax = max( ylim );
    yomin = min( ylim );

    % Build component near-field survey points
    for iseg=1:nseg
        xmin = min( xepts{iseg} );
        xmax = max( xepts{iseg} );

        ymin = min( yepts{iseg} );
        ymax = max( yepts{iseg} );

        xomax = max( xomax, xmax );
        xomin = min( xomin, xmin );
        yomax = max( yomax, ymax );
        yomin = min( yomin, ymin );

        if( ~props{iseg} )
            off = max( xmax - xmin, ymax - ymin ) * 0.2;

            xgrd = linspace( xmin - off, xmax + off, nsurvey );

            ygrd = linspace( ymin - off, ymax + off, nsurvey );

            [xg, yg] = meshgrid( xgrd, ygrd );

            xvi = reshape( xg, 1, [] );
            yvi = reshape( yg, 1, [] );

            yvi = yvi + rand(size(yvi))*.001-.0005;

            xv = [xv xvi];
            yv = [yv yvi];

            xinterest = [xepts{iseg}(1) xepts{iseg}(end)];
            yinterest = [yepts{iseg}(1) yepts{iseg}(end)];

            if ( abs( xinterest(1) - xinterest(2) ) < 1e-3 )
                imid = ceil( length( xepts{iseg} ) * 0.5 );
                xinterest(2) = xepts{iseg}(imid);
                yinterest(2) = yepts{iseg}(imid);
            end

            for i=1:length(xinterest)
                xgrd = linspace( xinterest(i) - off, xinterest(i) + off, nsurvey );
                ygrd = linspace( yinterest(i) - off, yinterest(i) + off, nsurvey );

                [xg, yg] = meshgrid( xgrd, ygrd );

                xvi = reshape( xg, 1, [] );
                yvi = reshape( yg, 1, [] );

                yvi = yvi + rand(size(yvi))*.001-.0005;

                xv = [xv xvi];
                yv = [yv yvi];
            end
        else
            off = ydisk{iseg} * 0.1;
            xgrd = linspace( xdisk{iseg} - off, xdisk{iseg} + off, nsurvey );
            ygrd = linspace( ydisk{iseg} - off, ydisk{iseg} + off, nsurvey );

            [xg, yg] = meshgrid( xgrd, ygrd );

            xvi = reshape( xg, 1, [] );
            yvi = reshape( yg, 1, [] );

            yvi = yvi + rand(size(yvi))*.001-.0005;

            xv = [xv xvi];
            yv = [yv yvi];
        end
    end

    % Build far-field survey points
    xgrd = linspace( xomin, xomax, 31 );
    ygrd = linspace( yomin, yomax, nsurvey );

    [xg, yg] = meshgrid( xgrd, ygrd );

    xvi = reshape( xg, 1, [] );
    yvi = reshape( yg, 1, [] );
    yvi = yvi + rand(size(yvi))*.001-.0005;

    xv = [xv xvi];
    yv = [yv yvi];

    % Replace points near propwash with aligned points
    for iseg=1:nseg
        if( props{iseg} )

            % First streamtube panel length
            len = sqrt( (xepts{iseg}(2) - xepts{iseg}(1))^2 + (yepts{iseg}(2) - yepts{iseg}(1))^2 );

            ymin = min( yepts{iseg} ) - len * 2.0;
            ymax = max( yepts{iseg} ) + len * 2.0;

            mask = ( xv > xepts{iseg}(1) ) & ( xv < xepts{iseg}(end) ) & ( yv > ymin ) & ( yv < ymax );
            xv = xv(~mask);
            yv = yv(~mask);

            xgrd = xepts{iseg};
            ygrd = linspace( ymin, ymax, nsurvey );

            [xg, yg] = meshgrid( xgrd, ygrd );

            xvi = reshape( xg, 1, [] );
            yvi = reshape( yg, 1, [] );
            yvi = yvi + rand(size(yvi))*.001-.0005;

            xv = [xv xvi];
            yv = [yv yvi];
        end
    end

    inv = zeros( size(xv) );

    for iseg=1:nseg
        if( ~props{iseg} )
            xep = xepts{iseg};
            yep = yepts{iseg};

            % Forces open polys to be closed
            % Signed minimum distance -- negative is inside.
            [dmin, x_d_min, y_d_min, is_vertex, idx_c] = p_poly_dist( xv, yv, xep, yep, true );

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

    [uv, vv] = planarvel( xv, yv, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props );

    xv = xv( ~inv );
    yv = yv( ~inv );
    uv = uv( ~inv );
    vv = vv( ~inv );

    Cedg = [];

    for iseg=1:nseg
        if( ~props{iseg} )
            nnext = length(xv) + 1;
            mcon = length(xcp{iseg}) - 1;

            xv = [xv xcp{iseg}];
            yv = [yv ycp{iseg}];

            uv = [uv gammas{iseg} .* cos( theta{iseg} ) .* sign(dx{iseg}) ];
            vv = [vv gammas{iseg} .* sin( theta{iseg} ) .* sign(dx{iseg}) ];

            Cedg = [Cedg [ nnext:nnext+mcon; nnext+1:nnext+mcon nnext ] ];
        end
    end

    if( isempty(Cedg) )
        DT = delaunayTriangulation( xv', yv' );
        tri = DT.ConnectivityList;
        IO = true( size( tri, 1 ), 1 );
    else
        DT = delaunayTriangulation( xv', yv', Cedg' );
        tri = DT.ConnectivityList;
        IO = ~isInterior(DT);
    end

    Vmagv = sqrt(uv.^2+vv.^2);
    Cpv = 1 - Vmagv.^2;


    for iseg=1:nseg
        if( props{iseg} )
            xst = xepts{iseg};
            yst = yepts{iseg};

            xst = [xst xst(end)+10 xst(end)+10];
            yst = [yst yst(end) 0];                    %%%%%%  Likely crap

            if ( xdisk{jseg} ~= xst(1) )
                xst = [xdisk{jseg} xdisk{jseg} xst];
                yst = [ 0  ydisk{jseg} yst];
            else
                xst = [xst(1) xst];
                yst = [ 0 yst];
            end


            % Forces open polys to be closed
            % Signed minimum distance -- negative is inside.
            [dmin, x_d_min, y_d_min, is_vertex, idx_c] = p_poly_dist( xv, yv, xst, yst, true );

            mask = dmin <= 0.0;

            Cpv(mask) = Cpv(mask) + deltaCP{iseg};

        end
    end
end  % if( drawplots )

% if ( drawplots )
%     verbose = 0;  % flag to report progress
%     maxits = 1e4; % maximum number of iterations
%     Ltol = 0.01;  % flowpath "out-of-triangle" tolerance
%     dLtol = 0.5;  % flowpath "curvature" tolerance
% 
%     ysl = linspace( min(ylim), max(ylim), nsl);
%     % ysl(1) = ysl(2) * 0.1;
%     xsl = min(xlim) * ones( size(ysl) ) + .01;
% 
%     FlowP = TriStream( tri, xv, yv, uv, vv, xsl, ysl, verbose, maxits, Ltol, dLtol );
% 
%     xcap = [];
%     ycap = [];
% 
%     for iseg=1:nseg
%         if( props{iseg} )
%             % Integrate mass flow at disk
%             x0 = xepts{iseg}(1);
%             mass0 = inf;
%             [ ydist, mdist ] = ode45( @planardm, [yemin{iseg}(1), yepts{iseg}(1)], 0, opts, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props, x0, mass0 );
%             mass0 = mdist(end);
% 
%             % Find capture streamtube by continuity
%             x0 = xlim(1) + .01;
%             [ydist, mdist, ye] = ode45( @planardm, ylim, 0, opts, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props, x0, mass0 );
%             y0 = ye;
% 
%             xcap = [xcap x0];
%             ycap = [ycap y0];
% 
%             FlowPStreamtube = TriStream( tri, xv, yv, uv, vv, xcap, ycap, verbose, maxits, Ltol, dLtol );
% 
%             % Build set of all geometry points 'aft' of xdisk
%             xtube = [];
%             ymax = 0;
%             for jseg = 1:nseg
%                 if( ~props{jseg} )  % Skip other props
%                     xtube = unique( [xtube xepts{jseg}] );
%                     ymax = max( ymax, max( yepts{jseg} ) );
%                 end
%             end
% 
%             ymaxtube = 1.1 * ymax * ones( size( xtube ) );
%             ymintube = zeros( size( xtube ) );
% 
%             for jseg = 1:nseg
%                 if( ~props{jseg} )  % Skip other props
% 
%                     xtrim = xepts{jseg};
%                     ytrim = yepts{jseg};
% 
%                     for i=1:length(xtube)
%                         p = InterX( [xtube(i) * [1 1]; [ymintube(i) ymaxtube(i)] ], [xtrim; ytrim] );
%                         x = p(1,:)';
%                         y = p(2,:)';
% 
%                         if( ~isempty(x) )
%                             if( ~kuttas{jseg} ) % Center line body
%                                 ymintube( i ) = max( max( y ), ymintube(i) );
%                             else  % Duct
%                                 ymaxtube( i ) = min( min( y ), ymaxtube(i) );
%                             end
%                         end
%                     end
%                 end
%             end
% 
%             xfarcap = FlowPStreamtube(1).x;
%             yfarcap = FlowPStreamtube(1).y;
% 
%             mask = ymaxtube > ymax;
% 
%             xtubecap = xtube( mask );
%             xtubeduct = xtube( ~mask );
% 
%             if ( ~isempty(xtubeduct) )
%                 mask2 = xfarcap < xtubeduct(1);
%                 xfarcap = xfarcap( mask2 );
%                 yfarcap = yfarcap( mask2 );
%             end
% 
%             ymaxtubecap =  interp1( xfarcap, yfarcap, xtube(mask) );
%             ymintubecap = ymintube( mask );
%             ymaxtube(mask) = nan;
% 
%             if ( ~isempty(xtubecap) )
%                 mask3 = xfarcap < xtubecap(1);
%                 xfarcap = xfarcap( mask3 );
%                 yfarcap = yfarcap( mask3 );
%             end
% 
%             Atubeup = pi * ( ymaxtube.^2 - ymintube.^2 );
%             Atubedn = pi * ( yepts{iseg}.^2 - yemin{iseg}.^2 );
%             Atubecap = pi * ( ymaxtubecap.^2 - ymintubecap.^2 );
%             Afarcap = pi * yfarcap.^2;
%         end
%     end
% end

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

        CDi{iseg} = - 2.0 * pi * sum( Cp{iseg} .* sin( theta{iseg} ) .* ycp{iseg} .* ds{iseg} ) / Sref;
        CD = CD + CDi{iseg};

        disp(['CD on ' names{iseg} ':  ' num2str(CDi{iseg})])
    end
end

disp(['CD total:  ' num2str(CD)])

if( drawplots )

    figure(5)
    quiver( xv, yv, uv, vv )
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, yepts{iseg} );
    end
    hold off
    axis equal


    figure(6)
    trisurf( tri(IO,:), xv, yv, zeros(size(xv)), Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, yepts{iseg} );
    end
    hold off
    axis equal
    view(0,90)
    title('v/Vinf')


    figure(7)
    trisurf( tri(IO,:), xv, yv, zeros(size(xv)), Cpv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, yepts{iseg} );
    end
    hold off
    axis equal
    view(0,90)
    title('Cp')

    figure(8)
%    PlotTriStream( FlowP );
    if( exist( 'FlowPStreamtube', 'var' ) )
        PlotTriStream( FlowPStreamtube );
    end
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, yepts{iseg} );
    end
    hold off
    axis equal

    figure(9)
    trisurf( tri(IO,:), xv, yv, zeros(size(xv)), Vmagv, 'LineStyle', 'none', 'FaceColor', 'interp' );
    hold on
%    PlotTriStream( FlowP, 'k' );
    hold on;
    if( exist( 'FlowPStreamtube', 'var' ) )
        PlotTriStream( FlowPStreamtube, 'k' );
    end
    hold on;  % PlotTriStream turns hold off.
    for iseg=1:nseg
        plot( xepts{iseg}, yepts{iseg},'k' );
    end
    hold off
    axis equal
    view(0,90)
    title('v/Vinf')


%     xsurvey = [-1 0 1 2 4];
%
%     for i=1:length(xsurvey);
%         mask = ( xv == xsurvey(i) );
%
%         r = rv(mask);
%         v = vv(mask);
%         u = uv(mask);
%
%         [r, indx] = sort( r );
%         v = v(indx);
%         u = u(indx);
%
%         figure(10)
%         plot( v, r )
%         hold on
%
%         figure(11)
%         plot( u, r );
%         hold on
%     end
%     figure(10)
%     hold off
%     figure(11)
%     hold off

    figure(12)
    trimesh( tri(IO,:), xv, yv, zeros(size(xv)) );
    hold on
    for iseg=1:nseg
        plot( xepts{iseg}, yepts{iseg},'k','LineWidth', 2.0 );
    end
    axis equal
    view(0,90)


    if( exist( 'xtube', 'var' ) )

        figure(13)
        plot( xtube, ymintube, xtube, ymaxtube, xtubecap, ymaxtubecap, xfarcap, yfarcap, 'LineWidth', lw )
        hold on
        for iseg = 1:nseg
             if ( props{iseg} )
                plot( xepts{iseg}, yepts{iseg}, xepts{iseg}, yemin{iseg}, 'LineWidth', lw )
            end
        end
        hold off
        axis equal

        figure(14)
        plot( xtube, Atubeup, 'LineWidth', lw)
        hold on
        plot( xtubecap, Atubecap, '--', 'LineWidth', lw )
        plot( xfarcap, Afarcap, '--', 'LineWidth', lw )
        for iseg = 1:nseg
            if ( props{iseg} )
                figure(14)
                plot( xepts{iseg}, Atubedn, '--', 'LineWidth', lw )
            end
        end
        hold off
        ax = axis;
        ax(3) = 0;
        axis(ax);
        ylabel('Streamtube Area')

    end

end
