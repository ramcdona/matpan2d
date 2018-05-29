function [u, v] = planarvel( x, y, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props )

% Handle x, y of any dimension (scalar, vector, matrix) by iterating
% through a linear index 1:numel(x).  This makes it impossible to get any
% speedup by treating x, y as vectors -- however it is assumed that the
% problem is made of many point vortices which may be treated as vectors.

nsurvey = numel( x );
u = W * ones( size(x) );
v = zeros( size(x) );

% Explicit isurvey loop over survey points
for isurvey = 1:nsurvey
    uv = W;
    vv = 0;

    for iseg=1:nseg

        mask = true( size( gammas{iseg} ) );
        if ( props{iseg} )
            mask(end/2) = false;
            mask(end) = false;
        end

        [ uj, vj ] = ptvortex( xcp{iseg}, ycp{iseg}, x(isurvey), y(isurvey) );
        uv = uv + sum( uj .* gammas{iseg}(mask) .* ds{iseg} .* sign(dx{iseg}) );
        vv = vv + sum( vj .* gammas{iseg}(mask) .* ds{iseg} .* sign(dx{iseg}) );

        if( props{iseg} )
            xv = xepts{iseg}(end/2) - xepts{iseg}(end/2 - 1);
            yv = yepts{iseg}(end/2) - yepts{iseg}(end/2 - 1);
            [ uj, vj ] = rayvortex( xepts{iseg}(end/2), yepts{iseg}(end/2), xv, yv, x(isurvey), y(isurvey) );
            uv = uv + uj * gammas{iseg}(end/2);
            vv = vv + vj * gammas{iseg}(end/2);

            xv = xepts{iseg}(end) - xepts{iseg}(end - 1);
            yv = yepts{iseg}(end) - yepts{iseg}(end - 1);
            [ uj, vj ] = rayvortex( xepts{iseg}(end), yepts{iseg}(end), xv, yv, x(isurvey), y(isurvey) );
            uv = uv + uj * gammas{iseg}(end);
            vv = vv + vj * gammas{iseg}(end);
        end
    end

    u(isurvey) = uv;
    v(isurvey) = vv;
end

end
