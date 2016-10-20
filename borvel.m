function [u, v] = borvel( x, r, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props )

% Handle x, r of any dimension (scalar, vector, matrix) by iterating
% through a linear index 1:numel(x).  This makes it impossible to get any
% speedup by treating x, r as vectors -- however it is assumed that the
% problem is made of many ring vortices which may be treated as vectors.

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
            mask(end) = false;
        end

        [ uj, vj ] = ringvortex( xcp{iseg}, rcp{iseg}, x(isurvey), r(isurvey) );
        uv = uv + sum( uj .* gammas{iseg}(mask) .* ds{iseg} .* sign(dx{iseg}) );
        vv = vv + sum( vj .* gammas{iseg}(mask) .* ds{iseg} .* sign(dx{iseg}) );

        if( props{iseg} )
            [ uj, vj ] = tubevortex( xepts{iseg}(end), repts{iseg}(end), x(isurvey), r(isurvey) );
            uv = uv + uj * gammas{iseg}(end);
            vv = vv + vj * gammas{iseg}(end);
        end
    end

    u(isurvey) = uv;
    v(isurvey) = vv;
end

end
