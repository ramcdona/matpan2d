function [u, v] = planarvel( x, y, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props );

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

        if ( props{iseg} )
            [ uj, vj ] = ptvortex( xupcp{iseg}, yupcp{iseg}, x(isurvey), y(isurvey) );

            uv = uv + sum( uj .* gammasup{iseg} .* dsup{iseg} .* sign(dxup{iseg}) );
            vv = vv + sum( vj .* gammasup{iseg} .* dsup{iseg} .* sign(dxup{iseg}) );

            xv = xuppts{iseg}(end) - xuppts{iseg}(end - 1);
            yv = yuppts{iseg}(end) - yuppts{iseg}(end - 1);
            [ uj, vj ] = rayvortex( xuppts{iseg}(end), yuppts{iseg}(end), xv, yv, x(isurvey), y(isurvey) );
            uv = uv + uj * gammasup{iseg}(end);
            vv = vv + vj * gammasup{iseg}(end);

            [ uj, vj ] = ptvortex( xlowcp{iseg}, ylowcp{iseg}, x(isurvey), y(isurvey) );

            uv = uv + sum( uj .* gammaslow{iseg} .* dslow{iseg} .* sign(dxlow{iseg}) );
            vv = vv + sum( vj .* gammaslow{iseg} .* dslow{iseg} .* sign(dxlow{iseg}) );

            xv = xlowpts{iseg}(end) - xlowpts{iseg}(end - 1);
            yv = ylowpts{iseg}(end) - ylowpts{iseg}(end - 1);
            [ uj, vj ] = rayvortex( xlowpts{iseg}(end), ylowpts{iseg}(end), xv, yv, x(isurvey), y(isurvey) );
            uv = uv + uj * gammaslow{iseg}(end);
            vv = vv + vj * gammaslow{iseg}(end);
        else
            [ uj, vj ] = ptvortex( xcp{iseg}, ycp{iseg}, x(isurvey), y(isurvey) );

            uv = uv + sum( uj .* gammas{iseg} .* ds{iseg} .* sign(dx{iseg}) );
            vv = vv + sum( vj .* gammas{iseg} .* ds{iseg} .* sign(dx{iseg}) );
        end

    end

    u(isurvey) = uv;
    v(isurvey) = vv;
end

end
