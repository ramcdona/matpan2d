function [u, v] = planarvel( x, y, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props )

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
            mask = true( size( gammasup{iseg} ) );
            mask(end) = false;

            [ uj, vj ] = ptvortex( xupcp{iseg}, yupcp{iseg}, x(isurvey), y(isurvey) );

            uv = uv + sum( uj .* gammasup{iseg}(mask) .* dsup{iseg} .* sign(dxup{iseg}) );
            vv = vv + sum( vj .* gammasup{iseg}(mask) .* dsup{iseg} .* sign(dxup{iseg}) );


            [ uj, vj ] = ptvortex( xlowcp{iseg}, ylowcp{iseg}, x(isurvey), y(isurvey) );

            uv = uv - sum( uj .* gammaslow{iseg}(mask) .* dslow{iseg} .* sign(dxlow{iseg}) );
            vv = vv - sum( vj .* gammaslow{iseg}(mask) .* dslow{iseg} .* sign(dxlow{iseg}) );



            xvu = xuppts{iseg}(end) - xuppts{iseg}(end - 1);
            yvu = yuppts{iseg}(end) - yuppts{iseg}(end - 1);
            dsu = sqrt( xvu.^2 + yvu.^2 );
            xvu = xvu ./ dsu;
            yvu = yvu ./ dsu;

            xvl = xlowpts{iseg}(end) - xlowpts{iseg}(end - 1);
            yvl = ylowpts{iseg}(end) - ylowpts{iseg}(end - 1);
            dsl = sqrt( xvl.^2 + yvl.^2 );
            xvl = xvl ./ dsl;
            yvl = yvl ./ dsl;

            xv = 0.5 * (xvu + xvl);
            yv = 0.5 * (yvu + yvl);
            ds = sqrt( xv.^2 + yv.^2 );
            xv = xv ./ ds;
            yv = yv ./ ds;


            [ uj, vj ] = rayptvortex( xuppts{iseg}(end), yuppts{iseg}(end), xv, yv, x(isurvey), y(isurvey) );

            uv = uv + uj * gammasup{iseg}(end) * sign(xv);
            vv = vv + vj * gammasup{iseg}(end) * sign(xv);

            [ uj, vj ] = rayptvortex( xlowpts{iseg}(end), ylowpts{iseg}(end), xv, yv, x(isurvey), y(isurvey) );

            uv = uv - uj * gammaslow{iseg}(end) * sign(xv);
            vv = vv - vj * gammaslow{iseg}(end) * sign(xv);
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
