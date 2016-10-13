function yprm = dm( r, mass, W, xv, rv, gammav, xt, rt, gammat, x, mass0 )

yprm = zeros( size( mass ) );

% i loop over vortex rings implied by . operations
[ uv, ~ ] = ringvortex( xv, rv, x, r );

% i loop over vortex tubes implied by . operations
[ ut, ~ ] = tubevortex( xt, rt, x, r );

uv = W + sum( gammav .* uv ) + sum( gammat .* ut );

yprm(1) = uv * 2.0 * pi * r;


end
