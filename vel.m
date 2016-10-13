function yprm = vel( t, y, W, xv, rv, gammav, xt, rt, gammat )

x = y(1);
r = y(2);

yprm = zeros( size( y ) );

% i loop over vortex rings implied by . operations
[ uv, vv ] = ringvortex( xv, rv, x, r );

% i loop over vortex tubes implied by . operations
[ ut, vt ] = tubevortex( xt, rt, x, r );

yprm(1) = W + sum( gammav .* uv ) + sum( gammat .* ut );
yprm(2) = sum( gammav .* vv ) + sum( gammat .* vt );

end
