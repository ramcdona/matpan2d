function dm = bordm( r, mass, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props, x, mass0 )

[u, ~] = borvel( x, r, W, nseg, xepts, repts, xcp, rcp, gammas, ds, dx, props );

dm = u * 2.0 * pi * r;

end
