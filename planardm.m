function dm = planardm( eta, mass, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props, x0, y0, xv, yv, len, mass0 )

[u, v] = planarvel( x0 + xv * eta, y0 + yv * eta, W, nseg, xepts, yepts, xcp, ycp, gammas, ds, dx, props );

dm = ( u * yv + v * xv ) * len;

end
