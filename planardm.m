function dm = planardm( eta, mass, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props, x0, y0, xv, yv, len, mass0 )

[u, v] = planarvel( x0 + xv * eta * len, y0 + yv * eta * len, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props );

dm = ( u * yv + v * xv ) * len;

end
