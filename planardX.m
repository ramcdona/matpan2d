function Xdot = planardX( t, X, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props, s0 )

x = X(1);
y = X(2);
s = X(3);

[u, v] = planarvel( x, y, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props );

Xdot = nan(3,1);
Xdot(1) = u;
Xdot(2) = v;
Xdot(3) = sqrt( u.^2 + v.^2 );

end
