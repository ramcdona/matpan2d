function [value, isterminal, direction] = planarslevt(  t, X, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props, s0 )

x = X(1);
y = X(2);
s = X(3);

value = s - s0;
isterminal = 1;
direction = 0;

end
