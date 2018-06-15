function [value, isterminal, direction] = planarmassevt( eta, mass, W, nseg, xepts, yepts, xuppts, yuppts, xlowpts, ylowpts, xcp, ycp, xupcp, yupcp, xlowcp, ylowcp, gammas, gammasup, gammaslow, ds, dsup, dslow, dx, dxup, dxlow, props, x0, y0, xv, yv, len, mass0  )

value = mass - mass0;
isterminal = 1;
direction = 0;

end
