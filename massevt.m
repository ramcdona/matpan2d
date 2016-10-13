function [value, isterminal, direction] = massevt( r, mass, W, xv, rv, gammav, xt, rt, gammat, x, mass0 )

value = mass - mass0;
isterminal = 1;
direction = 0;

end
