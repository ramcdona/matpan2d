function [value, isterminal, direction] = velevt( t, y, W, xv, rv, gammav, xt, rt, gammat )

x = y(1);
r = y(2);

xend = 3;

value = x - xend;
isterminal = 1;
direction = 0;

end
