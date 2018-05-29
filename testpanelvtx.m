clear all
format compact

xstart = 1.5;
ystart = 0.5;

xend = 0.75;
yend = -.75;

xdir = 1;
ydir = 0;
big = 100;
n=31;

xend = xstart + xdir * big;
yend = ystart + ydir * big;


xs = linspace(0,3,140);
ys = linspace(0,1,130);

[X,Y] = meshgrid(xs,ys);

xv = reshape( X, numel(X), 1 );
yv = reshape( Y, numel(Y), 1 );

% [ u, v ] = panelvortex( xstart, ystart, xend, yend, xv, yv );
 [ u, v ] = rayvortex( xstart, ystart, xdir, ydir, xv, yv );

% u = zeros(size(xv));
% v = u;
% for i = 1:n
%     [ du, dv ] = ptvortex( xstart + (i-1) * xdir * big / (n-1)  , ystart + (i-1) * ydir * big / (n-1) , xv, yv );
%     u = u + du;
%     v = v + dv;
% end

figure(1)
quiver( xv, yv, u, v )
hold on
%plot([xstart xend],[ystart yend])
plot([xstart xstart + big*xdir],[ystart ystart+big*ydir],'LineWidth',3)
hold off
ax=axis;

figure(2)
[c,h]=contourf(X,Y,reshape( sqrt(u.^2+v.^2), size(X)),31,'LineStyle','none');
clabel(c,h);
hold on
%plot([xstart xend],[ystart yend])
%plot([xstart xstart + big*xdir],[ystart ystart+big*ydir],'LineWidth',3)
hold off
colorbar

axis off


