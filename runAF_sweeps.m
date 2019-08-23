clear all
format compact
close all

% Cases of interest
% runAF_fun( 0.05, .16, 0, true );

% Airfoil type
% 1 NACA 4-digit
% 2 VR12
% 3 VKT
aftype = 2;

th = 5;

nx = [1 2 3]*25;
ny = [1 2 3]*26;

xmin = [-2 -1.5 -1.07];
xmax = [2 1.3 1.07];

ymin = [-2 -0.5 -0.13];
ymax = [2 0.5 0.13];

xv = [];
yv = [];
for i=1:length(xmin)
    xgrd = linspace( xmin(i), xmax(i), nx(i) );
    ygrd = linspace( ymin(i), ymax(i), ny(i) );

    [xg, yg] = meshgrid( xgrd, ygrd );

    xvi = reshape( xg, 1, [] );
    yvi = reshape( yg, 1, [] );

    xv = [xv xvi];
    yv = [yv yvi];
end
zv = -0.1 * ones(size(xv));

% yv = yv + rand( size(yv) ) * 0.001 - 0.0005;

DT = delaunayTriangulation( xv', yv' );
tri = DT.ConnectivityList;

for i=1:length(xv)

    [CLi, CDi] = runAF_fun( xv(i), yv(i), th, false, aftype );


    CLis(i,:) = cell2mat( CLi );
    CDis(i,:) = cell2mat( CDi );
end

[~, ~, xep, yep] = runAF_fun( 10, 10, th, false, aftype );

CLtot = sum(CLis,2);
CDtot = sum(CDis,2);



figure(1)
trisurf(tri,xv,yv,zv,CLtot,'EdgeColor','none','FaceColor','interp')
view(0,90)
axis equal
ax=axis;
hold on
if ( max(CLtot) > 0 && min(CLtot) < 0 )
    [c,h]=tricontour(tri,xv,yv,CLtot,[0 0]);
    set(h,'EdgeColor','k')
end
fill(xep,yep,'k')
hold off
view(0,90)
if ( th > 3 )
    caxis([0.5,1.5])
else
    caxis([-0.3 0.3])
end
colorbar
axis equal
axis(ax)
axis off
title('Total C_L')

figure(2)
trisurf(tri,xv,yv,zv,CLis(:,1),'EdgeColor','none','FaceColor','interp')
view(0,90)
axis equal
hold on
[c,h]=tricontour(tri,xv,yv,CLis(:,1),[0 0]);
set(h,'EdgeColor','k')
fill(xep,yep,'k')
hold off
view(0,90)
caxis([-6,6])
colorbar
axis equal
axis(ax)
axis off
title('C_L Fixed AF')

figure(3)
trisurf(tri,xv,yv,zv,CDis(:,1),'EdgeColor','none','FaceColor','interp')
view(0,90)
axis equal
ax=axis;
hold on
[c,h]=tricontour(tri,xv,yv,CDis(:,1),[0 0]);
set(h,'EdgeColor','k')
fill(xep,yep,'k')
hold off
caxis([-.3 .3])
colorbar
axis equal
axis(ax);
axis off
title('C_D Fixed AF')

figure(4)
trisurf(tri,xv,yv,zv,CDtot','EdgeColor','none','FaceColor','interp')
hold on
fill(xep,yep,'k')
hold off
view(0,90)
caxis([0 .0001])
colorbar
axis equal
axis(ax);
axis off
title('Total C_D')

figure(5)
trisurf(tri,xv,yv,zv,nan(size(xv)),'EdgeColor','none','FaceColor','interp')
hold on
plot(xep,yep,'k','LineWidth',1.5)
hold off
view(0,90)
caxis([0 .0001])
colorbar
axis equal
axis(ax);
axis off
title('Foil')


figure(6)
trisurf(tri,xv,yv,zv,CLtot-2*0.602829261927623,'EdgeColor','none','FaceColor','interp')
view(0,90)
axis equal
ax=axis;
hold on
if ( max(CLtot-2*0.602829261927623) > 0 && min(CLtot-2*0.602829261927623) < 0 )
    [c,h]=tricontour(tri,xv,yv,CLtot-2*0.602829261927623,[0 0]);
    set(h,'EdgeColor','k')
end
fill(xep,yep,'k')
hold off
view(0,90)
if ( th > 3 )
    caxis([0.5,1.5]-1.2)
else
    caxis([-0.3 0.3])
end
colorbar
axis equal
axis(ax)
axis off
title('Total C_L- 2 C_L^*')


for i=1:6
    figure(i)
    print('-dpng','-r300',['af_int_' num2str(i) '.png']);
end

