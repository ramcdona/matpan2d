clear all
format compact

th2 = 10.0;

drawplots = false;

nx = 55;
ny = 51;

xs = linspace( -3, 3, nx );
ys = linspace( -3, 3, ny );


for i=1:nx
    xoff2 = xs(i);

    for j=1:ny
        yoff2 = ys(j);
        [CLi, CDi] = runAF_fun2( xoff2, yoff2, th2, drawplots );
    
        CLis(i,j,:) = cell2mat( CLi );
        CDis(i,j,:) = cell2mat( CDi );
    end
end

[xep, yep] = setupellipse( 101, 0.5, 0.5, 0, 0 );


figure(1)
surf(xs,ys,sum(CLis,3)','LineStyle','none')
hold on
plot(xep,yep,'k')
hold off

view(0,90)
axis equal
axis off
colorbar
caxis([0 4])
