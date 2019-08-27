clear all
format compact

th2 = 10.0;

yoff2 = 0.5;


drawplots = true;

% n = 15;
% xs = linspace( -2, 2, n );
xs = -3:.1:3;
n = length(xs);


for i=1:n
    xoff2 = xs(i);

    [CLi, CDi] = runAF_fun2( xoff2, yoff2, th2, drawplots );
    
    CLis(i,:) = cell2mat( CLi );
    CDis(i,:) = cell2mat( CDi );

    figure(25)
    print('-dpng','-r300',['Cp_UT_Sweep_' num2str(i) '.png'])
    
    F(i) = getframe(25);
end

close all
figure
axis off
movie(flip(F))


v = VideoWriter('boom.avi','Uncompressed AVI');
v.FrameRate = 15;
open(v)
writeVideo(v,flip(F))
close(v)

figure(1)
plot(xs,CLis,xs,sum(CLis,2))

xlabel('Airfoil Position')
ylabel('C_L')
legend('Cylinder','Foil','Total')


figure(2)
plot(xs,CDis,xs,sum(CDis,2))

xlabel('Airfoil Position')
ylabel('C_D')
legend('Cylinder','Foil','Total')