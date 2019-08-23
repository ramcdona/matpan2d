clear all
format compact
close all



% Airfoil type
% 1 NACA 4-digit
% 2 VR12
% 3 VKT
aftype = 2;

R = 1.108;
C = 0.08;

r = 0.75 * R / C;

alpha = 8;

ysep = 0.73;
ysep = 0.3;

n=501;
phis = linspace( -90, 90, n );
xsep = -r * phis*pi/180;

xsep = [-2:.25:2];
phis = -( xsep / r ) * 180 / pi;
n=length(xsep);

for i = 1 : n

    [CLi, CDi] = runAF_fun( xsep(i), ysep, alpha, true, aftype );


    CLis(i,:) = cell2mat( CLi );
    CDis(i,:) = cell2mat( CDi );
    
    figure(25)
    print('-dpng','-r300',['Cp_UT_Sweep_' num2str(i) '.png'])
end



figure(11)

plot( phis, CLis, phis, sum(CLis,2)/2.0);



figure(12)
plot( phis, sum(CLis,2)/2.0/6.0/2.0, 'Color', [30 119 180]./255, 'LineWidth', 1.5 );

hold on
plot( phis, CLis(:,1)/6.0/2.0, 'Color',[255 127 15]./255, 'LineWidth', 1.5 )
plot( phis, CLis(:,2)/6.0/2.0, 'Color', [43 160 43]./255, 'LineWidth', 1.5  )
hold off

legend('Total','Lower','Upper')
xlabel('Index Angle [deg]')
ylabel('C_T/\sigma')

axis([-90 90 0.0 0.16])


