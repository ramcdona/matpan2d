clear all
format compact
close all

load('xcp.mat')
load('props.mat')
load('Cp0.mat')
load('Cp500.mat')
load('Cp700.mat')
load('Cp785.mat')

nseg = length(xcp);

figure(1)
clf
hold on
for iseg=1:nseg
    if( ~props{iseg} )
        plot( xcp{iseg}, -Cp0{iseg} );
    end
end
hold off
ylabel('-C_p Compressible')

figure(2)
clf
hold on
xsh = min(xcp{1});
dx = max(xcp{1})-min(xcp{1});

plot( (xcp{1}-xsh)./dx, Cp0{1}, '--', 'LineWidth', 1 );
plot( (xcp{1}-xsh)./dx, Cp500{1}, 'k--', 'LineWidth', 1 );
plot( (xcp{1}-xsh)./dx, Cp700{1}, 'g--', 'LineWidth', 1 );
plot( (xcp{1}-xsh)./dx, Cp785{1}, 'r--', 'LineWidth', 1 );
hold off
ylabel('C_p Compressible')
legend('M=0','M=0.5','M=0.7','M=0.785') %

ax=axis;
ax(3)=-1.25;
ax(4)=1.25;
axis(ax)


figure(3)
clf
hold on
plot( xcp{2}, Cp0{2}, '-', 'LineWidth', 1 );
plot( xcp{2}, Cp500{2}, 'k-', 'LineWidth', 1 );
plot( xcp{2}, Cp700{2}, 'g-', 'LineWidth', 1 );
plot( xcp{2}, Cp785{2}, 'r-', 'LineWidth', 1 );
hold off
ylabel('C_p Compressible')
legend('M=0','M=0.5','M=0.7','M=0.785') %

ax=axis;
ax(3)=-1.25;
ax(4)=1.25;
axis(ax)

figure(4)
clf
hold on

xsh = min(xcp{1});
dx = max(xcp{1})-min(xcp{1});

plot( (xcp{1}-xsh)./dx, Cp0{1}, '--', 'LineWidth', 1 );
plot( (xcp{1}-xsh)./dx, Cp500{1}, 'k--', 'LineWidth', 1 );
plot( (xcp{1}-xsh)./dx, Cp700{1}, 'g--', 'LineWidth', 1 );
plot( (xcp{1}-xsh)./dx, Cp785{1}, 'r--', 'LineWidth', 1 );

xsh = min(xcp{2});
dx = max(xcp{2})-min(xcp{2});

plot( (xcp{2}-xsh)./dx, Cp0{2}, '-', 'LineWidth', 1 );
plot( (xcp{2}-xsh)./dx, Cp500{2}, 'k', 'LineWidth', 1 );
plot( (xcp{2}-xsh)./dx, Cp700{2}, 'g', 'LineWidth', 1 );
plot( (xcp{2}-xsh)./dx, Cp785{2}, 'r', 'LineWidth', 1 );
hold off
ylabel('C_p Compressible')
legend('M=0','M=0.5','M=0.7','M=0.785') %

ax=axis;
ax(1)=0;
ax(2)=1;
ax(3)=-1.25;
ax(4)=1.25;
axis(ax)


