function [x,y]=naca4(dig1,dig2,dig34,n,ptspace)

if ptspace==1
    % Full period cosine spacing
    theta=linspace(0,-2*pi,n)';  % Cosine point spacing
    xc=(cos(theta)+1)/2;
elseif ptspace==2
    % Half period cosine spacing
    theta=linspace(0,-pi,n)';
    xc=(cos(theta-pi/2)+1);
elseif ptspace==3
    % Linear point spacing
    xc=linspace(1,0,ceil(n/2))';
    ytemp=-sqrt(1-(2*xc-1).^2);
    xc=[xc;flipud(xc(1:end-1))];
    ytemp=[ytemp;-flipud(ytemp(1:end-1))];
    theta=atan2(ytemp,2*xc-1);
    theta(floor(n/2)+1:end)=theta(floor(n/2)+1:end)-2*pi;
end


dig1=dig1/100;
dig2=dig2/10;
dig34=dig34/100;


% % Thickness function (original form)
% a0=1.4845;
% a1=-0.6300;
% a2=-1.7580;
% a3=1.4215;
% a4=-0.5075;
% ytc=dig34*(a0*sqrt(xc)+a1*xc+a2*xc.^2+a3*xc.^3+a4*xc.^4);

% % Classical coefficients * 5
% a0 = 0.2969;
% a1 = -0.1260;
% a2 = -0.3516;
% a3 = 0.2843;
% a4 = -0.1015;
% ytc=5.0*dig34*(a0*sqrt(xc)+a1*xc+a2*xc.^2+a3*xc.^3+a4*xc.^4);

% Closed TE coefficients * 5
a0 = 0.2983340510408059;
a1 = -0.1324560937789966;
a2 = -0.3285782439243375;
a3 = 0.2442055905428443;
a4 = -0.0815053038803161;
ytc=5.0*dig34*(a0*sqrt(xc)+a1*xc+a2*xc.^2+a3*xc.^3+a4*xc.^4);



%Camber line
t1=find(xc<dig2);
t2=find(xc>=dig2);
ycc=zeros(size(xc));
dyccdx=zeros(size(xc));
if dig1~=0
    ycc(t1)=dig1/dig2^2*(2*dig2*xc(t1)-xc(t1).^2);
    dyccdx(t1)=2*dig1/dig2^2*(dig2-xc(t1));

    ycc(t2)=dig1/(1-dig2)^2*(1-2*dig2+2*dig2*xc(t2)-xc(t2).^2);
    dyccdx(t2)=2*dig1/(1-dig2)^2*(dig2-xc(t2));
end

%Camber line slope
theta2=atan(dyccdx);

%Airfoil coordinates
t4=1:ceil(size(xc,1)/2);%find(theta<-pi);
t3=ceil(size(xc,1)/2)+1:size(xc,1); %find(theta>=-pi);
yc=zeros(size(xc));

xc(t3)=xc(t3)-ytc(t3).*sin(theta2(t3));
yc(t3)=ycc(t3)+ytc(t3).*cos(theta2(t3));

xc(t4)=xc(t4)+ytc(t4).*sin(theta2(t4));
yc(t4)=ycc(t4)-ytc(t4).*cos(theta2(t4));

% Rescale
x=xc/max(xc);
y=yc/max(xc);
