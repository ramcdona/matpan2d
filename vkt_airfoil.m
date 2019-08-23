% Karman-Trefftz airfoil
%
%   ell = quater Chord length (chord=4*ell)
%   alpha_ml = Angle of attack (radian)
%   epsilon = Thisckness
%   kappa = Camber
%   tau = Trailing edge angle (radian)
%   nx = # of nodes in the circumferential direction.
%
% Based loosely on work of Katate Masatsuka, Januray 2010.
% Available from http://www.cfdbooks.com

function [x, y, cp] = vkt_airfoil(nx, alpha, epsilon, kappa, tau )

%Chord length = 2*ell
ell = 0.25;


x=zeros(nx, 1);
y=zeros(nx, 1);
psi=zeros(nx, 1);
phi=zeros(nx, 1);
u=zeros(nx, 1);
v=zeros(nx, 1);
cp=zeros(nx, 1);


%Imaginray number
i = complex(0.0,1.0);

%Radius of the circle
a = ell.*sqrt((1.0+epsilon).^2 + kappa.^2 );

if(ell.*kappa./a > 1.0)
    disp( ' Camber parameter, kappa, is too large. Reduce kappa, and try again.');
end

%Angle of the TE location - Typo in Eq. (7.10.58)
beta = asin(ell.*kappa./a);
%Parameter related to TE angle
n = 2.0 - tau./pi;
%Center of the circle.
mu = complex(-ell.*epsilon,ell.*kappa);

% 1. Generate a KT-airfoil: nodes ordered clockwise from TE

for j = 1: nx
    theta = (j-1).*(2.0.*pi)./(nx-1);
    
    % Coordinates in the circle plane
    
    % Theta=0 corresponds to TE.
    xi = a.*cos(theta-beta) + real(mu);
    eta = a.*sin(theta-beta) + imag(mu);
    zeta = complex(xi,eta);
    
    % Coordinates in the airfoil plane
    
    temp =(zeta-ell).^n./(zeta+ell).^n;
    %Karman-Trefftz transformation
    z = n.*ell.*( 1.0 + temp )./( 1.0 - temp );
    x(j) = real(z);
    y(j) = imag(z);
    
    % f(zeta): Complex potential in the circle plane.
    
    [f ]=cmplx_potential(zeta, alpha,beta,a,mu,i);
    % Scalar potential (can be discontinuous for lifting flows)
    phi(j) = real(f);
    % Stream function
    psi(j) = imag(f);
    
    % w(zeta): Complex velocity in the circle plane (a flow around a cylinder)
    
    [w]=cmplx_velocity(zeta, alpha,beta,a,mu,i);
    
    % Compute the velocity in the airfoil plane: (u,v) = w/(dZ/dzeta)
    
    % Derivative of the Karman-Trefftz transformation.
    
    [dzdzeta]=derivative(zeta, ell,n);
    
    % Special treatment at the trailing edge(theta=zero).
    if( theta == 0.0 )
        %  (1)Joukowski airfoil (cusped trailing edge: tau = 0.0)
        if(tau == 0.0)
            u(j) =  real( ell./a.*exp(2.0.*i.*beta).*cos(alpha+beta) );
            v(j) = -imag( ell./a.*exp(2.0.*i.*beta).*cos(alpha+beta) );
            cp(j) =  1.0 -( u(j).^2 + v(j).^2 );
            %  (2)Karman-Trefftz airfoil (finite angle: tau > 0.0)
            %     TE must be a stagnation point.
        else
            u(j) =  0.0;
            v(j) =  0.0;
            cp(j) =  1.0;
        end
        % Elsewhere (anywhere not TE): (u,v) = w/(dZ/dzeta)
    else
        u(j) = real(w./dzdzeta);
        v(j) =-imag(w./dzdzeta);
        cp(j) = 1.0 -( u(j).^2 + v(j).^2 );
    end
    
end %  airfoil_boundary;

% Reverse point order to match pan3 expectations
x = flipud(x);
y = flipud(y);
cp = flipud(cp);


end

% Complex potential in the circle plane (zeta-plane):
function [cmplx_potentialresult]=cmplx_potential(zeta, alpha_ml,beta,a,mu,i)
cmplx_potentialresult = zeta.*exp(-i.*alpha_ml)+ i.*2.0.*a.*sin(alpha_ml+beta).*log(zeta-mu)+ a.*a.*exp(i.*alpha_ml)./(zeta-mu);
return
end

% Complex velocity in the circle plane (zeta-plane):
function [cmplx_velocityresult]=cmplx_velocity(zeta, alpha_ml,beta,a,mu,i)
cmplx_velocityresult = exp(-i.*alpha_ml)+ i.*2.0.*a.*sin(alpha_ml+beta)./(zeta-mu)- a.*a.*exp(i.*alpha_ml)./((zeta-mu).^2 );
return
end

% Derivative of the KT transformation:
function [derivativeresult]=derivative(zeta, ell,n)
derivativeresult = 4.0.*(n.*ell).^2.*((zeta+ell).*(zeta-ell)).^(n-1.0)./((zeta+ell).^n -(zeta-ell).^n ).^2;
return
end
