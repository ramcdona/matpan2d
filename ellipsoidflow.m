function [pot, u, v, w]=ellipsoidflow(x, y, z, V, a, b, c)

afun = @(t) 1.0./((a^2+t).*sqrt((a^2+t).*(b^2+t).*(c^2+t)));
bfun = @(t) 1.0./((b^2+t).*sqrt((a^2+t).*(b^2+t).*(c^2+t)));
cfun = @(t) 1.0./((c^2+t).*sqrt((a^2+t).*(b^2+t).*(c^2+t)));

alpha=a*b*c*quad(afun,0,1e14,1e-14);
beta=a*b*c*quad(bfun,0,1e14,1e-14);
gamma=a*b*c*quad(cfun,0,1e14,1e-14);

k1=alpha/(2.0-alpha);
k2=beta/(2.0-beta);
k3=gamma/(2.0-gamma);

k=[k1 k2 k3];
ABC=k+1;

Vmax=ABC.*V;

% Velocity potential.
% Only valid on surface.
pot=Vmax(1)*x + Vmax(2)*y + Vmax(3)*z;

% Normal vector
nx=2.0*x/a^2;
ny=2.0*y/b^2;
nz=2.0*z/c^2;

nmag=(nx.^2+ny.^2+nz.^2).^.5;
nx=nx./nmag;
ny=ny./nmag;
nz=nz./nmag;

% Vmax component in panel normal direction.
Vnormal=Vmax(1)*nx+Vmax(2)*ny+Vmax(3)*nz;

% Surface velocity is remainder.
u=Vmax(1)-Vnormal.*nx;
v=Vmax(2)-Vnormal.*ny;
w=Vmax(3)-Vnormal.*nz;