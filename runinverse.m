clear all
format compact
close all

% Freestream velocity
W = 1.0;

Sref = 1.0;

xlim = [-0.5 2];
rlim = [0 1.5];
nsurvey = 51;

nsl = 23;

drawplots = false;

ntstep = 1;


% Default to no kutta condition, turn on inside runcase setup.
kutta = false;
jtelow = -1;
jteup = -1;


% Construct example body
ncirc = 10*10;
nstraight = 23*2;

% ncirc = 5;
% nstraight = 11;

rad = 0.16;
len1 = 0.61;
len2 = rad / tan(15*pi/180);
xnose = 0.0;

% [xep, rep, rad, Vex] = setupbody( ncirc, nstraight, rad, len1, len2, xnose );
[xep, rep, rad, Vex] = setupsphere( ncirc, 1.0, 1.0 );
rep(1) = 0;
rep(end) = 0;

name = 'Body';


names{1} = name;
xepts{1} = xep;
repts{1} = rep;
kuttas{1} = kutta;
props{1} = false;
deltaCP{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=Vex;



% Execute script
run('bor')

xsol = xcp;
rsol = rcp;
gammasol = gammas{1};
clear gammas;




ncirc = length( xepts{1} );

% Construct an ellipse as starting point
rada = max( xepts{1} ) * 0.5;
radb = max( repts{1} );
radb = rada * 0.1;
xcen = rada;
%[xep, rep, ~, ~] = setupbody( ncirc, nstraight, 0.32, 0.3, 0.5, xnose );
%[xep, rep, ~, ~] = setupellipsoid( ncirc, rada, radb, xcen );

% reverse engineer matching ellipsoid points from body of revolution
% xep = rada * cos( thetacirc ) + xcen;
thetacirc = acos((xep - xcen)./rada);
rep = radb * sin( thetacirc );


rad = inf( 1, ncirc - 1 );

xcen = ( xep(2:end) + xep(1:end-1) ) * 0.5;

name = 'Designed';


%gammatgt{1} = interp1( xsol{1}./max(xsol{1}), gammasol, xcen./max(xcen), 'linear', 'extrap' );
gammatgt{1} = gammasol;


drawplots = true;

names{1} = name;
xepts{1} = xep;
repts{1} = rep;
kuttas{1} = kutta;
props{1} = false;
deltaCP{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=gammatgt{1};

% Execute script
run('bor_inv')
