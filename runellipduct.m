clear all
format compact
close all

% Freestream velocity
W = 1.0;

Sref = 1.0;

xlim = [-2 1];
rlim = [0 2];
nsurvey = 51;

nsl = 21;

drawplots = true;

ntstep = 1;

% Construct an ellipse
ncirc = 36;
rada = 1.0;
radb = .25;
xcen = 0.0;

% Construct a duct
naf = 51;

chord = 1.0;
alpha = 5;       % Degrees about LE
xoff = 0.75;      % X-Offset of LE
roff = 0.5;      % R-Offset of LE

% NACA 4-Digit airfoil parameters
dig1 = 4;
dig2 = 4;
dig34 = 12;

flipaf = false;

% Airfoil point spacing
% 1 -- Cosine cluster LE, TE
% 2 -- Cosine cluster LE
% 3 -- Uniform
spacing = 1;


% Actually generate bodies
[xep, rep, rad, Vex] = setupellipsoid( ncirc, rada, radb, xcen );

jtelow = -1;
jteup = - 1;
kutta = false;

names{1} = 'Ellipsoid';
xepts{1} = xep;
repts{1} = rep;
kuttas{1} = kutta;
props{1} = false;
deltaCP{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=Vex;


[xep, rep, rad, Vex] = setupNACAduct( naf, chord, alpha, xoff, roff, dig1, dig2, dig34, flipaf, spacing );

% Index for Kutta condition
jtelow = 1;
jteup = naf - 1;
kutta = true;


names{2} = ['NACA ' num2str(dig1) num2str(dig2) num2str(dig34) ' Duct'];
xepts{2} = xep;
repts{2} = rep;
kuttas{2} = kutta;
props{2} = false;
deltaCP{2} = 0.0;
jtels{2} = jtelow;
jteus{2} = jteup;
rads{2}=rad;
Vexs{2}=Vex;

% Execute script
run('bor')

