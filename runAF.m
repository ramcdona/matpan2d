clear all
format compact
close all

% Freestream velocity
W = 1.0;
Minf = 0.0;

Sref = 1.0;

xlim = [-2 2];
ylim = [-2 2];
nsurvey = 51;

nsl = 21;

drawplots = true;

ntstep = 1;

% Construct a duct
naf = 151;

chord = 1.0;
th = 5;       % Degrees about LE
xoff = 0.0;      % X-Offset of LE
yoff = 0.0;      % R-Offset of LE

% NACA 4-Digit airfoil parameters
dig1 = 0;
dig2 = 0;
dig34 = 12;

flipaf = false;

% Airfoil point spacing
% 1 -- Cosine cluster LE, TE
% 2 -- Cosine cluster LE
% 3 -- Uniform
spacing = 1;


[xep, yep, rad, Vex] = setupNACAduct( naf, chord, th, xoff, yoff, dig1, dig2, dig34, flipaf, spacing );

% Index for Kutta condition
jtelow = 1;
jteup = naf - 1;
kutta = true;


names{1} = ['NACA ' num2str(dig1) num2str(dig2) num2str(dig34) ' Duct'];
xepts{1} = xep;
yepts{1} = yep;
kuttas{1} = kutta;
props{1} = false;
deltaCP{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=Vex;

% Execute script
run('planar')

