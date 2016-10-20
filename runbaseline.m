clear all
format compact

% Freestream velocity
W = 1.0;

Sref = 1.0;

% Velocity survey vectors
xgrd = unique( [-10:2:-2 -2:.1:2 2:2:120 120:.1:135 135:2:140]);
rgrd = unique( [0:.1:1.5 1.5:.25:2.5 2.5:.1:4 4:.5:10 10:1:20] );
rgrd = rgrd + (pi-3)/1000;  % Offset by small non-round number.

nsl = 21;
rsl = linspace( 0, rgrd(end), nsl) + 0.3;
xsl = min(xgrd) * ones( size(rsl) ) + 1;
streamback = false;

drawplots = true;

ntstep = 1;

% Construct fuselage mesh

% Construct a duct
naf = 101;

chord = 7.2;
alpha = 0;       % Degrees about LE
xoff = 128.35-chord/2;   % X-Offset of LE
roff = 3.0;      % R-Offset of LE

% NACA 4-Digit airfoil parameters
dig1 = 3;
dig2 = 2.5;
dig34 = 10;

flipaf = false;

% Airfoil point spacing
% 1 -- Cosine cluster LE, TE
% 2 -- Cosine cluster LE
% 3 -- Uniform
spacing = 1;


% Actually generate bodies
[xfusegrid, rfusegrid] = fuse();

jtelow = -1;
jteup = - 1;
kutta = false;

names{1} = 'Fuselage';
xepts{1} = xfusegrid;
repts{1} = rfusegrid;
kuttas{1} = kutta;
props{1} = false;
deltaCP{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1} = inf(1,length(xfusegrid)-1);
Vexs{1} = nan(1,length(xfusegrid)-1);


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

