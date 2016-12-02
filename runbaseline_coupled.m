function [S, DLTAST] = runbaseline_coupled( Sin, DLTASTin, drawplots )

% Freestream velocity
W = 1.0;
Minf = 0.785;

Sref = 1.0;

xlim = [-10 140];
rlim = [0 20];
nsurvey = 51;

nsl = 21;

ntstep = 7;

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
[xfusegrid, rfusegrid] = fuse( 10 );

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


% Actuator disk
dCP = 177.889/204.99;

rdisk = repts{2}(1);
xstart = xepts{2}(1);
xdsk = xoff + chord * 0.5;

% End of contracting streamtube
xend = xstart + 10;

% Vortex ring spacing
dxring = rdisk * 0.1;

% Number of vortex ring panels
npan = ( xend - xstart ) / dxring;

% Set up initial geometry
xpts = linspace( xstart, xend, npan + 1 );
rpts = rdisk * ones( size( xpts ) );

names{3} = 'Disk';
xepts{3} = xpts;
repts{3} = rpts;
kuttas{3} = false;
props{3} = true;
deltaCP{3} = dCP;
xdisk{3} = xdsk;
jtels{3} = 0;
jteus{3} = 0;
rads{3} = 0;
Vexs{3} = nan;


% Execute script
run('bor')

run('dimensionalize')

system( './bl2d < lobli_bl2d.in > lobli_bl2d.out' )

[S, DLTAST] = read_bl2d( 'lobli_bl2d.out', drawplots );
