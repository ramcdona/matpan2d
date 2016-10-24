clear all
format compact
close all

% Freestream velocity
W = 1.0;

Sref = 1.0;

xlim = [-1 6];
rlim = [0 3];
nsurvey = 51;

nsl = 11;

drawplots = true;


ntstep = 5;

dCP = 5;


% Disk radius
rdisk = 1.0;

% Disk location
xstart = 0;
% End of contracting streamtube
xend = 4;

% Vortex ring spacing
dxring = rdisk * 0.1;

% Number of vortex ring panels
npan = ( xend - xstart ) / dxring;

% Set up initial geometry
xpts = linspace( xstart, xend, npan + 1 );
rpts = rdisk * ones( size( xpts ) );

% Set limits on streamtube dimensions
rmin = zeros( size( rpts ) );
rmax = 1.5 * ones( size( rpts ) );


names{1} = 'Disk';
xepts{1} = xpts;
repts{1} = rpts;
remin{1} = rmin;
remax{1} = rmax;
kuttas{1} = false;
props{1} = true;
deltaCP{1} = dCP;
jtels{1} = 0;
jteus{1} = 0;
rads{1} = 0;
Vexs{1} = 0;

% Execute script
run('bor')

