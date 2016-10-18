clear all
format compact
close all

% Freestream velocity
W = 1.0;

Sref = 1.0;

% Velocity survey vectors
xgrd = linspace( -1, 6, 71 );
rgrd = linspace( 0, 3, 31 );
rgrd = rgrd + (pi-3)/1000;  % Offset by small non-round number.

nsl = 11;
rsl = linspace( 0, rgrd(end), nsl);
rsl(1) = 1e-3;
xsl = min(xgrd) * ones( size(rsl) ) + 0.01;
streamback = false;

drawplots = true;


ntstep = 5;

CT = 5;

% Calculate vortex tube strength and eventual jet velocity
gammainf = -W * ( sqrt( CT + 1 ) - 1 );
Wjinf = W - gammainf;


% Disk radius
rdisk = 1.0;

% Disk location
xstart = 2;
% End of contracting streamtube
xend = 4;

% Vortex ring spacing
dxring = rdisk * 0.15;

% Number of vortex ring panels
npan = ( xend - xstart ) / dxring;

% Set up initial geometry
xpts = linspace( xstart, xend, npan + 1 );
rpts = rdisk * ones( size( xpts ) );

% Set limits on streamtube dimensions
rmin = zeros( size( rpts ) );
rmax = 1.5 * ones( size( rpts ) );

gammadisk = gammainf * ones( size( xpts ) );



names{1} = 'Disk';
xepts{1} = xpts;
repts{1} = rpts;
remin{1} = rmin;
remax{1} = rmax;
kuttas{1} = false;
props{1} = true;
gammaad{1} = gammadisk;
Wjad{1} = Wjinf;
jtels{1} = 0;
jteus{1} = 0;
rads{1} = 0;
Vexs{1} = 0;

% Execute script
run('bor')

