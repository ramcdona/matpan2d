clear all
format compact

% Freestream velocity
W = 1.0;

Sref = 1.0;

% Velocity survey vectors
xgrd = linspace( -1, 2, 131 );
rgrd = linspace( 0, 2, 121 );
rgrd = rgrd + (pi-3)/1000;  % Offset by small non-round number.

nsl = 11;
rsl = linspace( 0, rgrd(end), nsl);
rsl(1) = 1e-3;
xsl = min(xgrd) * ones( size(rsl) ) + 0.01;
streamback = false;

drawplots = true;

names{1} = 'Disk';
xepts{1} = 1.0 * [1 1];
repts{1} = 0.5 * [0 1];
kuttas{1} = false;
props{1} = true;
gammaad{1} = 5.0;
jtels{1} = 0;
jteus{1} = 0;
rads{1} = 0;
Vexs{1} = 0;

% Execute script
run('bor')

