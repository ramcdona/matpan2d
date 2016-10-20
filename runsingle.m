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

ntstep = 1;

% 1 -- Sphere
% 2 -- BOR with aft cone
% 3 -- Ellipsoid
% 4 -- NACA 4-Digit duct
% 5 -- Bontempo body
runcase = 4;
runcase = 1;

% Default to no kutta condition, turn on inside runcase setup.
kutta = false;
jtelow = -1;
jteup = -1;
if ( runcase == 1 )
    % Construct a circle
    ncirc = 36;
    rad = 0.5;
    xcen = 0.0;

    [xep, rep, rad, Vex] = setupsphere( ncirc, rad, xcen );
    name = 'Sphere';

elseif (runcase == 2)

    % Construct example body
    ncirc = 10;
    nstraight = 22;
    rad = 0.16;
    len1 = 0.61;
    len2 = rad / tan(15*pi/180);
    xnose = 0.0;

    [xep, rep, rad, Vex] = setupbody( ncirc, nstraight, rad, len1, len2, xnose );
    name = 'Body';

elseif (runcase == 3 )

    % Construct an ellipse
    ncirc = 101;
    rada = 5.0;
    radb = 0.5;
    xcen = 0.0;

    [xep, rep, rad, Vex] = setupellipsoid( ncirc, rada, radb, xcen );
    name = 'Ellipsoid';

elseif (runcase == 4 )
    naf = 51;

    chord = 1.0;
    alpha = 5;       % Degrees about LE
    xoff = 0.0;      % X-Offset of LE
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

    [xep, rep, rad, Vex] = setupNACAduct( naf, chord, alpha, xoff, roff, dig1, dig2, dig34, flipaf, spacing );

    % Index for Kutta condition
    jtelow = 1;
    jteup = naf - 1;
    kutta = true;
    name = ['NACA ' num2str(dig1) num2str(dig2) num2str(dig34) ' Duct'];

else
    npts = 151;

    [xep, rep, rad, Vex] = setupbontempobody( npts );
    name = 'Bontempo Body';
end

names{1} = name;
xepts{1} = xep;
repts{1} = rep;
kuttas{1} = kutta;
props{1} = false;
gammaad{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=Vex;

% Execute script
run('bor')

