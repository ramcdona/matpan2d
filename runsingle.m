clear all
format compact
close all

% Freestream velocity
W = 1.0;

Sref = 1.0;

xlim = [-1.5 7];
rlim = [0 3];
nsurvey = 51;

nsl = 23;

drawplots = true;

ntstep = 7;

% 1 -- Sphere
% 2 -- BOR with aft cone
% 3 -- Ellipsoid
% 4 -- NACA 4-Digit duct
% 5 -- Bontempo body
runcase = 4;

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
    ncirc = 1501;
    rada = 1.1;
    radb = 0.374;
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
deltaCP{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=Vex;


dCP = .75;

% Disk radius
rdisk = 1.0;
% rdisk = 0.5;
rdisk = 0.46;
rdisk = 0.486144640095513 - 0.05;  % 'on' duct lower surface.
% rdisk = 0.4128442572534;

% Disk location
% xstart = 1;
% xstart = 1.0;
xstart = 0.5;
%xstart = 0.99619469809175;
% xstart = 0.0;

% End of contracting streamtube
xend = 5;

% Vortex ring spacing
dxring = rdisk * 0.1;
dxring = .05;

% Number of vortex ring panels
npan = ( xend - xstart ) / dxring;

% Set up initial geometry
xpts = linspace( xstart, xend, npan + 1 );
rpts = rdisk * ones( size( xpts ) );

names{2} = 'Disk';
xepts{2} = xpts;
repts{2} = rpts;
kuttas{2} = false;
props{2} = true;
deltaCP{2} = dCP;
xdisk{2} = 0.5
jtels{2} = 0;
jteus{2} = 0;
rads{2} = 0;
Vexs{2} = 0;

% Execute script
run('bor')

