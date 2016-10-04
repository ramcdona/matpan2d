clear all
format compact

% Freestream velocity
W = 1.0;

drawplots = true;

% 1 -- Sphere
% 2 -- BOR with aft cone
% 3 -- Ellipsoid
% 4 -- NACA 4-Digit duct
runcase = 3;

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

elseif (runcase == 2)

    % Construct example body
    ncirc = 10;
    nstraight = 22;
    rad = 0.16;
    len1 = 0.61;
    len2 = rad / tan(15*pi/180);
    xnose = 0.0;

    [xep, rep, rad, Vex] = setupbody( ncirc, nstraight, rad, len1, len2, xnose );

elseif (runcase == 3 )

    % Construct an ellipse
    ncirc = 101;
    rada = 5.0;
    radb = 0.5;
    xcen = 0.0;

    [xep, rep, rad, Vex] = setupellipsoid( ncirc, rada, radb, xcen );

else
    naf = 51;

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

    [xep, rep, rad, Vex] = setupNACAduct( naf, alpha, xoff, roff, dig1, dig2, dig34, flipaf, spacing );

    % Index for Kutta condition
    jtelow = 1;
    jteup = naf - 1;
    kutta = true;

end

xepts{1} = xep;
repts{1} = rep;
kuttas{1} = kutta;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=Vex;

% Execute script
run('bor')

