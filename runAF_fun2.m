function [CLi, CDi, xep, yep ] = runAF_fun2( xoff2, yoff2, th2, drawplots )

%clear all
%format compact
%close all

% Freestream velocity
% W = 1.0;
Wfield = 1.0;
Minf = 0.0;

Sref = 1.0;

xlim = [-2 2];
ylim = [-2 2];
nsurvey = 51;

nsl = 21;

% drawplots = true;
drawmoreplots = false;

ntstep = 1;

% Construct a duct
naf = 151;

chord = 1.0;
th = 0;       % Degrees about LE
xoff = 0.0;      % X-Offset of LE
yoff = 0.0;      % R-Offset of LE

% th2 = 10;       % Degrees about LE
% xoff2 = 0.0;      % X-Offset of LE
% yoff2 = 0.75;      % R-Offset of LE


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


% [xep, yep, rad, Vex] = setupNACAduct( naf, chord, th, xoff, yoff, dig1, dig2, dig34, flipaf, spacing );
[xep, yep, rad, Vex] = setupellipse( naf, 0.5, 0.5, xoff, yoff );

% Index for Kutta condition
jtelow = 1;
jteup = naf - 1;
kutta = true;


names{1} = ['NACA ' num2str(dig1) num2str(dig2) num2str(dig34) ' Duct'];
W{1} = 0.0;
xuppts{1} = [];
yuppts{1} = [];
xlowpts{1} = [];
ylowpts{1} = [];
xepts{1} = xep;
yepts{1} = yep;
kuttas{1} = kutta;
props{1} = false;
deltaCP{1} = 0.0;
jtels{1} = jtelow;
jteus{1} = jteup;
rads{1}=rad;
Vexs{1}=Vex;

if( true )
    [xep2, yep2, rad, Vex] = setupNACAduct( naf, chord, th2, xoff2, yoff2, dig1, dig2, dig34, flipaf, spacing );
%    [xep2, yep2, rad2, Vex2] = setupellipse( naf, 0.5, 0.5, xoff2, yoff2 );

    names{2} = ['NACA ' num2str(dig1) num2str(dig2) num2str(dig34) ' Duct'];
    W{2} = 1.0;
    xuppts{2} = [];
    yuppts{2} = [];
    xlowpts{2} = [];
    ylowpts{2} = [];
    xepts{2} = xep2;
    yepts{2} = yep2;
    kuttas{2} = kutta;
    props{2} = false;
    deltaCP{2} = 0.0;
    jtels{2} = jtelow;
    jteus{2} = jteup;
    rads{2}=rad;
    Vexs{2}=Vex;
end


rep2s = sqrt( xep2.^2 + yep2.^2 ) - 0.5;

if ( min(rep2s) < 0.05 )
    
    CLi{1} = nan;
    CLi{2} = nan;
    CDi = CLi;
    
    return;
end


% Execute script
run('planar')

