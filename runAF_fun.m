function [CLi, CDi, xepsave, yepsave] = runAF_fun( xoff2, yoff2, th, drawplots, aftype )

% clear all
% format compact
% close all

% Freestream velocity
W = 1.0;
Minf = 0.0;

Sref = 1.0;

xlim = [-3 4];
ylim = [-3 4];
nsurvey = 51;

nsl = 21;

% drawplots = true;
drawmoreplots = false;

ntstep = 1;

% Construct a duct
naf = 151;  % 151

chord = 1.0;
%th = 5.0;       % Degrees about LE
xoff = 0.0;      % X-Offset of LE
yoff = 0.0;      % R-Offset of LE

th2 = th;
%th2 = 5.0;       % Degrees about LE
%xoff2 = 0.0;      % X-Offset of LE
%yoff2 = 0.5;      % R-Offset of LE

flipaf = false;


switch aftype
    case 1
        % Airfoil point spacing
        % 1 -- Cosine cluster LE, TE
        % 2 -- Cosine cluster LE
        % 3 -- Uniform
        spacing = 1;

        % NACA 4-Digit airfoil parameters
        dig1 = 0;
        dig2 = 0;
        dig34 = 12;
        [xep, yep, rad, Vex] = setupNACAduct( naf, chord, th, xoff, yoff, dig1, dig2, dig34, flipaf, spacing );
        afname=['NACA ' num2str(dig1) num2str(dig2) num2str(dig34) ' Airfoil'];
    case 2
        [xep, yep, rad, Vex, naf] = setupVR12duct( chord, th, xoff, yoff, flipaf );
        afname=['VR12'];
    case 3
        % Von Karman Trefftz airfoil for validation.
        % 0.06 -- epsion (thickness term)
        % 0.05 -- kappa (camber term)
        % 10 deg -- tau (TE angle)
        % Airfoil spacing results from transformation.
        [xep, yep, rad, Vex] = setupVKTduct( naf, chord, th, xoff, yoff, 0.06, 0.05, 10*pi/180, flipaf );
        afname=['VKT'];
end


xepsave = xep;
yepsave = yep;

% Index for Kutta condition
jtelow = 1;
jteup = naf - 1;
kutta = true;


names{1} = afname;
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
    
    switch aftype
        case 1
            % Airfoil point spacing
            % 1 -- Cosine cluster LE, TE
            % 2 -- Cosine cluster LE
            % 3 -- Uniform
            spacing = 1;

            % NACA 4-Digit airfoil parameters
            dig1 = 0;
            dig2 = 0;
            dig34 = 12;
            [xep2, yep2, rad2, Vex2] = setupNACAduct( naf, chord, th2, xoff2, yoff2, dig1, dig2, dig34, flipaf, spacing );
            afname2=['NACA ' num2str(dig1) num2str(dig2) num2str(dig34) ' Airfoil'];
        case 2
            [xep2, yep2, rad2, Vex2, naf] = setupVR12duct( chord, th2, xoff2, yoff2, flipaf );
            afname2=['VR12'];
        case 3
            % Von Karman Trefftz airfoil for validation.
            % 0.06 -- epsion (thickness term)
            % 0.05 -- kappa (camber term)
            % 10 deg -- tau (TE angle)
            % Airfoil spacing results from transformation.
            [xep2, yep2, rad2, Vex2] = setupVKTduct( naf, chord, th2, xoff2, yoff2, 0.06, 0.05, 10*pi/180, flipaf );
            afname2=['VKT'];
    end
    
    names{2} = afname2;
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
    rads{2}=rad2;
    Vexs{2}=Vex2;
    
    [d_min] = p_poly_dist(xep, yep, xep2, yep2, true );
    d_min = min(d_min);
    if ( d_min <= ( 0.25 * 0.12 ) )  % Quarter an airfoil thickness separation.
        CLi = {nan, nan};
        CDi = {nan, nan};
        xep = nan;
        yep = nan;
        return
    end
    
end

% Execute script
run('planar')

