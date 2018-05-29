clear all
format compact
close all

% Freestream velocity
W = 1.0;
Minf = 0.0;

Sref = 1.0;

xlim = [-1 6];
ylim = [-3 3];
nsurvey = 51;

nsl = 11;

drawplots = true;


ntstep = 5;

dCP = 5;


% Disk radius
rad = 1.0;

% Disk normal (thrust direction)
ndisk = [-1; 0];
% Ensure unit vector.
ndisk = ndisk ./ norm(ndisk);

% Unit direction from normal to lower lip of actuator disk.
lowdisk = [0 -1; 1 0] * ndisk;

% Disk location
xcen = 0;
ycen = 0;

% Lower and upper disk endpoints.
xl0 = xcen + lowdisk(1) * rad;
yl0 = ycen + lowdisk(2) * rad;

xu0 = xcen - lowdisk(1) * rad;
yu0 = ycen - lowdisk(2) * rad;


% Length of contracting streamtube
tubelen = 4;

% Initial streamtube endpoints.
xl1 = xl0 - ndisk(1) * tubelen;
yl1 = yl0 - ndisk(2) * tubelen;

xu1 = xu0 - ndisk(1) * tubelen;
yu1 = yu0 - ndisk(2) * tubelen;

% Vortex ring spacing
dlring = rad * 0.1;

% Number of vortex ring panels
npan = round( tubelen / dlring );

% Set up initial streamtube geometry
names{1} = 'Disk';
xepts{1} = [linspace( xl0, xl1, npan + 1 ) linspace( xu0, xu1, npan + 1 )];
yepts{1} = [linspace( yl0, yl1, npan + 1 ) linspace( yu0, yu1, npan + 1 )];
kuttas{1} = false;
props{1} = true;
deltaCP{1} = dCP;
jtels{1} = 0;
jteus{1} = 0;
rads{1} = 0;
Vexs{1} = nan;

% Execute script
run('planar')

