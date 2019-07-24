clear all
format compact
close all

% Freestream velocity
W = 1.0;
Minf = 0.0;

Sref = 1.0;

xlim = [-3 12];
ylim = [-4 4];
nsurvey = 51;

nsl = 41;

drawplots = true;


ntstep = 5;

dCP = 5;


% Disk radius
rad = 1.0;

% Disk normal (thrust direction)
% ndisk = [-1; 0]; % Propeller
% ndisk = [0; 1]; % Lift rotor
ndisk = [-1; 1];
% Ensure unit vector.
ndisk = ndisk ./ norm(ndisk);

% Streamtube eventually aligns with freestream.
if ( W == 0 ) % No freestream, 
    tubedir = ndisk;
else
    tubedir = [-1;0];
end
tubedir = tubedir ./ norm(tubedir);


% Unit direction from normal to lower lip of actuator disk.
lowdisk = [0 -1; 1 0] * ndisk;
lowtube = [0 -1; 1 0] * tubedir;


% Disk location
xcen = 0;
ycen = 0;

% Lower and upper disk endpoints.
xl0 = xcen + lowdisk(1) * rad;
yl0 = ycen + lowdisk(2) * rad;

xu0 = xcen - lowdisk(1) * rad;
yu0 = ycen - lowdisk(2) * rad;

dtube = [10, -1];

% % Length of contracting streamtube
% tubelen = 10;
% 
% 
% % Initial streamtube endpoints.
% xl1 = xcen - tubedir(1) * tubelen;
% yl1 = ycen - tubedir(2) * tubelen - rad ;
% 
% xu1 = xcen - tubedir(1) * tubelen;
% yu1 = ycen - tubedir(2) * tubelen + rad ;

% Vortex ring spacing
dlring = rad * 0.2;

tubelen = norm(dtube);
% Number of vortex ring panels
npan = round( tubelen / dlring );

quadbez = @(t,pt) kron((1.0-t).^2,pt(1)) + kron(2.0*(1.0-t).*t,pt(2)) + kron(t.^2,pt(3));
ts = linspace(0,1,npan+1);

xmax0 = max(xu0,xl0);

% Set up initial streamtube geometry
names{1} = 'Disk';
xuppts{1} = quadbez(ts,[xu0 (xu0 + xmax0+dtube(1))*0.5 xmax0 + dtube(1)]);
yuppts{1} = quadbez(ts,[yu0 yu0 + dtube(2)     yu0 + dtube(2)]);
xlowpts{1} = quadbez(ts,[xl0 (xl0 + xmax0+dtube(1))*0.5 xmax0 + dtube(1)]);
ylowpts{1} = quadbez(ts,[yl0 yl0 + dtube(2)     yl0 + dtube(2)]);
xepts{1} = [];
yepts{1} = [];
kuttas{1} = false;
props{1} = true;
deltaCP{1} = dCP;
jtels{1} = 0;
jteus{1} = 0;
rads{1} = 0;
Vexs{1} = nan;

% Execute script
run('planar')

