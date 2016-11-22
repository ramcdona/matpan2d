
Tamb = 390.53;  % deg R
Pamb = 475.35;  % lbf/ft^2

Neta = 301;      % Mesh in eta-direction
Ns = 5000;       % Mesh in S-direction
Nprofile = Ns/5; % Number of profiles to output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



R = 1716;
rhoamb = Pamb/(R*Tamb);

spdsnd = sqrt( gammagas * R * Tamb );
Vinf = Minf * spdsnd;

qinf = 0.5 * rhoamb * Vinf.^2;


Pt = Pamb * (1.0 + 0.5 * ( gammagas - 1.0 ) * Minf.^2 ).^( gammagas / ( gammagas - 1.0 ) );
Tt = Tamb * (1.0 + 0.5 * ( gammagas - 1.0 ) * Minf.^2 );

%for iseg=1:nseg
for iseg=1:1

    P = ( CpNLR{iseg}.* qinf + Pamb ) ./ Pt;
    P = [ 1.0 P(1:end-1) 1.0 ];

    x = [xepts{iseg}(1) xcp{iseg}(1:end-1) xepts{iseg}(end)];
    r = [repts{iseg}(1) rcp{iseg}(1:end-1) repts{iseg}(end)] ./ beta;

    dxs = x(2:end) - x(1:end-1);
    drs = r(2:end) - r(1:end-1);

    dss = sqrt( dxs.^2 + drs.^2 );

    s = cumsum(dss);
    s = [0 s];

end

Nbody = length(P);


fid = fopen( 'lobli_bl2d.in', 'w' );


fprintf( fid, '  BL2D Input from LOBLI\n' );
fprintf( fid, '  $NAM1\n\n' );

fprintf( fid, '  XMA     = %f\n', Minf );
fprintf( fid, '  PT1     = %f\n', Pt );
fprintf( fid, '  TT1     = %f\n', Tt );

fprintf( fid, '  PRNTINC = %.16E\n', s(end)/Nprofile );
fprintf( fid, '  PROINC  = %.16E\n', s(end)/Nprofile );


fprintf( fid, '  JI      = 1\n' );
fprintf( fid, '  JK      = %d\n', Neta + 2 );
fprintf( fid, '  JL      = 1\n' );
fprintf( fid, '  JM      = %d\n', Nprofile );
fprintf( fid, '  JN      = %d\n', Nprofile );
fprintf( fid, '  JH      = %d\n', Ns + 1 );
fprintf( fid, '  JJ      = %d\n', Nbody );

fprintf( fid, '  IE      = %d\n', Neta );
fprintf( fid, '  IEND1   = %d\n', Ns );


fprintf( fid, '  CONV   = 0.1E-06\n' );
fprintf( fid, '  CONVE  = 0.1E-01\n' );
fprintf( fid, '  DETA1  = 0.E+00\n' );
fprintf( fid, '  FT     = 0.1E+01\n' );
fprintf( fid, '  G      = 0.14E+01\n' );
fprintf( fid, '  GLAR   = 0.E+00\n' );
fprintf( fid, '  IBODY  =           1\n' );
fprintf( fid, '  IENTRO =           1\n' );
fprintf( fid, '  IGAS   =           1\n' );
fprintf( fid, '  IGEOM  =           1\n' );
fprintf( fid, '  IPRNT  =           0\n' );
fprintf( fid, '  IPRO   =           0\n' );
fprintf( fid, '  ITMAX  =           3\n' );
fprintf( fid, '  IYINT  =           1\n' );
fprintf( fid, '  J      =           1\n' );
fprintf( fid, '  KODAMP =           2\n' );
fprintf( fid, '  KODE   =           0\n' );
fprintf( fid, '  KODPRT =           1\n' );
fprintf( fid, '  KODUNIT=           0\n' );
fprintf( fid, '  KODVIS =           2\n' );
fprintf( fid, '  KODWAL =           2\n' );
fprintf( fid, '  KTCOD  =           2\n' );
fprintf( fid, '  NAUXPRO=           2\n' );
fprintf( fid, '  NITMAX =          30\n' );
fprintf( fid, '  NUMB1  =           0\n' );
fprintf( fid, '  PHII   = 0.9E+02\n' );
fprintf( fid, '  PR     = 0.7E+00\n' );
fprintf( fid, '  PRNTVAL= 0.E+00\n' );
fprintf( fid, '  PROVAL = 0.E+00\n' );
fprintf( fid, '  PRT    = 0.95E+00\n' );
fprintf( fid, '  PRTAR  = 0.E+00\n' );
fprintf( fid, '  R      = 0.1716E+04\n' );
fprintf( fid, '  SMXTR  = 0.1E+09\n' );
fprintf( fid, '  SST    = 0.E+00\n' );
fprintf( fid, '  TLNGTH = 0.2E+01\n' );
fprintf( fid, '  VELEDG = 0.99995E+00\n' );
fprintf( fid, '  VIS1C1 = 0.227E-07\n' );
fprintf( fid, '  VIS1C2 = 0.1986E+03\n' );
fprintf( fid, '  VIS2C1 = 0.E+00\n' );
fprintf( fid, '  VIS2C2 = 0.E+00\n' );
fprintf( fid, '  W      =           1\n' );
fprintf( fid, '  WAVE   = 0.E+00\n' );
fprintf( fid, '  XEND   = 0.3E+03\n' );
fprintf( fid, '  XK     = 0.104E+01\n' );
fprintf( fid, '  XT1    = 0.4E+00\n' );
fprintf( fid, '  XT2    = 0.168E-01\n' );
fprintf( fid, '  XT3    = 0.5E+01\n' );
fprintf( fid, '  XT4    = 0.78E+00\n' );
fprintf( fid, '  XT5    = 0.108E+00\n' );
fprintf( fid, '  XT6    = 0.26E+02\n' );
fprintf( fid, '  NZT    =         100\n' );
fprintf( fid, '  IOX    =           1\n' );

fprintf( fid, '/\n' );

fprintf( fid, '  $NAM2\n' );
fprintf( fid, '  L      =           1\n' );
fprintf( fid, '  SS     =        -1.0\n' );

fprintf( fid, '  NUMBER = %d\n', Nbody );
fprintf( fid, '  QW     = %d*0.,\n', Nbody );
fprintf( fid, '  RVWALD = %d*0.,\n', Nbody );
fprintf( fid, '  RADC   = %d*1.E31,\n', Nbody );

namelist_vec( fid, 'PE', P );

namelist_vec( fid, 'RMI', r );

namelist_vec( fid, 'S', s );

namelist_vec( fid, 'Z', x );

fprintf( fid, '/\n' );
fclose( fid );

