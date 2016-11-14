
Tamb = 390.53;  % deg R
Pamb = 475.35;  % lbf/ft^2

R = 1716;

rhoamb = Pamb/(R*Tamb);

spdsnd = sqrt( gammagas * R * Tamb );
Vinf = Minf * spdsnd;

qinf = 0.5 * rhoamb * Vinf.^2;


Pt = Pamb + qinf;
Tt = Tamb * (1.0 + 0.5 * ( gammagas - 1.0 ) * Minf.^2 );

%for iseg=1:nseg
for iseg=1:1

    P = Cp{iseg} .* qinf + Pamb;

    P = [Pt P(1:end-1) Pt];

    x = [xepts{iseg}(1) xcp{iseg}(1:end-1) xepts{iseg}(end)];
    r = [repts{iseg}(1) rcp{iseg}(1:end-1) repts{iseg}(end)];

    dxs = x(2:end) - x(1:end-1);
    drs = r(2:end) - r(1:end-1);

    dss = sqrt( dxs.^2 + drs.^2 );

    s = cumsum(dss);
    s = [0 s];

end

fid = fopen( 'bornam2.in', 'w' );

namelist_vec( fid, 'PE', P./Pt );

namelist_vec( fid, 'RMI', r );

namelist_vec( fid, 'S', s );

namelist_vec( fid, 'Z', x );

fclose( fid );

