clear all
format compact

nsurvey = 1001;

xepts = [ 1.0 1.1 ];
repts = [ 1.0 1.0 ];

xsurvey = linspace(0.5,1.5,nsurvey);
rsurvey = 1.001 * ones(size(xsurvey));

% rsurvey = linspace(0.01, 3, nsurvey);
% xsurvey = 1.049 * ones(size(xsurvey));


uv = zeros( size(xsurvey) );
vv = uv;
uvs = uv;
vvs = uv;

for i=1:length(xsurvey)
    [ui, vi] = ringsubvortex( xsurvey(i), rsurvey(i), xepts(1), repts(1), xepts(2), repts(2) );
    
    uvs(i) = ui;
    vvs(i) = vi;
    
    [ui, vi] = ringvortex( (xepts(1) + xepts(2) ) *0.5, (repts(1) + repts(2) ) * 0.5, xsurvey(i), rsurvey(i) );

    uv(i) = ui;
    vv(i) = vi;
end

figure(1)
plot( xsurvey, uv, '--', xsurvey, uvs )

figure(2)
plot( xsurvey, vv, '--', xsurvey, vvs )


figure(3)
plot( rsurvey, uv, '--', rsurvey, uvs )

figure(4)
plot( rsurvey, vv, '--', rsurvey, vvs )
