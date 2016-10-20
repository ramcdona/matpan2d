function [xep, rep, rad, Vex] = setupbontempobody( npts )

M = [-1.0482	0
    -1.046	0.0264
    -1.0395	0.0557
    -1.0287	0.0829
    -1.0136	0.1096
    -0.9944	0.1354
    -0.9709	0.1602
    -0.9435	0.1837
    -0.9121	0.2057
    -0.877	0.226
    -0.8381	0.2444
    -0.7958	0.2607
    -0.7501	0.2745
    -0.7012	0.2861
    -0.6494	0.2974
    -0.5948	0.3088
    -0.5376	0.3199
    -0.4781	0.3305
    -0.4165	0.3404
    -0.3531	0.3495
    -0.2881	0.3575
    -0.2217	0.3641
    -0.1543	0.3692
    -0.0861	0.3726
    -0.0173	0.374
    0.0518	0.3735
    0.1209	0.3707
    0.1897	0.3658
    0.2579	0.3585
    0.3254	0.349
    0.3917	0.3372
    0.4568	0.3232
    0.5202	0.3072
    0.5818	0.2893
    0.6412	0.2698
    0.6984	0.2489
    0.753	0.2269
    0.8048	0.2041
    0.8537	0.1808
    0.8994	0.1573
    0.9417	0.1343
    0.9806	0.112
    1.0158	0.0909
    1.0471	0.0713
    1.0746	0.0535
    1.098	0.0378
    1.1173	0.0246
    1.1323	0.014
    1.1432	0.0063
    1.1497	0.0016
    1.1518	0];

xbody = M(:,1)';
rbody = M(:,2)';


% dx = xbody(2:end)-xbody(1:end-1);
% dr = rbody(2:end)-rbody(1:end-1);
%
% ds = sqrt( dx.^2 + dr.^2 );
% s = [0 cumsum(ds) ];
%
% sint = linspace(0, s(end), npts );
%
% xep = interp1( s, xbody, sint, 'spline' );
% rep = interp1( s, rbody, sint, 'spline' );

xep = xbody;
rep = rbody;

rad = inf( 1, length( xep ) - 1 );
Vex = nan( 1, length( xep ) - 1 );

end