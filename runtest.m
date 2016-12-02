clear all
format compact

[S, DLTAST] = runbaseline_coupled( [], [], false );

% Hold final value constant
S = [0 S 2*S(end)];
DLTAST = [0 DLTAST DLTAST(end)];

[S2, DLTAST2] = runbaseline_coupled( S, DLTAST, false );

S2 = [0 S2 2*S2(end)];
DLTAST2 = [0 DLTAST2 DLTAST2(end)];

[S3, DLTAST3] = runbaseline_coupled( S2, DLTAST2, false );

S3 = [0 S3 2*S3(end)];
DLTAST3 = [0 DLTAST3 DLTAST3(end)];

[S4, DLTAST4] = runbaseline_coupled( S3, DLTAST3, true );

S4 = [0 S4 2*S4(end)];
DLTAST4 = [0 DLTAST4 DLTAST4(end)];

n1 = length(DLTAST);
n2 = length(DLTAST2);
n3 = length(DLTAST3);
n4 = length(DLTAST4);

figure(123)
semilogy(S(1:min(n1,n2)),(DLTAST(1:min(n1,n2))-DLTAST2(1:min(n1,n2))),S(1:min(n1,n2)),-(DLTAST(1:min(n1,n2))-DLTAST2(1:min(n1,n2))),...
    S(1:min(n2,n3)),(DLTAST2(1:min(n2,n3))-DLTAST3(1:min(n2,n3))),S(1:min(n2,n3)),-(DLTAST2(1:min(n2,n3))-DLTAST3(1:min(n2,n3))),...
    S(1:min(n3,n4)),(DLTAST3(1:min(n3,n4))-DLTAST4(1:min(n3,n4))),S(1:min(n3,n4)),-(DLTAST3(1:min(n3,n4))-DLTAST4(1:min(n3,n4))) )
