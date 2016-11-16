clear all
format compact


%fid = fopen( 'baseline.out', 'r' );

fid = fopen( 'lobli_bl2d.out', 'r' );


MAXSTA = 1000;

fposchar( fid, '=', 3 );
IE = fscanf( fid, '%f' );



ETA = nan( IE, MAXSTA );
Y_YE = ETA;
Y_dim = ETA;
U_UE = ETA;
U_dim = ETA;
T_TE = ETA;
TT_TTE = ETA;
CROCCO = ETA;
PT_PTR = ETA;
M_ME = ETA;
YPLUS = ETA;
UPLUS = ETA;
UDEF = ETA;
VISEFF = ETA;

next = true;

ista = 0;
while( next )

    next = fposstr( fid, 'Y/YE' );

    if( next )
        ista = ista + 1;

        M = fscanf( fid, '%f', [12 IE] );

        len = size( M, 2 );

        ETA(1:len,ista) = M(1,:);
        Y_YE(1:len,ista) = M(2,:);
        U_UE(1:len,ista) = M(3,:);
        T_TE(1:len,ista) = M(4,:);
        TT_TTE(1:len,ista) = M(5,:);
        CROCCO(1:len,ista) = M(6,:);
        PT_PTR(1:len,ista) = M(7,:);
        M_ME(1:len,ista) = M(8,:);
        YPLUS(1:len,ista) = M(9,:);
        UPLUS(1:len,ista) = M(10,:);
        UDEF(1:len,ista) = M(11,:);
        VISEFF(1:len,ista) = M(12,:);

        % Position file at S
        fposchar( fid, '=', 2 );
        S(ista) = fscanf( fid, '%f' );

        fposchar( fid, '=', 3 );
        CFW(ista) = fscanf( fid, '%f' );

        fposchar( fid, '=', 2 );
        YE(ista) = fscanf( fid, '%f' );

        Y_dim(1:len,ista) = YE(ista) * M(2,:);

        jnk = fgetl( fid );
        fposchar( fid, '=', 1 );
        RMI(ista) = fscanf( fid, '%f' );

        fposchar( fid, '=', 1 );
        PE(ista) = fscanf( fid, '%f' );

        jnk = fgetl( fid );
        fposchar( fid, '=', 1 );
        Z(ista) = fscanf( fid, '%f' );

        fposchar( fid, '=', 2 );
        DLTAST(ista) = fscanf( fid, '%f' );

        jnk = fgetl( fid ); % Finish line
        fposchar( fid, '=', 3 );
        THETA(ista) = fscanf( fid, '%f' );

        jnk = fgetl( fid ); % Finish line
        fposchar( fid, '=', 2 );
        UE(ista) = fscanf( fid, '%f' );

        U_dim(1:len,ista) = UE(ista) * M(3,:);

        jnk = fgetl( fid ); % Finish line
        fposchar( fid, '=', 3 );
        TAUD(ista) = fscanf( fid, '%f' );

        jnk = fgetl( fid ); % Finish line
        fposchar( fid, '=', 3 );
        CFE(ista) = fscanf( fid, '%f' );

    end
end
fclose( fid );

NSTA = ista;

figure(1)
plot( U_UE, Y_dim );

figure(2)
plot( U_UE, Y_YE );

figure(3)
plot( U_dim, Y_dim );

figure(4)
plot( UPLUS, YPLUS );



figure(5)
plot( S, CFE, S, CFW )
xlabel('S (ft)')
legend('c_{f,e}','c_{f,w}')

figure(6)
plot( S, TAUD )
xlabel('S (ft)')
ylabel('\tau_w (lbf/ft^2)')

figure(7)
plot( S, DLTAST, S, THETA, S, YE )
xlabel('S (ft)')
ylabel('BL Thickness (ft)')
legend( '\delta^* (Displacement)', '\theta (Momentum)', 'Y_e' )



[xfusegrid, rfusegrid] = fuse( 1 );

figure(8)
plot( xfusegrid, rfusegrid, 'Color', [0.5 0.5 0.5] )
hold on

kv = .001;

for ista = 1:NSTA
    plot( Z(ista) + kv*U_dim(:,ista), RMI(ista) + Y_dim(:,ista),'k', 'LineWidth', 1.0 )
end
hold off
axis equal
ax=axis;
ax(1) = 0;
ax(3) = 0;
ax(4) = 10;
axis(ax);
axis off


figure(9)
plot( xfusegrid, rfusegrid, 'k' )
hold on
plot( [0 Z] , [0 RMI + DLTAST], 'Color', [0.5 0.5 0.5] )
plot( [0 Z] , [0 RMI + THETA], '--', 'Color', [0.5 0.5 0.5] )
plot( [0 Z] , [0 RMI + YE], 'k:' )
hold off
axis equal
ax=axis;
ax(1) = 0;
ax(3) = 0;
ax(4) = 10;
axis(ax);
axis off



figure(10)
plot( xfusegrid, rfusegrid, 'k' )
hold on
plot( [0 Z] , [0 RMI + DLTAST], 'Color', [0.5 0.5 0.5] )
plot( [0 Z] , [0 RMI + THETA], '--', 'Color', [0.5 0.5 0.5] )
plot( [0 Z] , [0 RMI + YE], 'k:' )

for ista = 1:NSTA
    plot( Z(ista) + kv*U_dim(:,ista), RMI(ista) + Y_dim(:,ista),'k', 'LineWidth', 1.0 )
end
hold off
axis equal
ax=axis;
ax(1) = 0;
ax(3) = 0;
ax(4) = 10;
axis(ax);
axis off


figure(11)
plot( xfusegrid, rfusegrid, 'k' )
hold on
plot( Z , PE )
