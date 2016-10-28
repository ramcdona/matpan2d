

xorig = xepts;
rorig = repts;

nseg = length( xepts );

figure(1)
plot(xepts{1},repts{1},'k')
hold on

for iter=1:50

    npan = 0;
    % Process geometry
    for iseg=1:nseg
        xep = xepts{iseg};
        rep = repts{iseg};

        dx{iseg} = xep(2:end) - xep(1:end-1);
        dr{iseg} = rep(2:end) - rep(1:end-1);

        ds{iseg} = sqrt( dx{iseg}.^2 + dr{iseg}.^2 );        % Panel arclength
        theta{iseg} = atan2( dr{iseg}, dx{iseg} );           % Panel slope angle
        xcp{iseg} = 0.5 * ( xep(2:end) + xep(1:end-1) );     % Panel x center point
        rcp{iseg} = 0.5 * ( rep(2:end) + rep(1:end-1) );     % Panel r center point

        if ( ~props{iseg} ) % 'Normal' components
            npans{iseg} = length( xcp{iseg} );
            npan = npan + npans{iseg};
        else % Actuator disk doesn't contribute panels
            npans{iseg} = 0;
        end
    end


    u = nan( 1, npan );
    v = nan( 1, npan );

    k = 1;
    for iseg=1:nseg
        for isurvey=1:npans{iseg}
            uv = W;
            vv = 0;

            for jseg=1:nseg
                if ( iseg == jseg )
                    mask = true( 1, npans{jseg} );

                    istart = max( 1, isurvey - 3 );
                    iend = min( npans{jseg}, isurvey + 3 );

                    mask(istart:iend) = false;

                    for jpan=istart:iend
                        [ uj, vj ] = ringsubvortex2( xcp{iseg}(isurvey), rcp{iseg}(isurvey), xepts{jseg}(jpan), repts{jseg}(jpan), xepts{jseg}(jpan+1), repts{jseg}(jpan+1) );
                        uv = uv + sum( uj .* gammatgt{jseg}(jpan) .* ds{jseg}(jpan) .* sign(dx{jseg}(jpan)) );
                        vv = vv + sum( vj .* gammatgt{jseg}(jpan) .* ds{jseg}(jpan) .* sign(dx{jseg}(jpan)) );
                    end

                    [ uj, vj ] = ringvortex( xcp{jseg}(mask), rcp{jseg}(mask), xcp{iseg}(isurvey), rcp{iseg}(isurvey) );
                    uv = uv + sum( uj .* gammatgt{jseg}(mask) .* ds{jseg}(mask) .* sign(dx{jseg}(mask)) );
                    vv = vv + sum( vj .* gammatgt{jseg}(mask) .* ds{jseg}(mask) .* sign(dx{jseg}(mask)) );

                else
                    [ uj, vj ] = ringvortex( xcp{jseg}, rcp{jseg}, xcp{iseg}(isurvey), rcp{iseg}(isurvey) );
                    uv = uv + sum( uj .* gammatgt{jseg} .* ds{jseg} .* sign(dx{jseg}) );
                    vv = vv + sum( vj .* gammatgt{jseg} .* ds{jseg} .* sign(dx{jseg}) );
                end
            end

            u(k) = uv;
            v(k) = vv;

            k = k + 1;
        end
    end

    Vmag = sqrt( u.^2 + v.^2 );

    dri = v.*ds{1}./Vmag;

    n = length(dri) + 1 - 2;
    M = eye( n, n - 1 );
    M(n,:) = ones( 1, n - 1 );

    dri(1) = rsol{1}(1)*2;
    dri(end) = -rsol{1}(end)*2;

    b = [dri(2:end-1) -(dri(1)+dri(end))]';

    drnew = (M\b)';

    dri(2:end-1) = drnew;

    rep = abs( [0 cumsum(dri)] );

    relax = 0.5;
    repts{1} = ( 1 - relax ) * repts{1} + relax * rep;


    figure(1)
    plot(xepts{1},repts{1},'k')

end

axis equal
hold off


figure(2)
plot(xepts{1},repts{1},'k')
hold on
for iseg=1:length( xsol )
    plot( xsol{iseg}, rsol{iseg},'k--' )
end

for iseg=1:length( xorig )
    plot( xorig{iseg}, rorig{iseg},'k:' )
end

hold off
axis equal

figure(3)
plot(xcp{1},2*Vmag);
hold on
plot(xcp{1},Vexs{1});
