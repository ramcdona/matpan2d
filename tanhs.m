% Hyperbolic tangent stretching function
function x=tanhs(DS0,DS1,xii)

% J.F. Thompson, Z.U.A. Warsi, and C.W. Mastin, NUMERICAL GRID
% GENERATION, North-Holland, New York, pp. 307-308 (1985).
%
% M. Vinokur, On One-Dimensional Stretching Functions for
% Finite Difference Calculations, J. of Comp. Phys., Vol. 50,
% pp. 215-234 (1983).

epsb = 1e-3;

a     = sqrt(DS1/DS0);
b     = 1.0/sqrt(DS0*DS1);
u1    = 0.5;
u2    = 1.0;
u3    = 2.0;
u4    = 0.5;

%
%    Generate hyperbolic tangent arc length distribution for N number
%    of points and report maximum stretching ratio.
%
oneminus  = 1.0 - epsb;
oneplus  = 1.0 + epsb;

x = xii;  % initialize to size and first/last entries.
% x(1)  = xii(1);
% x(end)  = xii(end);

N=length(x);

if (b <= oneminus)
    delta = asinc( b );
    hdelta = 0.5 * delta;
    tnh2 = tan( hdelta );
    for i = 2:N-1
        x(i) = u1 * ( u2 + tan( hdelta * (xii(i) / u1 - u2 ) ) / tnh2 );
    end
elseif (b >= oneplus)
    delta = asinhc(b);
    hdelta = 0.5*delta;
    tnh2 = tanh(hdelta);
    for i = 2:N-1
        x(i) = u1 * ( u2 + tanh( hdelta * ( xii(i) / u1 - u2 ) ) / tnh2 );
    end
else
    ubm = u3 * ( 1.0 - b );
    for i = 2:N-1
        x(i) = xii(i) * ( 1.0 + ubm * ( xii(i) - u4 ) * ( xii(i) - 1.0) );
    end
end

%      Rescale coordinates
am = 1.0 - a;
for i=2:N-1
    x(i) = x(i) / ( a + am * x(i) );
end

return
end

function [delta, approx] = asinhc(b)
% Consider algorithm from
% http://mathforum.org/kb/message.jspa?messageiD=449151

maxiter = 4; % Maximum number of iterations
delmin=5.0E-5;
tol=1.0E-6;

% use series expansions to get initial guess for delta where
% sinh(delta)/delta = b
% Then use Newton iterations to converge delta. Since b is never
% below (1+epsb) in the calling sequence, delta should remain
% well above delmin.

if ( b <= 2.7829681178603 )
    a1= -0.15;
    a2=  0.0573214285714;
    a3= -0.024907294878;
    a4=  0.0077424460899;
    a5= -0.0010794122691;

    bb = b - 1.0;
    delta = sqrt( 6.0 * bb ) * ( ( ( ( ( a5 * bb + a4 ) * bb + a3 ) * bb + a2 ) * bb + a1 ) * bb + 1.0);
else
    c0= -0.0204176930892;
    c1=  0.2490272170591;
    c2=  1.9496443322775;
    c3= -2.629454725241;
    c4=  8.5679591096315;

    v = log( b );
    w = 1.0 / b - 1.0 / 35.0539798452776;
    delta = v + log( 2.0 * v ) * ( 1.0 + 1.0 / v ) + ( ( ( c4 * w + c3 ) * w + c2) * w + c1 ) * w + c0;
end

approx = delta;

%    Newton iterations
for i=1:maxiter
    if ( abs( delta) >= delmin )
        sdelta = sinh( delta );
        cdelta = cosh( delta );
        f = sdelta / delta - b;
        fp = ( delta * cdelta - sdelta ) / ( delta * delta );
        dd = -f / fp;
    else
        disp( 'delta fell below delmin.' );
        delta = 0.0;
        return;
    end

    if ( abs( f ) < tol )
        %disp( strcat('converged delta=', num2str(delta),'   f=',num2str(f),'   ddelta=',num2str(dd),'   i=',num2str(i) ) );
        %semilogy(1:i,errhist,'x-')
        return;
    end

    delta = delta + dd;
end

%disp( 'Exceeded max number of iterations.' );
%semilogy(1:i,errhist,'x-')
%disp( strcat('delta=',num2str(delta),'   f=',num2str(f),'   ddelta=',num2str(dd),'   i=',num2str(i)) );

return;
end

% inverse of sinc on [ 0, 1 ] only.
function [delta, approx] = asinc(b)

maxiter=4;
delmin=5.0E-5;
tol=1.0E-6;

%
%    use series expansions to get initial guess for delta where
%    sin(delta)/delta = b
%    Then use Newton iterations to converge delta. Since b is never
%    above (1-epsb) in the calling sequence, delta should remain
%    above delmin.
%
if ( b <= 0.2693897165164 )
    c3= -2.6449340668482;
    c4= 6.7947319658321;
    c5=-13.2055008110734;
    c6= 11.7260952338351;

    delta = pi * ( ( ( ( ( ( c6 * b + c5 ) * b + c4) * b + c3 ) * b + 1.0 ) * b - 1.0 ) * b + 1.0 );
else
    a1= 0.15;
    a2= 0.0573214285714;
    a3= 0.0489742834696;
    a4= -0.053337753213;
    a5= 0.0758451335824;

    bb = 1.0 - b;
    delta = sqrt( 6.0 * bb ) * ( ( ( ( ( a5 * bb + a4 ) * bb + a3 ) * bb + a2 ) * bb + a1 ) * bb + 1.0 );
end

approx = delta;

%    Newton iterations

for i = 1:maxiter
    if ( abs( delta ) >= delmin )
        sdelta = sin( delta );
        cdelta = cos( delta );
        f  = sdelta / delta - b;
        fp = ( delta * cdelta - sdelta ) / ( delta * delta );
        dd = -f / fp;
    else
        disp( 'delta fell below delmin.' );
        delta = 0.0;
        return
    end

    if ( abs( f ) < tol )
        %disp( strcat('converged delta=', num2str(delta),'   f=',num2str(f),'   ddelta=',num2str(dd),'   i=',num2str(i) ) );
        %semilogy(1:i,errhist,'x-')
        return
    end

    delta = delta + dd;
end

%disp( 'Exceeded max number of iterations.' );
%semilogy(1:i,errhist,'x-')
%disp( strcat('delta=',num2str(delta),'   f=',num2str(f),'   ddelta=',num2str(dd),'   i=',num2str(i)) );

return
end
