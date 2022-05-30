function [xR, flag] = myNewtonRaphson( functionHandle , x0 , xLim , ...
    Tol_x , Tol_y , MaxIter , ccSlope , Voc )

% Implementation of the Newton-Raphson root-finding algorithm, intended 
% for the numerical "solution" of the characteristic equation of slab (1D)
% multilayer waveguides.
%
% functionHandle is a handle to a function (written by user) designed to 
% return the value of a function y(x0) *AND* its dy/dx value at x0. In our
% case, the function is something like: 
%
%         [XE,dXE]=MLSWG_CharEq(ne0,nL,nR,ns,wl,ts,ModePol);
%
% See the example/test function written for how to write this function.
% x0 is the starting "guess" for the root (obligatory input) and xLim are
% the [lower,upper] x-boundaries where we're allowed to search for a root.
%
% Alexandros Pitilakis / Thessaloniki, Greece
% 2015 Sept: Original version

% Test inputs
if nargin == 0
   functionHandle = @testFunct; % Check testFunct below
   x0 = -2;
   xLim = [ -3 4 ];
end

% Error check
if nargin == 1
    error( ' ## Newton-Raphson: Not starting guess provided!' );
end

% Default parameters
if nargin == 2
    xLim = [-Inf, +Inf];
end
if nargin < 4
    MaxIter = 100;
end
if nargin < 5
    Tol_x = 1e-10;
end
if nargin < 6
    Tol_y = 1e-10;
end
if nargin < 7
    ccSlope = 0.3; % Set to zero to "bypass" this feature of myNR solver
end
if nargin < 8
    Voc = 1;
end

% Begin NR iterations, until k
xR = [];
flag = [];
for ii = 1 : MaxIter
    
    % Find the function- and derivative-value at x=ne0, using the
    % functionHandle input. In our case, these are acquired by a function
    % of the sort [XE,dXE]=MLSWG_CharEq(ne0,nL,nR,ns,wl,ts,ModePol);
    [ y, dy ] = functionHandle( x0 );
    
    % Controlled convergense (stepping) of the method. The default value of
    % cc==1, but here we start it from something smaller and gradually
    % increase it as we approach MaxIter. The paramter ccSlope controls the
    % shape of this increase. If ccSlope>1 then cc increases fast and goes
    % into "saturation" (like a log function). If ccSlope<1, then cc
    % incrases slowly at first and rapidly near MaxIter (like an exp
    % funct). ccSlope>1 typically leads to faster convergence for "regular"
    % functions, but it may miss some roots where the y-slope is steep.
    cc = (ii/MaxIter)^ccSlope; 
    x0 = x0 + cc*( -y./dy );
    
    % Check the algorithm
    if abs(y./dy) < Tol_x && abs(y) < Tol_y && ...
            real(x0) > min(xLim) && real(x0) < max(xLim)
        if Voc == 1
            fprintf( ' ** Newton-Raphson: Root Found x = %f (after %d iters)\n' , x0 , ii );
        end
        xR = x0;
        flag = 1;
        break;
    elseif real(x0) < min(xLim) || real(x0) > max(xLim)
        if Voc == 1
            fprintf( ' ## Newton-Raphson: Off x-Limits => Aborting.\n' );
        end
        flag = 2;
        break;
    elseif ii == MaxIter,
        if Voc == 1
            fprintf( ' ## Newton-Raphson: MaxIterations reachead.\n' );
        end
        flag = 3;
    end
    
end

% Test outputs
if nargin == 0 && all( ~isinf(xLim) )
    close
    x = linspace( min(xLim) , max(xLim) , 1e3 );
    y = functionHandle( x );
    plot( x , y , 'b-' ); hold on;
    plot( xR*[1 1],get(gca,'YLim'),'r-')
    plot( x , 0*x , 'k:' );
end

end

% Test function 
function [f,df]=testFunct(x)
    f = 3*x.^2 + 2*x - 8;
    df = 6*x + 2 ;
end