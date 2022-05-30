function [ x0 ] = interpinv(x,y,y0,tangent_Tol)

% "Graphically" find all the x0s where a y(x) curve crosses the y==y0 line.
%
% Given the x,y vectors, where y=function(x), and the value y0, this
% mini-function computes using linear interpolation the value(s) of
% x0 where y=y0, i.e. y(x0)=y0. It segments the y-y0 function in
% monotonically invariant parts, and uses the built-in interp1 
% matlab function in a piece-wise approach. The roots are return
% sorted in ascending order. Tangent_Tol is a limit for tangent-like 
% discontinuities, that can be confused for roots. Set to ~10.
% 
% Alexandros Pitilakis / Thessaloniki, Greece
%  2015 Nov : Original Version

% Test
if nargin ==0 && nargout==0    
   x = linspace(0,1,100);
   y = tan(x*10)-1; %sin(2*pi*x).*exp(3*x);
   y0 = 0;
   tangent_Tol = 10;
end

% Set a limit for the tangent-like discontinuities.
if nargin == 3
    tangent_Tol = 10;
end

% -------------------------------------------------------------------------
ya = y-y0; %Look for zero-crosses in the aux function ya=y-y0;
sch = sign(ya(1:end-1)) .* sign( ya(2:end) ); %-1 where sign-change

% Look for tangeant-like infinties who can be confused for roots.
% So, set this flag to 1, where the function is continuous... 
tanch = abs( ya(1:end-1) - ya(2:end)) < tangent_Tol; 

% Where sch==-1 AND tanch==1, it means that we have a "valid" sign-change
% between ya(is) and ya(is+1)
is = find( (sch == -1) & (tanch == 1) ); 

x0 = NaN*zeros( size(is) );
for k = 1 : length(is)
    x0(k) = interp1( y(is(k):is(k)+1) , x(is(k):is(k)+1) , y0 , 'pchip' );
end
% -------------------------------------------------------------------------

% Test
if nargin ==0 && nargout==0
    clc; close all
    plot(x,y,'bo-'); hold on;
    stem( x0 , y0*ones(size(x0)) , 'ro' )
end