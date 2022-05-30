function TheLine = flwcs( chars , L )

% Mnemonic : Fill Line With CharacterS
%
% This functions returns a string of length (L) with a periodic repetion of
% the chars supplied. Default inputs: chars= '|-<>-|', L=CommandWindow
% This function fills the length of the current command window with a
% periodic repetition of the string characters given in chars.
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2009 December

if nargin == 0
    clc; chars = '|-<>-| ';
end
if nargin <= 1
    L  = get(0,'CommandWindowSize')-2;
end


dL = length(chars);
N  = ceil( L(1)/dL ); 
TheLine = repmat( chars , 1 , N );
TheLine = TheLine(1:L-1); %Crop a little from the end, cuz sometimes it leaks

if nargout == 0 && nargin == 0
    disp( TheLine(1:L-1) );
end



