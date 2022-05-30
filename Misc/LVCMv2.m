function [ J ] = LVCMv2( CmapSelect , M , CNNDs )

% Produces various fancy colormaps. Check the comments for more info.
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2009 circa : Original version

%Input-Check/Set
if nargin == 2
    CNNDs = []; %Default Value (Equally Spaced Colors)
elseif nargin == 1
    CNNDs = [];
    M = 128 + 1;
elseif nargin == 0
    CNNDs = []; %Color-Node Normalized Distances
    M = 128 + 1; 
    CmapSelect = 7;
end

%If CmapSelect is a Nc-by-3 vector, then use the colors in each row of this
%matrix to create a Linearly_Varying Colormap, with these colors. Else if
%CmapSelect is a scalar (1..11 etc), this creates a preset colormap.
if ~isscalar( CmapSelect )
    cs = CmapSelect ;
else

    switch CmapSelect
        
        %Good for Abs Values ( 0 ... 1 )
        case 2 , cs = [ 0 0 0 ; 1 0 0 ; 1 1 0 ; 1 1 1 ]; %Hot
        case 3 , cs = [ 1 1 1 ; 1 1 0 ; 1 0.5 0 ; 1 0 0 ]; %White-Heat
        case 4 , cs = [ 0 0 0 ; 0 1 0 ; 1 1 0 ; 1 1 1 ]; %Mojito
        case 5 , cs = [ 0 0 0 ; 0 0 1 ; 0 1 1 ; 1 1 1 ]; %SapFire
        case 11, cs = [ 1 1 1 ; 0 1 1 ; 0 0 1 ; 0 0 0]; %inv SapFire
        case 12, cs = [ 0 0 0 ; .5 0 1 ; .8 0 0 ; 1 1 0 ]; % EyePat:Black/Purple/Dark-Red/Yellow
        case 12.1,cs= [ 0 0 0 ; .5 0 .5 ; 1 0 0 ; 1 1 0 ; 1 1 1]; %EyePat Mod v1
        case 12.2,cs= [ 0 0 0 ; 0 0 1 ; 1 0 0 ; 1 1 0 ; 1 1 1 ]; CNNDs = [ 0.5 1 1 0.5 ];%EyePat Mod v2
        case 12.3,cs= [ 1 1 1 ; .5 0 .5 ; 1 0 0 ; 1 1 0 ]; %EyePat Mod v3 (white)
        case 13, %Black-Jet
            cs = [ 0 0 0 ; 0 0 0.5 ; 0 0 1 ; 0 1 1 ; 1 1 0 ; 1 0 0 ; 0.5 0 0 ]; 
            if isempty( CNNDs )
                CNNDs = [1 2 2 2 2 2]; 
            end
        case 13.1 %White-Jet
            cs = [ 1 1 1 ; 0 0 1 ; 0 1 1 ; 1 1 0 ; 1 0 0 ; 0.5 0 0 ]; 
            if isempty( CNNDs )
                CNNDs = [1 2 2 2 2]; 
            end
        case 14, cs = [ 0 0 0 ; 0.5 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 1 ; 0 0 1 ; 0 0 0.5 ]; %Inv-Black-Jet
        case 15, cs = [ 0 0 0 ; 0.5 0 0 ; 1 0 0 ; 1 1 0 ; 0 1 1 ; 1 1 1 ]; %Mod-Inv-Blk-Jet
        case 16, cs = [ 1 1 1 ; 1 1 0 ; 1 0 0 ; 0 0 1 ; 0 1 1 ]; %Inverted-Jet
        case 17, cs = [0 .7 .9; 1 1 0 ; 1 0 0; .5 0 0]; % Cyan -> Yellow -> Red
    
        %Good for Real Values ( -1 ... 0 ... 1 )
        case 1 , cs = [ 0 0 0.5 ; 0 0 1 ; 0 1 1 ; 1 1 0 ; 1 0 0 ; 0.5 0 0 ]; %Jet
        case 6 , cs = [ 0 0 0.5 ; 0 0 1 ; 0 1 1 ; 1 1 1 ; 1 1 0 ; 1 0 0 ; 0.5 0 0 ]; %AltJet
        case 7 , cs = [ 0 0 0.5 ; 0 0 1 ; 1 1 1 ; 1 0 0 ; 0.5 0 0 ]; %Deep-LaFrance
        case 7.1, cs = [ 0 0 0.5 ; 0 0 1 ; .9 .9 .8 ; 1 0 0 ; 0.5 0 0 ]; %Deep-LaFrance (Cream-White)
        case 8 , cs = [ 0 0 1 ; 1 1 1 ; 1 0 0 ]; %LaFrance
        case 9 , cs = [ 0 1 1 ; 0 0 1 ; 0 0 0 ; 1 0 0 ; 1 1 0 ]; %Black AltJet
        case 10, cs = [ 0 1 1 ; 0 0 1 ; 1 1 1 ; 1 0 0 ; 1 1 0 ]; %White AltJet
    end

end
Nc = size(cs,1); %Number of Color-Nodes Requested

%Default CNNDs: Equally Spaced Colors
if isempty( CNNDs )
    CNNDs = ones(1,Nc-1) / (Nc-1); 
end

%Checks!
if length( CNNDs ) ~= (Nc-1)
    error( ' ** Invalid Number of CNNDs! (Number-Color-Nodes minus one!)...' )
end

%Normalize (0,1) the CNNDs
if sum( CNNDs ) ~= 1, 
    CNNDs = CNNDs / sum( CNNDs );
end

%Generate Lin-Var-Color-Map: The values of [u1] are ranged in [1,Nc]. For example if u1(45)=2.78, 
%this means that the color with index 45 will be 78% that of cs(3,:) mixed with 22% that of cs(2,:).
Mcs  = [ 0 cumsum(CNNDs*(M-1)) ]; %At which "continuous" M-index is each color of [cs] matrix
u1 = interp1( Mcs+1 , 1:Nc , 1:M )'; %Between which Color-index each M-index falls
u1 = max( min( u1 , Nc ) , 1 ); % min/max it, in [1,Nc], so that it doesn't fall off-range
u2 = [ floor( u1 ) , ceil( u1 ) , 1-mod( u1 , 1 ) , mod( u1 , 1 ) ]; % floor, ceil, 1-mod, mod
J = cs( u2(:,1) , : ) .* (u2(:,3)*ones(1,3)) + ...
    cs( u2(:,2) , : ) .* (u2(:,4)*ones(1,3)) ;

%Test Plot
if nargin == 0 && nargout == 0
    close all; %clc
    u = peaks; u = u / max(abs(u(:)));
    
    subplot(2,1,1)
    surf(u); view(2); colorbar; shading flat; caxis( [-1 +1] );
    colormap(J);
    axis off; title(' real/imag , min/max')
    
    subplot(2,1,2)
    surf(abs(u)); view(3); colorbar; shading interp; caxis( [0 +1] );
    colormap(J);
    axis off; title( 'abs , zero/max' )
    
end


