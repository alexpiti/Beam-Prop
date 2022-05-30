function fmfp( Par )

% This mini-function fixes the positions of multiple figures
% for a MATLAB/Monitor configuration like the one used at the
% AUTH-Photonics Lab during my PhD
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2010 October
%  2011 June
%  2011 Sept
%  2015 Nov

% % Tests...
% close all
% for k = 1:12, figure; end
% close(3);

% Position-Fix param:{[0],1}={Cascade,Array}
if nargin == 0, Par = 1; end

%Get all figure handles
Fs = get( 0 , 'Children' );
Fs = flipud( Fs );
Fs = sort(Fs,'ascend');

%Cascade them...
if Par == 0    
    
    %dd is the "cascading-distance" of consecutive Figures.
    dd = round( 480 / (length(Fs)-1) );
    
    for k = 2 : length(Fs)
        auxi = get( Fs(k-1) , 'Position' );
        set( Fs(k) , 'Position' , [auxi(1)+dd,auxi(2)-dd,auxi(3:4)] );
    end
    
end

%Array them...
if Par == 1
    
    Nf = length(Fs);
    Nx = ceil( sqrt(Nf*1.3) ); % Number of figs along x-dim (horiz) // 1.3 is a standard A4 aspect ratio
    Ny = ceil( Nf/Nx ); % Number of figs along y-dim (vert)
    
    ScrDims = get(0,'ScreenSize'); % [ 1 1 W_pix H_pix]
    
    gap  = 10; % pixel gap between figures (i.e. between OS windows)
    bPL  = 50; % Y // bottom PixelsLost (from Windows taskbar) >= 31
    mBPL = 80; % Y // menubar PixelsLost (for each matlab figure) >=69
    
    dx = round( ((ScrDims(3)- 5 )-Nx*(gap + 5 )) / Nx );
    dy = round( ((ScrDims(4)-bPL)-Ny*(mBPL+gap)) / Ny );
    
    kk = 0;
    for ky = 1 : Ny
        for kx = 1 : Nx
            kk = kk + 1;
            if kk > Nf, break; end
            set( Fs(kk) , 'Position' , [5+(kx-1)*(dx+10), bPL+(Ny-ky)*(dy+(mBPL+gap)), dx-gap, dy] );
        end
    end
        
end




