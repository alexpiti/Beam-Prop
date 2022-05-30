function [ ns , zAxisVector , zVarFlags , zxLines , sx ] = ...
    BPMFD2D_PreProcLayout( Struct_Layout , x , PML_params )

% Calculates various Structure_Layout (SL) parameters to be used in a 
% Finite Difference (FD) Beam Propagation Method (BPM). Check DoProp
% function for more details.
%
% === Inputs ===
%  - Struct_Layout : NM-by-? CELL array that describes the structure
%      (refer to DoProp routine help for more details)
%  - x : 1-by-Nx vector with the FD discretization of the x-axis. This must
%      be in the format returned by MATLAB's linspace() routine
%  - PML_params : 1-by-4 array with the Thickness [up,low] & absorption
%      Strength [up,low] of the Perfectly Matched Layers used to truncate
%      the compuatational window on the two sides of the x-axis.
% 
% === Outputs ===
%  - ns : Nx-by-Nz 2D array, with the (complex) refr.indices of the SL
%  - zAxisVector : 1-by-Nz FD discretization of the z-axis (longitudinal)
%  - zVarFlags : 1-by-zMod [boolean] determines if a zModule is "straight"
%  - zxLines : Set of curves in zx-plane defining the sidewalls of WGs
%  - sx : Nx-by-Nz 2D array, with the PML tensors
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2010 Sep : Original version
%  2015 Nov : Revised version

% Test inputs
if nargin == 0
    
    % Sample Y-combiner
    XBR1 = [ -2 0 2 2   0  2 ];
    XBR2 = [  2 0 2 2   0  2 ];
    XBRX = [  0 0 2 1 -0.1 2 ];
    Struct_Layout(1,:) = { [1 2] , 10 , 20 , XBR1 , XBR2 };
    Struct_Layout(2,:) = { [1 2] ,  5 , 11 , XBRX , []  };
    x = linspace( -5 , +5 , 101 );    
    
    PML_params = [ 1*[1 1] , 1*[1 1] ]; % Thickess [up,low], Str [up,low]
    
end

%Even for Imaginary Distance BPM, these should be turned real!
for kkm = 1 : size( Struct_Layout , 1 )
    Struct_Layout{kkm,2} = abs( Struct_Layout{kkm,2} );
end

%For easier access
SL = Struct_Layout;

% ========================================================================
% "Easy-to-Get" outputs, directly from SL
% ========================================================================

%Total Steps
StepsTotal = sum( cell2mat(SL(:,3)) );

%Total z-Distance
LzTotal = sum( cell2mat(SL(:,2)) );

%Number of z-Modules
NzMs = size( SL , 1 );

%Grid-Size in the lateral (cross-sectional) dimension
Nx = length( x );

% ========================================================================
% Scan z-Modules, for the module-dependent variables
% ========================================================================

%Init module-dependent output vects
zAxisVector = 0;
zVarFlags = NaN * ones(1,NzMs);
ns = NaN * zeros( Nx , StepsTotal+1 ); 
zxLines = cell( NzMs , size(SL,2) - 3 ); % [N-Modules by (max)-N-Branches 

% Scan Modules
for kkm = 1 : NzMs

    %Primary Information for zx-Layout from SL array
    NBs = size(SL,2) - 3; %(max) Number of x-branches (1=single,2=coupler etc)
    for kkb = NBs : -1 : 1 %Reduce NBs for this module, appropriately
        if isempty( SL{kkm,kkb+3} ), NBs = NBs-1 ; end
    end
    Lz  = SL{kkm,2}; %z-Length in [um]
    Nzs = SL{kkm,3}; %Number of z-steps
    
    %Get the info of the total zAxisVector
    zStart = sum( cell2mat( Struct_Layout(1:kkm-1,2) ) );
    zAxisVector = [ zAxisVector , zStart + (1:Nzs)/Nzs*Lz ];
    
    %z-Varying check:
    for kkb = 1 : NBs %Scan x-Branches
        %Check if collumns 1=2 (x-centers) and 3=4 (x-widthts)
        Bi = SL{kkm,3+kkb};
        if abs(Bi(1)-Bi(2))>1e-10 || abs(Bi(3)-Bi(4))>1e-10
            zVarFlags( kkm ) = 1; break
        end
        zVarFlags( kkm ) = 0;
    end; clear Bi;
    
    %Init RefrIndx array for this module:
    n_sub = SL{kkm,1}(1);
    n_gui = SL{kkm,1}(2);
    nMs = ones( Nx , Nzs ) * n_sub;
        
    % ========================================================================
    % Scan x-Branches to fill the RefrIndx array
    % ========================================================================
    for kkb = 1 : NBs
        
        %Get zx-Layout info for this x-branch
        Bi = SL{kkm,3+kkb};

        %Get I/O x-centers and x-widths for this {Branch,Layer}
        xi = Bi(1); xo = Bi(2); wi = Bi(3); wo = Bi(4);
        
        %Apply possible Index-Shift in this Module/Branch
        n_gui_this_xBr = n_gui + Bi(5); 

        %Sample kkk (intra-module step-counter)
        kkk = 1 : Nzs;
        
        %Calculate Tapering (Linear):
        wkkk = wi + (wo-wi)*kkk/Nzs;

        %Calculate x-centers and x-widths
        Dx = (xo-xi);
        switch Bi(6) % Bend-Type
            case 1 %Straight
                xc = xi + Dx*kkk/Nzs;
                xw = wkkk / cos(atan(Dx/Lz));
            case 2 %S-Bend
                zn = ( kkk/Nzs*Lz - Lz/2 ) * pi / Lz;
                xc = xi + Dx/2 + Dx/2*sin(zn);
                xw = wkkk ./ cos(atan(Dx/Lz*pi/2*cos(zn)));
            otherwise
                error( ' ## BPMFD2D_PreProLayout: Incorrect BendType!' );
        end
    
        %Calculate Up/Down Wall-Arrays for this z-Module
        DoWA = ones( Nx , 1 ) * ( xc(:)' - xw(:)'/2 );
        UpWA = ones( Nx , 1 ) * ( xc(:)' + xw(:)'/2 );
        
        %Set indices in the RefrIndx array of this module
        xxs = x(:) * ones( 1 , Nzs );
        iis = UpWA > xxs & DoWA < xxs ;
        nMs( iis ) = n_gui_this_xBr ;
        
        %Fill the zxLines cell-array:
        dz = Lz/Nzs; %step-size
        auxi = [ zStart+(dz:dz:Lz) ; xc-xw/2 ; xc+xw/2 ];
        auxi = [ [zStart;auxi(2:3,1)] , auxi ];
        zxLines(kkm,kkb) = { auxi }; clear auxi
        
    end
    
    %Add this module's RefrIndx array in the total one:
    NstepsUpToPrev = sum( cell2mat( SL(1:(kkm-1),3) ) );
    ns( : , (1:Nzs)+1+NstepsUpToPrev ) = nMs;
    
    clear nMs iis Lz Nzs NBs 

end

% The indices of the global-input x-Slice
ns( : , 1 ) = ns( : , 2 );

% =========================================================================
% PML Termination
% =========================================================================

% PMLs take-up the start/end parts of the x-vector defined above (not extra!)

PML_Thk = PML_params(1:2); % [Upper Lower] PML thickess // low: x<0, up: x>0
PML_Str = PML_params(3:4); % [Upper Lower] PML strength

sx = ones(size(ns)); %Initialize PML tensor sx=1-j*tan(delta)

UpPML = x > (x(end)-PML_Thk(1)) ; %Logical "1s" for upper-PML
LoPML = x < (x(1)  +PML_Thk(2)) ; %Logical "1s" for lower-PML

sx( UpPML , : ) = ( 1 - 1j*PML_Str(1)*linspace( 0,1,sum(UpPML) ).^2 ).' ...
    * ones(1,StepsTotal+1); %Quadratic Profile
sx( LoPML , : ) = ( 1 - 1j*PML_Str(2)*linspace( 1,0,sum(LoPML) ).^2 ).' ...
    * ones(1,StepsTotal+1); %Quadratic Profile


% =========================================================================
% Test Ouput Vars
% =========================================================================
if nargin == 0
   close all; clc;
   fprintf( ' Number of z-Modules : %d\n' , NzMs );
   fprintf( ' Total z-Steps       : %d\n' , StepsTotal );
   fprintf( ' Total z-Distance    : %d\n' , LzTotal );
   disp(zVarFlags)
   disp(zAxisVector)
   imagesc( zAxisVector , x , ns ); colorbar
end
