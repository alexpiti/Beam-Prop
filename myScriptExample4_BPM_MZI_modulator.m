% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis

clc; close all; clear all;

% This ExampleScript uses the FD-BPM to study a SOI-waveguide MZI
% amplitude modulator, consisting of a Y-splitter/combiner and phase-tuning
% (arbitrary refr. index perturbation) of one arm of the MZI.

addpath( [pwd, '\BPM'] ); % Required path to access BPM-related m-files
addpath( [pwd, '\Solver'] ); % Required path to access solver m-files
addpath( [pwd, '\Misc'] ); % Some more misc functions

% =========================================================================
% Structure Parameters
% =========================================================================

% Fixed parameters
nsg = [1.45 3.2]; % [] refr.indices of substrate & guiding-layer 
wl  = 1.55; % [um] operating wavelength
win = 0.25; % [um] x-width of waveguide (input plane)
wot = win;  % [um] x-width of waveguide (output plane)
DSB = 4;    % [um] distance of MZI-arms in x-plane
LSB = 20;   % [um] z-length of S-Bendss
LMZI= 50;   % [um] z-length of MZI arms

% Finite-Difference (FD) discretization of x-axis
Nx = 1000; % [.] number of x-samples -- Affects FD-BPM accuracy!

% Modulator controlling parameter
dn = -0.05; % [.] refr.index change in one MZI arm -- See below!
OnOff = 1; % [.] 0=Off, 1=On, other=use_previous dn -- See below!

% PML (Perfectly Match Layers) for reflectionless absorption -- Be careful!
PML_Thk = 1*[ 1 1 ] ; % [um] PML Thickness for up/down layer, in x-dim
PML_Str = 1*[ 1 1 ] ; % [.] PML Absorption "strength"
PML_params = [ PML_Thk , PML_Str ]; 

% =========================================================================
% Structure Layout (SL) -- TopView
% =========================================================================

% The entire structure is made up of "modules", i.e. blocks in the xz-plane
% that are connected in cascade, one after the other. The output field
% profile from one module serves as an excitation/input field for the next
% module. These modules are the *ROWS* in the *CELL* array SL, which is
% defined like this:
%
%    SL(1,:) = { nsg1 , L1 , N_steps1 , XBR1a , XBR1b , []    };
%    SL(2,:) = { nsg2 , L2 , N_steps2 , XBR2a , []    , []    };
%    SL(3,:) = { nsg3 , L3 , N_steps3 , XBR3a , XBR3b , XBR3b };
%    SL(4,:) = ... and so on ...
%
% where "nsg" is a 2x1 matrix with the refractive indices of substrate and
% guiding (core) layer materials, "L" is the z-Length of this module, 
% "N_steps" is the number of BPM steps it is broken into (typically we need
% steps in the order of dz~1/4*wavelegth or less for longitudinally-varying
% modules, and dz~1/2*wavelength for straight modules) and "XBRs" are the 
% "x-branches" that are the various waveguides that compose this module.
%
% These "x-branches" are pairs-of-curves the xz-plane that reresent the 
% side-walls of each waveguide. They are defined by a 6x1 array: 
%
%    XBR = [ xc, xco, wi, wo, dn, BT ] ;
%
% where "xci" (xco) is the x-center of the waveguide in the input (output)
% planes, "wi" (wo) is the x-width of the waveguide in the input (output) 
% plane, "dn" is the refr.index difference of the waveguide core wrt to 
% default value [that is nsg(2)] and "BT" stands for "Bend Type" and its 
% only valid values are {1,2}={Linear,Sinusoidal}, where 2 is the default.

% EXAMPLE: MZI On-Off modulator (single-mode waveguide)

% 1st Module: Input Y-Splitter
XBR1 = [ 0 -DSB/2 win wot 0 2 ]; % [ xc, xco, wi, wo, dn, BT ] 
XBR2 = [ 0 +DSB/2 win wot 0 2 ];
SL1 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };

% 2nd Module:  MZI Arms
XBR1 = [ -DSB/2 -DSB/2 win wot dn 2 ]; % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MZI index-change here!
XBR2 = [ +DSB/2 +DSB/2 win wot  0 2 ];
SL2 = { nsg , LMZI , 2*round(LMZI/wl) , XBR1 , XBR2 };

% 3rd Module: Output Y-Combiner
XBR1 = [ -DSB/2 0 win wot 0 2 ];
XBR2 = [ +DSB/2 0 win wot 0 2 ];
SL3 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };

% 4th Module: Some more straight waveguide
Laux = LSB;
SL4 = { nsg , Laux , 2*round(Laux/wl) , [0 0 win win 0 2 ] , [] };

% Lets put everything together, in the ROWS of SL
SL = [ SL1 ; SL2 ; SL3 ; SL4 ];

% =========================================================================
% Theory/Mode Analysis -- Accurately estimate the "dn" required for MZI
% =========================================================================

% We will use the MLSWG function (MultiLayer Slab WaveGuide -- Check the 
% function's help for more details on the inputs) to calculate the 
% delta-beta (db) required to get pi=db*LZMI. 
dns = linspace( 0 , -0.03 , 6 ); % we will only try a few dn samples
neffs = NaN*ones(size(dns)); % initialize matrix
for ii = 1:length(dns) % repeated calls to MLSWG
    neffs(ii) = MLSWG( 'TE' , wl , nsg(1)*[1 1] , nsg(2)+dns(ii) , win );
end
dPhi = 2*pi/wl * abs( neffs(1)-neffs ) * LMZI; % calc phase

% Let's get a visual inspection of how that looks
figure;
plot( dns, dPhi*180/pi , 'bo-' ); hold on; % plot phase
plot( dns , 180*ones(size(dns)) , 'r' );
set( gca, 'YTick' , 0:45:360 );
xlabel( '\Delta n_{gui}' );
ylabel( '\Delta\Phi = k_0 * \Delta n_{eff} * L_{MZI} (deg)' );

% Use linear-interpolation to get a better estimate
dn = interp1( dPhi , dns , pi );

% Correct the SL matrix with the accurate dn value
if OnOff == 0 % Off ==> dn nonzero
    SL{ 2 , 4 }(5) = dn;
elseif OnOff == 1 % On ==> dn==0
    SL{ 2 , 4 }(5) = 0;
else
    disp( ' ** Using PreSet dn' );
end

% =========================================================================
% Finite-Difference (FD) Discretization / TopView
% =========================================================================

% Discretize structure's cross-section plane (x-axis)
x = linspace( -4 , +4 , Nx ) ; % [um] x-axis vector

% Get index profile and z-axis discretized vector
% -- [ns,z,zVarFlags,zxLines,StepsTotal]=BPMFD2D_PreProcLayout( SL,x,PML_params );
[ ns_zxProf, z ]=BPMFD2D_PreProcLayout( SL,x,PML_params );

% Draw the FD-discretized SL
figure;
pcolor( z , x , ns_zxProf ); shading flat; colorbar;

% Draw (on same axes) the side-walls of all waveguides in the SL
LinCol = 'w' ; % Color of lines
BPMFD2D_DrawLayout( gca, SL, x, PML_params, LinCol ); %axis equal tight;

% Set axes cosmetics
title( 'Structure TopView (refractive index)' );
xlabel( 'z-axis (um)' ); 
ylabel( 'x-axis (um)' );

% =========================================================================
% Generate BPM excitation (on input plane)
% =========================================================================

% We will use the MLSWG (MultiLayer Slab WaveGuide) function to get the
% fundamental TE mode supported by the input waveguide. Check the function
% help for more details on the inputs. We stress that the last input (that
% is the x-vector over which the mode profile will be calculated) needs to
% be offset by +width/2, because in the MLSWG x=0 is NOT in the middle of
% the guiding-layer, but on its left (low-x) side. The effective index of
% the mode will be used as the REFERENCE index of the BPM.
[ neffs , modeProfs ] = MLSWG( 'TE' , wl ,nsg(1)*[1 1] , nsg(2) , win , x+win/2 );

% In case the structure supported more modes:
nref = neffs(1);
Excitation = modeProfs(1,:);

% =========================================================================
% Beam Propagation Method (BPM)
% =========================================================================

% Call a single-routine to do the propagation, from the input of first
% module of SL to the output of its last module. This returns the 2D matrix
% "FA" that contains the field values in ALL the nodes of the FD-mesh.
FA = BPMFD2D_DoProp( SL, x , PML_params , Excitation , nref , wl  );

% =========================================================================
% Post-processing
% =========================================================================

% Post-process
FA = FA / max(abs( FA(:) )); % normalize abs-value to 1
Phase = angle( FA ); % calc phase
Phase( abs(FA).^2<0.001 ) = NaN; % don't calc/plot phase where |Ey|^2<-30dB

% Calculate the overlap between the input and output fields to determine
% the extinction of the modulator (abs & angle)
Fi = FA(:, 1 ); % input excitation
Fo = FA(:,end); % output port
OI = sum( Fi .* conj(Fo) ) / sum( Fi .* conj(Fi) ); % sum(..)==Integral along x-axis (in our case)

% =========================================================================
% Plots
% =========================================================================

% TopView -- Intensity
figure;
pcolor(z , x, 10*log10( abs(FA).^2 ) ); shading flat;
BPMFD2D_DrawLayout( gca, SL, x, PML_params, 'w' );
colorbar;
caxis( [-30 0] );
title( 'BPM calculated : |Ey| (dB)' );
xlabel( 'z-axis (um)' ); 
ylabel( 'x-axis (um)' );

% TopView -- Phase
figure;
pcolor(z , x, Phase/pi*180); shading flat;
BPMFD2D_DrawLayout( gca, SL, x, PML_params, 'k' );
colorbar; colormap( hsv );
title( 'BPM calculated : angle(Ey) (deg)' );
xlabel( 'z-axis (um)' );  
ylabel( 'x-axis (um)' );

% Input and output profiles (along x-axis)
figure;
plot( x , abs(FA(:, 1 )).^2 , 'b' ); hold on;
plot( x , abs(FA(:,end)).^2 , 'r' );
legend( 'input' , 'output' );
xlabel( 'x-axis (um)' ); 
ylabel( '|Ey|' );

% Array the figures in the computer screen
fmfp;