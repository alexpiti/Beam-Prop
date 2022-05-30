% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis

clc; close all; clear all;

% This ExampleScript uses the FD-BPM to study a SOI-waveguide 10dB tap,
% including S-bend sections to realistically distance/decouple the I/O.

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
xop = 4;    % [um] x-position of tap's output ports (input waveguide is at x=0);
LSB = 20;   % [um] z-length of S-Bendss

% Tap parameters
Tap = -10;   % [dB] tapped power ratio
Lcp =  20;  % [um] z-length of tap-coupler -- It will be re-calculated
ga  = 0.5;  % [um] x-gap of tap-coupler waveguides

% Finite-Difference (FD) discretization of x-axis
Nx = 1000; % [.] number of x-samples -- Affects FD-BPM accuracy!

% PML (Perfectly Match Layers) for reflectionless absorption -- Be careful!
PML_Thk = 1*[ 1 1 ] ; % [um] PML Thickness for up/down layer, in x-dim
PML_Str = 1*[ 1 1 ] ; % [.] PML Absorption "strength"
PML_params = [ PML_Thk , PML_Str ];

% =========================================================================
% Theory/Mode Analysis -- Estimate coupling-length
% =========================================================================

% We will use the MLSWG function (MultiLayer Slab WaveGuide -- Check the 
% function's help for more details on the inputs) to calculate the coupling
% length where the tapping takes place, via the "supermodes".
neffs = MLSWG( 'TE' , wl , nsg(1)*[1 1] , nsg([2 1 2]) , [win ga win] );
Lc = 0.5*wl / abs(diff(real(neffs))); % [um] coupling length

% Use sync'd waveguide theory to estimate length required for given tap:
%               P_tap = |sin( kappa * Lcp )|^2       ==> 
%                 Lcp = arcsin( sqrt(P_tap) )/kappa,
% where kappa = pi/2/Lc.
kappa = pi/2/Lc; % [rad/m] coupling coefficient
Lcp = asin( sqrt( 10^(Tap/10) ) )/kappa; % [um] coupling length  

% ATTENTION:
disp( ' ** The I/O S-Bends will also contribute some coupling! ** ');

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

% EXAMPLE: Waveguide tap / Asymmetric

% 1st Module: Input S-Bends waveguide
XBR1 = [ 0       0 win wot 0 2 ]; % [ xc, xco, wi, wo, dn, BT ]
XBR2 = [xop ga+win win wot 0 2 ];
SL1 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };

% 2nd Module: Coupling region
XBR2 = [ ga+win ga+win win wot  0 2 ];
SL2 = { nsg , Lcp , 2*round(Lcp/wl) , XBR1 , XBR2 };

% 3rd Module: Bend "tapped" waveguide away (decoupling)
XBR2 = [ ga+win xop win wot 0 2 ];
SL3 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };

% Lets put everything together, in the ROWS of SL
SL = [ SL1 ; SL2 ; SL3 ];

% =========================================================================
% Finite-Difference (FD) Discretization / TopView
% =========================================================================

% Discretize structure's cross-section plane (x-axis)
x = linspace( -2 , xop+2 , Nx ) ; % [um] x-axis vector

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

% Alternatively -- We can focus excitation on upper input port:
% [ neffs , modeProfs ] = MLSWG( 'TE' , wl ,nsg(1)*[1 1] , nsg(2) , win , x+win/2 - xop );

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

% Measure the device's Insertion Losses (IL), from simulation
Ey_inp = FA(:, 1 ); % [V/m] Ey field at input plane
Ey_out = FA(:,end); % [V/m] Ey field at output plane
dx = mean(diff(x)); % [um] x-step
PT_inp = sum( abs( Ey_inp ).^2 )*dx;
PT_out = sum( abs( Ey_out ).^2 )*dx;
fprintf( ' == Insertion Losses : %+4.2f dB \n' , 10*log10( PT_out/PT_inp ) );

% Measure the tapped power ratio from simulation
P1 = sum( abs( Ey_out( x < 2 ) ).^2 )*dx; % bus waveguide
P2 = sum( abs( Ey_out( x > 2 ) ).^2 )*dx; % tap waveguide
fprintf( ' == Measured Tap : %+4.2f dB \n' , 10*log10( P2/(P1+P2) ) );

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

% Input and output profiles (along x-axis)
figure;
plot( x , 10*log10( abs(Ey_inp).^2 ) , 'b' ); hold on;
plot( x , 10*log10( abs(Ey_out).^2 ) , 'r' );
plot( x , -10+0*x , 'k:' );
legend( 'input' , 'output' );
xlabel( 'x-axis (um)' ); 
ylabel( '|Ey|' );
set( gca , 'YLim' , [-30 0] );

% Array the figures in the computer screen
fmfp;