clc; close all; clear all;

% This ExampleScript uses the FD-BPM to study a simple straight coupler,
% composed of two parallel (straight) waveguides, on SOI technology.

addpath( [pwd, '\BPM'] ); % Required path to access BPM-related m-files
addpath( [pwd, '\Solver'] ); % Required path to access solver m-files
addpath( [pwd, '\Misc'] ); % Some more misc functions

% =========================================================================
% Structure Parameters
% =========================================================================

% Fixed parameters
nsg = [1.45 3.2]; % [] refr.indices of substrate & guiding-layer 
wl  = 1.55; % [um] operating wavelength
wid = 0.25; % [um] x-width of waveguide (input plane)

% Coupler parameters
LDC =  100;  % [um] z-length of Directional Coupler
gap = 0.30;  % [um] x-gap, between coupler waveguides
dn  = 0.02;  % [.] refr.index change in one waveguide core (e.g. 0 or 0.05)

% Finite-Difference (FD) discretization of x-axis 
% *** This affects considerably  the FD-BPM accuracy! Checking your 
%     system's free RAM, it should range between [1e3,1e5] for 6GB systems.
Nx = 1000; % [.] number of samples on x-axis

% PML (Perfectly Match Layers) for reflectionless absorption -- Be careful!
PML_Thk = 1*[ 1 1 ] ; % [um] PML Thickness for up/down layer, in x-dim
PML_Str = 1*[ 1 1 ] ; % [.] PML Absorption "strength"
PML_params = [ PML_Thk , PML_Str ];

% =========================================================================
% Theory/Mode Analysis -- Estimate coupling-length
% =========================================================================

% We will use the MLSWG function (MultiLayer Slab WaveGuide -- Check the 
% function's help for more details on the inputs) to calculate the coupling
% length of the coupler.
%
% The MLSWG might NOT converge in a solution. In this case, we restart the
% simulation until exactly two modes are found. If more that 2 modes are
% found, this means that our system is not single-mode (the lone
% waveguide), so we typically need to reduce guiding layer width, or
% decrease its refr.index.
N_iterations = 0;
while N_iterations < 10;
    neffs = MLSWG( 'TE' , wl , nsg(1)*[1 1] , nsg([2 1 2]) , [wid gap wid] );
    if length(neffs) == 2
        Lc = 0.5*wl / abs(diff(real(neffs))); % [um] coupling length
        disp( ' ** OK: MLSWG found EXACTLY 2 modes :)' );
        break;
    elseif length(neffs) >= 3
        disp( ' ## ERROR: MLSWG solver found MORE than 2 supermodes!' );
        return;
    else
        disp( ' ## WARNING: MLSWG solver MISSED a supermode! ReStarting...' );
        N_iterations = N_iterations + 1;
        if N_iterations == 10
            disp( ' ## ERROR: MLSWG failed to find 2 supermodes.' );
            return;
        end        
    end
end

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

% EXAMPLE: A single-module waveguide coupler
xc1 = -(gap+wid)/2; % x-center of 1st waveguide (same for input/output)
xc2 = +(gap+wid)/2; % x-center of 2nd waveguide (same for input/output)
N_steps = 4*round(LDC/wl); % BPM-step size ~wl/4
XBR1 = [ xc1 xc1 wid wid  0 2 ]; % [ xc, xco, wi, wo, dn, BT ]
XBR2 = [ xc2 xc2 wid wid dn 2 ]; % [ xc, xco, wi, wo, dn, BT ]
SL = { nsg , LDC , N_steps , XBR1 , XBR2 };

% =========================================================================
% Finite-Difference (FD) Discretization / TopView
% =========================================================================

% Discretize structure's cross-section plane (x-axis)
x = linspace( -2-(gap/2+wid) , (gap/2+wid)+2 , Nx ) ; % [um] x-axis vector

% Get index profile and z-axis discretized vector
% -- [ns,z,zVarFlags,zxLines,StepsTotal]=BPMFD2D_PreProcLayout( SL,x,PML_params );
[ ns_zxProf, z, ~, ~, sx_zxProf ]=BPMFD2D_PreProcLayout( SL,x,PML_params );

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
x_offset = (wid+gap/2); % Offset required to focus beam on bottom waveguide
[ neffs , modeProfs ] = MLSWG( 'TE' , wl ,nsg(1)*[1 1] , nsg(2) , wid , x+x_offset );

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

% Measure the power in each side of the coupler
dx = mean(diff(x)); % [um] x-step
P_Left  = sum( abs( FA( x > 0 , : ) ).^2 , 1 )*dx; % bus waveguide
P_Right = sum( abs( FA( x < 0 , : ) ).^2 , 1 )*dx; % tap waveguide
P_Input = P_Left(1)+P_Right(1);

% =========================================================================
% Plots
% =========================================================================

% Power in L/R waveguides along the waveguide direction
figure;
plot( z , P_Left /P_Input , 'b' ); hold on;
plot( z , P_Right/P_Input , 'r' );
legend( 'Left WG (x>0)' , 'Right WG (x>0)' );
xlabel( 'z-axis (um)' ); 
ylabel( 'Normalized Power' );
 
% TopView -- Intensity
figure;
pcolor(z , x, 10*log10( abs(FA).^2 ) ); shading flat;
colorbar; 
caxis( [-30 0] );
title( 'BPM calculated : |Ey| (dB)' );
xlabel( 'z-axis (um)' ); 
ylabel( 'x-axis (um)' );

% Plot SL lines over intensity:
BPMFD2D_DrawLayout( gca, SL, x, PML_params, 'w' );  

% TopView -- Phase
figure;
pcolor(z , x, Phase/pi*180); shading flat;
colorbar; colormap( hsv );
title( 'BPM calculated : angle(Ey) (deg)' );
xlabel( 'z-axis (um)' );  
ylabel( 'x-axis (um)' );

% Plot SL lines over intensity // Black color!
BPMFD2D_DrawLayout( gca, SL, x, PML_params, 'k' );  

% Array the figures in the computer screen
fmfp;