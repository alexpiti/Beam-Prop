clear all; close all; clc;

% This ExampleScript studies the supermodes of a SOI-slab WG coupler using
% a characteristic-equation solver (MLSWG)

addpath( [pwd, '\Solver'] ); % Required path to access solver m-files
addpath( [pwd, '\Misc'] ); % Some more misc functions

% =========================================================================
% Parameters
% =========================================================================
wl  = 1.55; % [um] wavelength
nLR = [ 1.45 1.45 ]; % [.] refr.index of L/R semi-inf slabs (substrate/cladding)
ns  = [ 3.20 1.45 3.20 ]; % [.] refr.index of intermediate layers
wid = 0.30; % [um] x-width of "core" layers
gap = 0.50; % [um] x-gap (substrate) separating the cores
ts  = [ wid gap wid ]; % [um] thickness of intermediate layers 

% =========================================================================
% MLSWG - Characteristic equation solver // NewtonRaphson
% =========================================================================

% Define an x-space for displaying the mode profiles
xPlot = linspace( -1 , sum(ts)+1 , 400 ); % [um] cross-section x-axis

% Solve char-eq using Newton-Raphson method and calculate modes on xPlot
[ neXE , modeProfsXE ] = MLSWG( 'TE' , wl , nLR , ns , ts , xPlot );

% Plot all the mode Profiles 
figure;
linCols = hsv( length(neXE) );
myLegends = cell( 1,length(neXE) );
for ii = 1 : length(neXE)
    plot( xPlot , real(modeProfsXE(ii,:)) , 'Color' , linCols(ii,:) ); hold on;
    myLegends{ii} = sprintf( ' TE%d | neff = %6.4f' , ii , real(neXE(ii)) );
end
legend( myLegends );
set( gca, 'XLim' , xPlot([1 end]) ); % set axis limits
title( 'MLSWG (Newton-Raphson) modes');

