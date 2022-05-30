function BPMFD2D_DrawLayout( AxHandle, SL, x, PMLs, LinCol, xLims )

% Plots the zx-curves the represent the side-walls of all waveuides (also
% called "x-Branches") that make the TopView of the Structure_Layout (SL) 
% to be used in BPM-FD-2D. LinCol (1-by-3) is the color to be used.
%
%  *** Refer to PreProcLayout for more details on Inputs
%  *** xLims input is optional.
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2010 Sep : Original version
%  2015 Nov : Revised version

% Test inputs
if nargin == 0
   AxHandle = gca; cla;
   SL(1,1:5) = { [1 2] , 10 , 25 , [-2 +2 1 3 0 2 ] , [-4 -4 1 2 0 2 ] };
   SL(2,1:5) = { [1 2] , 10 , 25 , [+2  0 3 2 0 2 ] , [-4 -4 2 1 0 2 ] };
   x = linspace( -8,8, 100);
   PMLs = [ 1 1 1 1 ];
   LinCol = [1 0 0];
   xLims = x([1 end]);   
end

% Set xLimits from discretized x-vector
if nargin == 5,
    xLims = x([1 end]);
end

%Set Axes for plot
axes( AxHandle );
hold on;

% Get zxLines from PreProcLayout routine.
[~,~,~,zxLines]=BPMFD2D_PreProcLayout( SL,x,PMLs );

%Number of Modules and max x-Branches
NzM = size( zxLines , 1 );
NxB = size( zxLines , 2 );

%Scan z-Modules:
for kkm = 1 : NzM
    
    %Scan all x-Branches of this z-Module. Each x-Branch has two curves
    for kkb = 1 : NxB
        
        TheLines = zxLines{kkm,kkb};
        
        if isempty( TheLines ), break; end
                
        plot( TheLines(1,:) , TheLines(2,:) , 'Color' , LinCol ); % one wall
        plot( TheLines(1,:) , TheLines(3,:) , 'Color' , LinCol ); % other wall
        
    end
    
    % Plot zModule-separating lines (vertical lines)
    TheLines = zxLines{kkm,1};
    plot( TheLines(1,end)*[1,1] , xLims , 'Color' , LinCol , 'LineStyle' , ':' ); %z=const separation line
    
end

% Plot PML-edges (horizontal lines)
plot( get(gca,'XLim') , (min(x)+PMLs(1))*[1 1] , 'g:' );
plot( get(gca,'XLim') , (max(x)-PMLs(2))*[1 1] , 'g:' );
