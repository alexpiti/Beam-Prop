function Fi = BPMFD2D_DoProp( Struct_Layout , x , PML_params , xInpProf , nref , wl )

% Implementation of Finite-Difference Beam Propagation Method (FD-BPM).
%
% === Output ===
%  - Fi: The complex field intensity at each FD-node in the zx-plane 
%
% === Inputs ===
%  - SL: The structure is described by "Struct_Layout" (SL) array, which 
%    contains information of the topview, in the xz-plane where x is
%    cross-sectional transverse axis and z is the longitudinal axis. 
%    For more details on how SL is defined, see below.
%  - x : 1-by-Nx vector with the FD discretization of the x-axis. This must
%    be in the format returned by MATLAB's linspace() routine
%  - PML_params: 1-by-4 array with the Thickness [up,low] & absorption
%    Strength [up,low] of the Perfectly Matched Layers used for the 
%    reflectionless absorption of the radiation (truncates the 
%    compuatational window on the two sides of the x-axis).
%  - xInpProf: The BPM excitation at the input plane 
%  - nref: Reference index of the BPM. This should be approximately equal 
%    to the *effective* mode index of the guided mode in-question.
%  - wl: Operating wavelength.
%
% -------------------------------------------------------------------------
% Struct_Layout Description
% -------------------------------------------------------------------------
%
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
%
% *** EXAMPLE: The SL of a MZI on-off modulator
%
%         % 1st Module: Input Y-Splitter
%         XBR1 = [ 0 -DSB/2 win wot 0 2 ];
%         XBR2 = [ 0 +DSB/2 win wot 0 2 ];
%         SL1 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };
% 
%         % 2nd Module:  MZI Arms
%         XBR1 = [ -DSB/2 -DSB/2 win wot dn 2 ]; % Index-change here!
%         XBR2 = [ +DSB/2 +DSB/2 win wot  0 2 ];
%         SL2 = { nsg , LMZI , 2*round(LMZI/wl) , XBR1 , XBR2 };
% 
%         % 3rd Module: Output Y-Combiner
%         XBR1 = [ -DSB/2 0 win wot 0 2 ];
%         XBR2 = [ +DSB/2 0 win wot 0 2 ];
%         SL3 = { nsg , LSB , 8*round(LSB/wl) , XBR1 , XBR2 };
% 
%         % Lets put everything together, in the ROWS of SL
%         SL = [ SL1 ; SL2 ; SL3 ];
%
% *** ATTENTION: This m.FILE contains two internal functions (bottom).
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2015 Nov : Original Version

% =========================================================================
% Default parameters -- Be careful with these!
% =========================================================================

% This displays a message every X steps
MonitorStep = Inf;

% Formulation
BPM = 0    ; % {0,1}={Scalar,Semi-Vector} BPM Formulation to use
MSP = 1    ; % [.] Multi-Step Param (1..4), default (1) is "WideAngle", 0 would be "Paraxial"
aCN = 0.501; % [.] Crank-Nicolson Stability param (must be >=0.5)

% =========================================================================
% Pre-Process Layout
% =========================================================================

% x-axis params
Nx = length(x); % discretization size of cross-section
dx = mean(diff(x)); % Finite-difference descretization size

% Generates various several variables of zBPM
[ns, zAx, zVarFlags, ~, sx] = BPMFD2D_PreProcLayout( Struct_Layout,x,PML_params );

% Number of z-Modules
NzMs = size( Struct_Layout , 1 );

% Number of z-BPM-steps
StepsTotal = length(zAx)-1;

% =========================================================================
% BPM Algorithm
% =========================================================================

fprintf( ' -- Propagating... ');

% Initialize field matrix
Fi = NaN*zeros( Nx , StepsTotal );
Fi( : , 1 ) = xInpProf;

% Propage through all the zModules of the layout of this structure
kkg = 0; % Global step-counter
for kkm = 1 : NzMs
    
    % z-Step for this zModule
    dz = Struct_Layout{kkm,2} / Struct_Layout{kkm,3} ; %[um]
    
    % Propage through current zModule
    for kkk = 1 : Struct_Layout{kkm,3}
        
        kkg = kkg + 1; % increase global step-counter
        
        % Monitor the BPM progress
        if mod(kkg,round(StepsTotal*MonitorStep))==0,
            fprintf( ' -- %3d%% Completed @ z = %5.1f [um] // Step %4d of %4d :' , ...
                round( 100*( kkg/StepsTotal ) ) , zAx(kkg) , kkg , StepsTotal-1 ); tic;
        end
        
        % Calculate Propagation Matrices only for zVar-Module or 1st step of non-zVar
        if kkk<=2 || zVarFlags(kkm) == 1
            % Propagation Matrices (Tri-Diagonal) for Next and Curr step,
            % using given Pade-MultiStepping accuracy
            ns_curr = ns(:,kkg );  ns_next = ns(:,kkg+1);
            sx_curr = sx(:,kkg );  sx_next = sx(:,kkg+1);
            switch BPM
                case 0, BPMstr = 'Scalar';
                case 1, BPMstr = 'SemiVector';
            end
            [ NE1,CU1,NE2,CU2,NE3,CU3,NE4,CU4 ] = zBPM( BPMstr , MSP , aCN , wl , nref , ...
                dx , dz , ns_next , ns_curr , sx_next , sx_curr );
            clear BPMstr
        end
        
        % Do the z-Stepping (BPM)
        Ecurr = Fi( : , kkg );
        switch MSP % Multi-Step Parameters
            case 0, % Simplest Paraxial
                Enext = NE1 \ (CU1 * Ecurr);
            case 1, % Simple Wide-Angle
                Enext = NE1 \ (CU1 * Ecurr);
            case 2, % More Wide-Angle
                Ehalf = NE1 \ (CU1 * Ecurr);
                Enext = NE2 \ (CU2 * Ehalf);
            case 3, % Extra Wide-Angle
                Ehalf = NE1 \ (CU1 * Ecurr);
                Ehalf = NE2 \ (CU2 * Ehalf);
                Enext = NE3 \ (CU3 * Ehalf);
            case 4, % Ultra Wide-Angle
                Ehalf = NE1 \ (CU1 * Ecurr);
                Ehalf = NE2 \ (CU2 * Ehalf);
                Ehalf = NE3 \ (CU3 * Ehalf);
                Enext = NE4 \ (CU4 * Ehalf);
        end
        Fi( : , kkg+1 ) = Enext;
        clear Ecurr Enext;
        
        
        % Monitor BPM progress
        if mod( kkg , round(StepsTotal*MonitorStep) ) == 0,
            fprintf( ' %2d [msec]\n' , round(1e3*toc) );
        end
        
    end
    
end
fprintf( 'Done!\n' );

end


% =========================================================================
% Internal Function: BPM Sparse arrays
% =========================================================================
function [ NE1,CU1,NE2,CU2,NE3,CU3,NE4,CU4 ] = ...
    zBPM( BPM , MSP , aCN , wl , nref , dx , dz , ...
    ns_next , ns_curr , sx_next , sx_curr )

% This function returns the propagation matrices, Next and Current Step,
% of a 2D BPM (x- and z-axis) based on finite difference (FD) scheme. It
% also uses Crank-Nicolson stabilization across z-propagation and a
% multi-stepping [Hadley] algorithm for WA
%
% Alexandros Pitilakis
% Thessaloniki, June 2009
% Update: Oct 2015 -- Included PadeCoeffs as internal function

%Test Input Vars
if nargin == 0,
    BPM = 'Scalar'; MSP = 1; aCN = 0.5;
    wl = 1.55;     nref = 3.2;     dx = 0.1; dz = 0.5;
    ns_next = 3 + 0.5*rand(1,100)';   ns_curr = 3 + 0.5*rand(1,100)';
    sx_next = ones(size(ns_next)); sx_curr = sx_next;
end

%Restore Units to [m]
wl = wl*1e-6;
dx = dx*1e-6;
dz = dz*1e-6;

%Aux Variables:
k0 = 2*pi/wl;
kr = 2*pi/wl*nref;
Nx = length( ns_next );

% -------------------------------------------------------------------------
% Shift RefrInx & PML (sx) matrices, for each set of 3-neighbouring nodes
% -------------------------------------------------------------------------

%epsilon matrices (e=n^2), shifted +/- in the x-dim, for the next and curr step
es_next_pshift = [ ns_next(2:end) ; ns_next(end) ].^2;
es_next_mshift = [ ns_next(1) ; ns_next(1:end-1) ].^2;
es_curr_pshift = [ ns_curr(2:end) ; ns_curr(end) ].^2;
es_curr_mshift = [ ns_curr(1) ; ns_curr(1:end-1) ].^2;

%sx PML matrix, shifted likewise:
sx_next_pshift = [ sx_next(2:end) ; sx_next(end) ];
sx_next_mshift = [ sx_next(1) ; sx_next(1:end-1) ];
sx_curr_pshift = [ sx_curr(2:end) ; sx_curr(end) ];
sx_curr_mshift = [ sx_curr(1) ; sx_curr(1:end-1) ];

% -------------------------------------------------------------------------
% Elementary "Trans" and "Refl" coeffs for each inter-node-interface
% -------------------------------------------------------------------------

if strcmp( BPM , 'SemiVector' ) %Semi-Vectorial formulation (ri-discontinuities taken into account)
    
    %Trans Coeff, plus and mins (p,m), for the next and current step
    Tp_next = 2 * es_next_pshift ./ ( es_next_pshift + ns_next.^2 ) .* 2 ./ sx_next_pshift ./ ( sx_next_pshift + sx_next );
    Tp_curr = 2 * es_curr_pshift ./ ( es_curr_pshift + ns_curr.^2 ) .* 2 ./ sx_curr_pshift ./ ( sx_curr_pshift + sx_curr );
    Tm_next = 2 * es_next_mshift ./ ( es_next_mshift + ns_next.^2 ) .* 2 ./ sx_next_mshift ./ ( sx_next_mshift + sx_next );
    Tm_curr = 2 * es_curr_mshift ./ ( es_curr_mshift + ns_curr.^2 ) .* 2 ./ sx_curr_mshift ./ ( sx_curr_mshift + sx_curr );
    
    %Reflect Coeff, plus and minus (p,m), for the next and current step
    Rp_next = 1 + 2 * ns_next.^2 ./ ( es_next_pshift + ns_next.^2 ) .* 2 ./ sx_next ./ ( sx_next_pshift + sx_next );
    Rp_curr = 1 + 2 * ns_curr.^2 ./ ( es_curr_pshift + ns_curr.^2 ) .* 2 ./ sx_curr ./ ( sx_curr_pshift + sx_curr );
    Rm_next = 1 + 2 * ns_next.^2 ./ ( es_next_mshift + ns_next.^2 ) .* 2 ./ sx_next ./ ( sx_next_mshift + sx_next );
    Rm_curr = 1 + 2 * ns_curr.^2 ./ ( es_curr_mshift + ns_curr.^2 ) .* 2 ./ sx_curr ./ ( sx_curr_mshift + sx_curr );
    
elseif strcmp( BPM , 'Scalar' ) %Scalar Approach
    
    %Trans Coeff, plus and mins (p,m), for the next and current step
    Tp_next = 2 ./ sx_next_pshift ./ ( sx_next_pshift + sx_next );
    Tp_curr = 2 ./ sx_curr_pshift ./ ( sx_curr_pshift + sx_curr );
    Tm_next = 2 ./ sx_next_mshift ./ ( sx_next_mshift + sx_next );
    Tm_curr = 2 ./ sx_curr_mshift ./ ( sx_curr_mshift + sx_curr );
    
    %Reflect Coeff, plus and minus (p,m), for the next and current step
    Rp_next = 1 + 2 ./ sx_next ./ ( sx_next_pshift + sx_next );
    Rp_curr = 1 + 2 ./ sx_curr ./ ( sx_curr_pshift + sx_curr );
    Rm_next = 1 + 2 ./ sx_next ./ ( sx_next_mshift + sx_next );
    Rm_curr = 1 + 2 ./ sx_curr ./ ( sx_curr_mshift + sx_curr );
    
end

% -------------------------------------------------------------------------
% Propagator Matrices // Multi-Stepping
% -------------------------------------------------------------------------

%Calculate Pade-Coefficient values for Multi-Stepping Algorightm
bCN = 1-aCN;
[aN, aD] = Pade_Coeffs( MSP , kr , aCN , bCN , abs(dz) ); %abs(dz)->Imaginary Distance BPM

%Init Propagation matrices, up to Pade(3,3)
NE1 = []; CU1 = []; NE2 = []; CU2 = []; NE3 = []; CU3 = []; NE4 = []; CU4 = [];

for kms = 1 : MSP
    
    % -------------------------------------------------------------------------
    % Form the 3 Diagonals-Vectors of the Sparse Matrices:
    % -------------------------------------------------------------------------
    
    %Assign the proper values for the mutli-stepping params:
    constN = aN(kms); %Nominator   -> Curr
    constD = aD(kms); %DeNominator -> Next
    % constN =  ( 1 - j * dz * kr ) / ( 4 * kr^2 ); %Nominator   -> Curr
    % constD =  ( 1 + j * dz * kr ) / ( 4 * kr^2 ); %DeNominator -> Next
    
    %Assign the next/curr index-step values appropriately
    switch MSP
        case 0, %Simplest, Paraxial
            Rpn=Rp_next; Rmn=Rm_next; nsn=ns_next; Tpn=Tp_next; Tmn=Tm_next;
            Rpc=Rp_curr; Rmc=Rm_curr; nsc=ns_curr; Tpc=Tp_curr; Tmc=Tm_curr;
        case 1, %Pade(1,1), Simple Wide Angle
            Rpn=Rp_next; Rmn=Rm_next; nsn=ns_next; Tpn=Tp_next; Tmn=Tm_next;
            Rpc=Rp_curr; Rmc=Rm_curr; nsc=ns_curr; Tpc=Tp_curr; Tmc=Tm_curr;
        case 2, %Pade(2,2), More Wide Angle
            switch kms
                case 1,
                    Rpn=Rp_curr; Rmn=Rm_curr; nsn=ns_curr; Tpn=Tp_curr; Tmn=Tm_curr;
                    Rpc=Rp_curr; Rmc=Rm_curr; nsc=ns_curr; Tpc=Tp_curr; Tmc=Tm_curr;
                case 2,
                    Rpn=Rp_next; Rmn=Rm_next; nsn=ns_next; Tpn=Tp_next; Tmn=Tm_next;
                    Rpc=Rp_next; Rmc=Rm_next; nsc=ns_next; Tpc=Tp_next; Tmc=Tm_next;
            end
        case 3, %Pade(3,3), Extra Wide Angle
            switch kms
                case 1,
                    Rpn=Rp_curr; Rmn=Rm_curr; nsn=ns_curr; Tpn=Tp_curr; Tmn=Tm_curr;
                    Rpc=Rp_curr; Rmc=Rm_curr; nsc=ns_curr; Tpc=Tp_curr; Tmc=Tm_curr;
                case 2,
                    Rpn=Rp_next; Rmn=Rm_next; nsn=ns_next; Tpn=Tp_next; Tmn=Tm_next;
                    Rpc=Rp_curr; Rmc=Rm_curr; nsc=ns_curr; Tpc=Tp_curr; Tmc=Tm_curr;
                case 3,
                    Rpn=Rp_next; Rmn=Rm_next; nsn=ns_next; Tpn=Tp_next; Tmn=Tm_next;
                    Rpc=Rp_next; Rmc=Rm_next; nsc=ns_next; Tpc=Tp_next; Tmc=Tm_next;
            end
        case 4, %Pade(4,4), Ultra Wide Angle
            switch kms
                case 1,
                    Rpn=Rp_curr; Rmn=Rm_curr; nsn=ns_curr; Tpn=Tp_curr; Tmn=Tm_curr;
                    Rpc=Rp_curr; Rmc=Rm_curr; nsc=ns_curr; Tpc=Tp_curr; Tmc=Tm_curr;
                case 2,
                    Rpn=Rp_curr; Rmn=Rm_curr; nsn=ns_curr; Tpn=Tp_curr; Tmn=Tm_curr;
                    Rpc=Rp_curr; Rmc=Rm_curr; nsc=ns_curr; Tpc=Tp_curr; Tmc=Tm_curr;
                case 3,
                    Rpn=Rp_next; Rmn=Rm_next; nsn=ns_next; Tpn=Tp_next; Tmn=Tm_next;
                    Rpc=Rp_next; Rmc=Rm_next; nsc=ns_next; Tpc=Tp_next; Tmc=Tm_next;
                case 4,
                    Rpn=Rp_next; Rmn=Rm_next; nsn=ns_next; Tpn=Tp_next; Tmn=Tm_next;
                    Rpc=Rp_next; Rmc=Rm_next; nsc=ns_next; Tpc=Tp_next; Tmc=Tm_next;
            end
    end
    
    
    next_Diag0 = ( (2-Rpn-Rmn)/dx^2  + (k0*nsn).^2 - kr^2  ) * constD + 1;
    next_Diagp = Tpn/dx^2 * constD ;
    next_Diagm = Tmn/dx^2 * constD ;
    
    curr_Diag0 = ( (2-Rpc-Rmc)/dx^2  + (k0*nsc).^2 - kr^2  ) * constN + 1;
    curr_Diagp = Tpc/dx^2 * constN ;
    curr_Diagm = Tmc/dx^2 * constN ;
    
    % -------------------------------------------------------------------------
    % Assemble Sparse-Matrices
    % -------------------------------------------------------------------------
    
    %Some Minor Correction for the p/m matrices: A shifting to work w/ spdiags:
    next_Diagp = circshift( next_Diagp , +1 );
    curr_Diagp = circshift( curr_Diagp , +1 );
    next_Diagm = circshift( next_Diagm , -1 );
    curr_Diagm = circshift( curr_Diagm , -1 );
    
    %Form the Tri-Diagonal Sparse Propagation Matrices:
    N_Matrix = spdiags( [ next_Diagm next_Diag0 next_Diagp ] , [-1 0 +1] , Nx , Nx );
    C_Matrix = spdiags( [ curr_Diagm curr_Diag0 curr_Diagp ] , [-1 0 +1] , Nx , Nx );
    
    %Assign the final propagator matrices (next/curr)
    %for each "sub-step"
    switch kms
        case 1, %Pade(1,1)
            NE1 = N_Matrix;
            CU1 = C_Matrix;
        case 2, %Pade(2,2)
            NE2 = N_Matrix;
            CU2 = C_Matrix;
        case 3, %Pade(3,3)
            NE3 = N_Matrix;
            CU3 = C_Matrix;
        case 4, %Pade(4,4)
            NE4 = N_Matrix;
            CU4 = C_Matrix;
    end
    clear N_Matrix C_Matrix
    
end

% -------------------------------------------------------------------------
% Paraxial Approximation
% -------------------------------------------------------------------------
if MSP == 0
    
    %Assign the proper values for the mutli-stepping params:
    constN = +bCN*dz/(2*1j*kr) ; %Nominator   -> Curr
    constD = -aCN*dz/(2*1j*kr) ; %DeNominator -> Next
    
    next_Diag0 = ( (2-Rp_next-Rm_next)/dx^2  + (k0*ns_next).^2 - kr^2  ) * constD + 1;
    next_Diagp = Tp_next/dx^2 * constD ;
    next_Diagm = Tm_next/dx^2 * constD ;
    
    curr_Diag0 = ( (2-Rp_curr-Rm_curr)/dx^2  + (k0*ns_curr).^2 - kr^2  ) * constN + 1;
    curr_Diagp = Tp_curr/dx^2 * constN ;
    curr_Diagm = Tm_curr/dx^2 * constN ;
    
    %Some Minor Correction for the p/m matrices: A shifting to work w/ spdiags:
    next_Diagp = circshift( next_Diagp , +1 );
    curr_Diagp = circshift( curr_Diagp , +1 );
    next_Diagm = circshift( next_Diagm , -1 );
    curr_Diagm = circshift( curr_Diagm , -1 );
    
    %Assemble the Sparse Tri-Diagonal Propagation Matrices:
    NE1 = spdiags( [ next_Diagm next_Diag0 next_Diagp ] , [-1 0 +1] , Nx , Nx );
    CU1 = spdiags( [ curr_Diagm curr_Diag0 curr_Diagp ] , [-1 0 +1] , Nx , Nx );
    
end


%Test Plot:
if nargin == 0, spy( CU1 ); end

end


% =========================================================================
% Internal Function: Pade Coefficients
% =========================================================================
function [aN, aD] = Pade_Coeffs( MSP , kr , aCN , bCN , dz )

% Calculates the Pade(MSP,MSP) coefficients with a recursive procedure,
% for use in a Multi-Step Beam-Propagation Method to introduce a
% wide-angle approximation of the Propagator (2nd dz deriv).
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2010 June

%Assign the proper polynomials of "P" operator (= A^-1 * B , see my notes)
%for the requested Pade approximation, according to the recursive
%approximation of the dz^2 derivative:
%                     diff(F)/dz = -j*N/D * F ,
%where N,D stand for the Nominator/Denominator polynomials of "P" and
%F is the (transformed) field envelope
% switch MSP
%     case 1,
%         NomPoly = [ 2*kr 0      ] ;
%         DeNPoly = [ 1    4*kr^2 ] ;
%     case 2,
%         NomPoly = [ 4*kr 8*kr^3  0       ] ;
%         DeNPoly = [ 1    12*kr^2 16*kr^4 ] ;
%     case 3,
%         NomPoly = [ 6*kr 32*kr^3 32*kr^5 0       ] ;
%         DeNPoly = [ 1    24*kr^2 80*kr^4 64*kr^6 ] ;
% end

%Recursive Calculation with the relationship:
%d/dz(n+1) * F = -P/(  d/dz(n)  - j*2*kr   ) * F
%resulting in: N(n+1)=P*D(n), D(n+1)=N(n)+2*kr*D(n)

N = 0; D = 1;%Initial Operator d/dz=0
for n = 1 : (2*MSP)
    
    %N,D Polynomial-sizes must be the same (N increases before D)
    if length( D ) < length(N), auxD = [ 0 D ]; else  auxD = D; end
    
    %Recursive Algorithm
    NN = [ D 0 ];
    DD = N + 2*kr*auxD ;
    
    %It is customary to normalize all coeffs of the polynomials with
    %max coeff of D (odd iteration) or N (even iteration)
    if mod(n,2) == 1, auxDiv = DD(1); else auxDiv = NN(1); end
    N = NN / auxDiv;
    D = DD / auxDiv;
    
end
NomPoly = N;
DeNPoly = D;

%Find the Coefficients needed to factorize the "propagator", after the
%application of the CN-scheme and the step-dz-derivative:
%              (Fn-Fc)/dz = j*N/D*( aCN*Fn + bCN*Fc ) =>
%              Fn = (D-j*aCN*dz*N)/(D+j*bCN*dz*N) * Fc
%where Fn,Fc are the field cross-sections at the next (unknown) and
%current (known) z-step. aCN>=0.5 for stability
aN = -1./roots( DeNPoly - 1j*dz*bCN* NomPoly ); %Nominator-> "known field"
aD = -1./roots( DeNPoly + 1j*dz*aCN* NomPoly ); %DeNominator-> "unknown field"

%Sort them in some decent way (roots returns them abs-wise sorted)
re = real( aN ); [~,i] = sort( re ); aN = aN(i);
re = real( aD ); [~,i] = sort( re ); aD = aD(i);


% % ......................................................................
% % Obsolte PML Implementation
% % ......................................................................
% %"Add" the PML tensor to the Refr.Indx. It can be easily separated afterwards
% ns = ns .* sx;
% %Separate PML from Refr.Indx. Profile:
% sx_next = 1 + j*imag(ns_next)./real(ns_next);
% sx_curr = 1 + j*imag(ns_curr)./real(ns_curr);
%
% %Keep only the real part in ns_next/curr:
% ns_next = real(ns_next);
% ns_curr = real(ns_curr);
end


