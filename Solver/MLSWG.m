function [ neffs , Fy , x_out ] = MLSWG( ModePol , wl , nLR , ns , ts , ...
    x , neSL , DoPlotCharEq , DoPlotModeProfiles , DoVocalize )

% MultiLayer Slab WaveGuide (MLSWG) solver.
%
%     Used to find the modes of MLSWGs with arbitrary number of guiding 
%     layers, that can have arbitrary thicknesses and complex refractive 
%     indices. This function can find the guided modes by solving the 
%     characteristic equation. It might not find all modes at once, but 
%     repeated calls with different neff search limits (neSL) and can
%     make it happen.
%
% ==== Outputs ====
%  - neffs : vector containing the effective refr indices of modes found
%  - Fy    : corresponding mode profiles (Ey or Hy, for TE and TM modes)
%  - x_out : cross-section space, used to plot Fy (if input x is not given)
%
% ==== Inputs / Obligatory ====
%  - ModePol : 'TE' or 'TM' (string) - Polarization of modes
%  - wl  : wavelength [Same units as ts!]
%  - nLR : refr. indices of [Left,Right] semi-inf layers (2x1 array)
%  - ns  : refr. indices of guiding layers (NLx1 array)
%  - ts  : thicknesses of guiding layers (NLx1 array) [Same units as wl!]
%
% ==== Inputs / Optional ====
%  - neSL: n-effective search limits (2x1 array, [Low,High])
%  - x   : WG 1D cross-section (Nx-length vector) [Same units as wl & ts!]
%
% ==== Inputs / Monitoring (optional) ====
%  - DoPlotCharEq       : 1 or 0 --> Plot Characteristic Equation
%  - DoPlotModeProfiles : 1 or 0 --> Plot mode field profiles found
%  - DoVocalize         : 1 or 0 --> Disp info on MATLAB's command window
%
%  * Note: 1st guiding-layer extends from x=[0,ts(1)], 
%          2nd from x=ts(1)+[0,ts(2)] and so on.
%
%  * Note: Feel free to change the parameters inside this function and
%          maybe activate other available solvers, like:
%     1. FIDODIES: Numerical Finite-Difference Eigenmode solver (SPTARN)
%     2. INTERIPV: Graphical solution to the Characteristic-Equation
%     3. FZERO   : MATLAB's solvers
%
% Alexandros Pitilakis / Thessaloniki, Greece
%  2015 Nov: Original version

% PRIVATE FUNCTION: Post a message in MATLAB's command window...
    function postMsg( msgString )
        if DoVocalize == 1 % ...controlled by this variable
            fprintf( msgString );
        end
    end

% Optional Inputs
if nargin < 6 , x    = []; end % for field's calculation
if nargin < 7 , neSL = []; end % neff search limits (speed-up!)

% Monitoring params
if nargin < 8 , DoPlotCharEq = 0; end % Plot CharEq
if nargin < 9 , DoPlotModeProfiles  = 0; end % Plots Fields of modes found
if nargin < 10, DoVocalize   = 0; end % Show stuff on command window

% Test inputs
if nargin == 0
        
    ModePol = 'TE'; % Mode-Polarization    
    wl = 1.55; % Wavelength [um]
        
    % Indices (guiding layers are from L->R)
    nLR = [ 1.45, 1.45 ]; % Indices of L/R layers
    Nbr = 2; % number of "core" waveguides (1=Single, 2=Coupler, 3=Tripler etc)
    nco = 3.2 - 0.00j; % core index
    nbu = nLR(1); % buffer (between cores) index
    ns = [ nco , repmat( [nbu nco] , 1 , Nbr-1 ) ];
    
    % Thicknesses of guiding-layers (from L->R)
    wid = 0.25; % core-widths
    gap = 0.42; % gap between cores (buffer-filled)
    ts = [ wid , repmat( [gap wid] , 1 , Nbr-1 ) ]; % [um]
    
    % Perturbations
    if length(ns) > 1 && Nbr >= 2
        ns(3) = ns(3) - 0.025*1;
    end 
    
    % Monitoring params
    DoPlotCharEq = 1; % Plot CharEq
    DoPlotModeProfiles  = 1; % Plots Fields of modes found
    DoVocalize   = 1; % Show stuff on command window

end

% Which n-effective set to use for mode-profile calculation & plotting
whichNeffsToOutput = 0; % 0=myNewtonRaphson, 1=FIDODIES, 2=INTERPINV, 3=FZERO 

% -------------------------------------------------------------------------
% Main Solver + params
% -------------------------------------------------------------------------

% Root-Finding-Algo (RFA): {0=BiSection,1=Secant,2=RegulaFalsi,3=NR,4=Advanced-NR}
RFA = 4; % Default =4

% Method Parameters
roundTol = 1e-6; % <<1 Round-off tolerance (for comparing roots)
Tol_neff = roundTol; % tol/accuracy in neff=root (horiz-axis)
Tol_XE   = roundTol; % tol/accuracy in XE(neff)=0 (vert-axis)
MaxIter  = 50 ; % max iterations
ccSlope  = 3; % Controlled Convergence of NR (>1: slow-fast / <1: fast-slow / ==0: flat)

% -------------------------------------------------------------------------
% Aux Solvers + Params
% -------------------------------------------------------------------------

% Alt method: Eigenmode FD solution (FIDODIES) + params ----> TE modes!
useAlsoFIDODIES = 0 + whichNeffsToOutput==1; 
xL = max(ts); % Size of semi-inf space (Left)
xR = max(ts); % Size of semi-inf space (Right)
Nx = 1e4; % Number of descretization points

% Alt method: Graphical (INTERPINV) + params ----> Uses neSL !
useAlsoINTERPIV = 0 + whichNeffsToOutput==2;  
nNeff = 1e5; % Number of neff-samples for graphical solution
tanTol = 1e0; % Tangent-like discontinuity (at pi/2) tolerance.

% Alt method: MATLAB's own FZERO (or FSOLVE) ----> Uses neSL !
useAlsoFZERO = 0 + whichNeffsToOutput==3; 
Nx0 = 50; % steps of starting-guesses (in neSL range)
Tol_XE_FZ = 1e-8; % f(x) tolerance (we want this close to zero)
roundTol_FZ = 1e-8; % round-off tolerance

% -------------------------------------------------------------------------
% Preparations
% -------------------------------------------------------------------------

% x-space for calc/plot fields of modes found
if isempty(x)
    subsSiz = 1; % substrate size, in multiples of max(ts)
    x = linspace( -subsSiz*max(ts) , sum(ts)+subsSiz*max(ts) , 1e3 );
end

% Return x-vector, for plotting the fields externally
if nargout == 3,
    x_out = x;
end

% Acquire input-params and assign them to local variables
nL = nLR(1);
nR = nLR(2);

% The min/max limits of the valid neff-range of guided modes in
% photonic AND plasmonic waveguides, lossy or lossless.
neVL = [ max([nL nR]) max(ns) ] + 1e-10 * [ +1 -1 ];
neVL = sort( sqrt( abs( real( neVL.^2 ) ) ) );

% n-effective search limits
if isempty(neSL)
    neSL = neVL;
end

% -------------------------------------------------------------------------
% Solve w FIDODIES
% -------------------------------------------------------------------------
if useAlsoFIDODIES == 1
    
    postMsg( sprintf('%s\n FIDODIES: Numerical FD Eigenmode Solver (SPTARN) \n%s\n', ...
        flwcs('-') , flwcs('-') ) );
    
    if strcmp( ModePol , 'TE' )
        
        % x-vector (FDM space) -- The denser, the higher the accuracy
        xF = linspace( -xL , sum(ts)+xR , Nx );
        
        % Form: PMLs tensors (sx) at each x-point
        sx = ones(size(xF)); % initialize -- no anisotropy/absorption
        sPML = 1; % strength of PML
        wPML = (max(xF)-min(xF))/10; % thickness of PML
        if sPML ~= 0
            iL = xF<=(xF( 1 )+wPML); % indices of Left PML
            iR = xF>=(xF(end)-wPML); % indices of Right PML
            sx( iL ) = sx( iL ) - 1j*sPML*( ( sum(iL==1):-1:1 )/sum(iL==1) ).^2 ; % Parabolic profile
            sx( iR ) = sx( iR ) - 1j*sPML*( ( 1:+1:sum(iR==1) )/sum(iR==1) ).^2 ; % Parabolic profile
        end
        
        % Form: Index-profile (nsF) at each  x-point
        nsF = -1 * ones(size(xF)); % initialize -- all at n_substrate
        hn = [0 cumsum(ts)];
        nsF(  xF<=0 ) = nL;
        nsF ( xF>hn(end) ) = nR;
        for jj = 1:length(ns)
            nsF( xF>hn(jj) & xF<=hn(jj+1) ) = ns(jj);
        end
        
        % Call FIDODIES solver (SPTARN-based)
        neffsEI = FIDODIESv2( wl , xF , nsF , sx );
        postMsg( sprintf( ' >> Roots = %12.10f (real)\n' , real(neffsEI) ) );
        
    else
        postMsg( ' ## FIDODIES: Can''t Solve for TM modes... Sorry!' );
    end
    
end

% -------------------------------------------------------------------------
% Solve w INTERPOLATION
% -------------------------------------------------------------------------
if useAlsoINTERPIV == 1
    
    % Solver -- Neff range, and samples in that range. Sometimes two modes
    % can be very close (almost degenerate, e.g. at couplers) which makes
    % discerning them quite hard.
    
    postMsg( sprintf('%s\n INTERPINV: Graphical Solution to XE \n%s\n', ...
        flwcs('-') , flwcs('-') ) );
    
    % Candidate n_eff range:
    ne = linspace( neSL(1)+eps , neSL(2)-eps , nNeff );
    
    % Get values of CharEq at ne-space
    XE = MLSWG_CharEq( ne , nL , nR , ns , wl , ts , ModePol );
    
    % Solve graphically, with interpolation
    neffsIP = flipud( interpinv( ne , real(XE) , 0 , tanTol ) );
    
    postMsg( sprintf( ' >> Roots = %12.10f (real)\n' , real(neffsIP) ) );
    
end

% -------------------------------------------------------------------------
% Solve w MATLAB's own FZERO
% -------------------------------------------------------------------------
if useAlsoFZERO == 1
    
    postMsg( sprintf('%s\n MATLAB''s F-ZERO: \n%s\n', ...
        flwcs('-') , flwcs('-') ) );
    
    % Solve with FZERO --> single point
    neStart = linspace( neSL(1)+eps , neSL(2)-eps , Nx0 );
    for kk = 1:length(neStart) % multiple runs w diff starting guess, to find multiple roots
        [ ne0(kk) fv(kk) errFlag(kk)] = fzero( ...
            @(ne) MLSWG_CharEq( ne, nL, nR, ns, wl, ts, ModePol, true ) , neStart(kk) );
    end
    
    ne0 = unique( ne0( abs(fv) < Tol_XE_FZ & errFlag==1 ) ); % keep only roots (|fv|<<Tol, exitFlag=1)
    ne0 = unique( round( ne0 / roundTol_FZ  ) * roundTol_FZ  ); % round-off to 10'th decimal
    neffsFZ = ne0( ne0 > neStart(1) & ne0 < neStart(end) ); % keep in-range roots
    
    postMsg( sprintf( ' >> Roots = %12.10f (real)\n' , real(neffsFZ) ) );
    
end

% -------------------------------------------------------------------------
% Solve w BISECTION / SECANT / REGULA-FALSI / NEWTON-RAPHSON
% -------------------------------------------------------------------------

% ** Infinities-search (NUMERICAL), based on diff(XE) function might help.

postMsg( sprintf('%s\n SOLVE: Characteristic Equation \n%s\n', flwcs('-') , flwcs('-') ) );

% functionHandle (for NR), that returns both y|@ne0 and dy/dx|@ne0
fH = @(ne) MLSWG_CharEq(ne,nL,nR,ns,wl,ts,ModePol);

% Root-Finding-Algo (RFA) Used
RFA_Name = {'BiSection','Secant','Regula Falsi','Newton-Raphson','Advanced Newton-Raphson'};
postMsg( sprintf(' ** Method used: %s\n', RFA_Name{RFA+1} ) );

% --- BISECTION / SECANT / REGULA-FALSI ---
if any( RFA == [0,1,2] )
    
    % Runs from two designated start-points (in ne-Limits)
    a = neSL(1); % 1st of the two extremities of the neff-range
    b = neSL(2); % 2nd of the two extremities of the neff-range
    xR = []; % this will hold the root
    
    % Iterate
    for knr = 1 : MaxIter
        
        % Values of functions f(x) at points a & b
        fa = MLSWG_CharEq( a , nL , nR , ns , wl , ts , ModePol );
        fb = MLSWG_CharEq( b , nL , nR , ns , wl , ts , ModePol );
        
        % Special conditions for BiSection & RegulaFalsi methods.
        if any( RFA == [0,2] )
            
            if any( abs(imag([fa,fb])) > roundTol )
                error( ' ## Bisection / Regula Falsi: Only for REAL functions!' );
            end
            
            if sign(fa) == sign(fb)
                error( ' ## Bisection / Regula Falsi: Starting points f(a) and f(b) must have opposite sign' );
            end
            
        end
        
        % Next root-estimate via secant or midpoint (for Bisection)
        if RFA==1 || RFA==2
            x0 = b - fb*(b-a)/(fb-fa);
        else
            x0 = (a+b)/2;
        end
        postMsg( sprintf( ' Iteration %3d --> x0 = %1.10f\n' , knr , x0 ) );
        
        % Calc funct value f(x0) at new estimate point
        fx0 = MLSWG_CharEq( x0 , nL , nR , ns , wl , ts , ModePol );
        
        % Checks for:
        if abs(x0-b) < Tol_neff && abs(fx0) < Tol_XE ... % Root-found?
                && real(x0) > min(neSL) && real(x0) < max(neSL)
            xR = [ xR , x0 ];
            break;
        elseif abs(a-b) < Tol_neff && abs(fx0) > Tol_XE... % INF convergence?
                && real(x0) > min(neSL) && real(x0) < max(neSL)
            postMsg( sprintf( ' ## Secant/RegulaFalsi: Convergence at INF?\n' ) );
            break;
        elseif knr == MaxIter, % iterations exceeded?
            postMsg( sprintf( ' ## Secant/RegulaFalsi: MaxIterations Exceeded.\n' ) );
        else % prepare for next iteration
            
            if RFA == 1, % Secant Method
                a = b;
                b = x0;
            else % Regula-Falsi / Bisection
                if sign(real(fa)) == sign(real(fx0))
                    a = x0;
                else
                    b = x0;
                end
            end
            
        end
    end
    
    % Round-off and concat roots
    xR = round(1/roundTol*xR)*roundTol;
    neffsNR = fliplr( unique(xR) );
    
end

% --- NEWTON-RAPHSON (Basic): Multiple evenly spaced guesses ---
if RFA==3
    % Run NR method from various starting estimates
    a = neSL(1); % 1st of the two extremities of the neff-range
    b = neSL(2); % 2nd of the two extremities of the neff-range
    xR = []; % this will hold the roots
    for x0 = linspace( a, b, 1 )
        
        for knr = 1 : MaxIter
            [ XEa, dXEa ] = MLSWG_CharEq( x0 , nL , nR , ns , wl , ts , ModePol );
            
            om = (knr/MaxIter)^ccSlope ; % increases from 0->1 (for faster convergence) as knr approaches MaxIter
            x0 = x0 + om*( -XEa./dXEa );
            
            postMsg( sprintf( ' Iteration %3d --> x0 = %f\n' , knr , x0 ) );
            
            if abs(XEa./dXEa) < Tol_neff && abs(XEa)< Tol_XE ...
                    && real(x0) > min(neSL) && real(x0) < max(neSL)
                xR = [ xR , x0 ];
                break;
            elseif real(x0) < min(neSL) || real(x0) > max(neSL)
                postMsg( sprintf( ' ## Newton-Raphson: Solver outside neff-range. Aborting.\n' ) );
                break;
            elseif knr == MaxIter,
                postMsg( sprintf( ' ## Newton-Raphson: MaxIterations Exceeded.\n' ) );
            end
        end
        
        % Or, replacing the above for-loop:
        % xR = myNewtonRaphson( fH, x0, neL, Tol_neff , Tol_XE , MaxIter , ccSlope , 0 );
        
    end
    
    % Compare roots, and decide if you should continue:
    xR = round(1/roundTol*xR)*roundTol;
    neffsNR = fliplr( unique(xR) );
end

% --- NEWTON-RAPHSON (Advanced): Leapfrog between roots (guesses) ---
if RFA==4
    
    DoProceed = 1; % This contols the "solver"
    
    % First NR, from lowest possible neff
    postMsg( '  = RUN#1: From lowest possible neff.\n' );
    xR = myNewtonRaphson( fH, neSL(1), neSL, Tol_neff , Tol_XE , MaxIter , ccSlope , 0 );
    
    % Check if a valid root was found. If not, don't despair (yet).
    xR = round(1/roundTol*xR)*roundTol;
    if isempty(xR),
        postMsg( '  * Nothing found from here.\n' );
    else
        postMsg( sprintf( '  * Root Found = %12.10f\n' , xR ) );
    end
    
    % Second NR, from highest possible neff
    postMsg( '  = RUN#2: From highest possible neff.\n' );
    xRc = myNewtonRaphson( fH, neSL(2), neSL, Tol_neff , Tol_XE , MaxIter , ccSlope , 0 );
    
    % Check the root found. If it's equal to previous, quit solver
    xRc = round(1/roundTol*xRc)*roundTol;
    if isempty(xRc),
        postMsg( '  * RUN#2: No root found... Strange!\n' );
        DoProceed = 0;
    elseif xRc == xR
        postMsg( '  * RUN#2: Root identical to 1st run ==> Single-root.\n' );
        DoProceed = 0;
    else
        postMsg( sprintf( '  > Root Found = %12.10f\n' , xRc ) );
        xR = [xR xRc];
    end
    
    % Now, start "leapfrogging" around all the means
    newRoots = xR;
    cc = 0;
    while DoProceed == 1
        
        APs = unique( [ neSL , newRoots ] ); % all points
        [a,b]=meshgrid(APs,APs); r=rand; c=(1-r)*a+r*b;
        c=unique(c(:));
        c = min(max(c,neSL(1)),neSL(2));
        AMs=setdiff(c,APs); % all means
        
        cc = cc + 1;
        postMsg( sprintf( '  = RUN#3.%d: Leap-frogging around (size=%d)\n' , cc , length(AMs) ) );
        
        % Scan all-means and re-fill newRoots
        newRoots = [];
        for ne0 = AMs(:)'
            
            xRc = myNewtonRaphson( fH, ne0, neSL, Tol_neff , Tol_XE , MaxIter , ccSlope , 0 );
            newRoots = [newRoots, xRc];
            
        end
        newRoots = round(1/roundTol*newRoots)*roundTol;
        
        % Check if any there's any root in newRoots that is not in xR
        if ~isempty( setdiff( newRoots , xR ) )
            postMsg( sprintf( '  > New Root Found = %12.10f\n' , setdiff( newRoots , xR ) ) );
            xR = unique( [xR,newRoots] ); % merge the two
            newRoots = xR; % replace, for next iteration
        else
            DoProceed = 0;
            postMsg( ' ** END: No More New Roots Found. Completed.\n' );
        end
        
    end
    
    % Report results
    neffsNR = fliplr(xR);
    
    
end

% Show roots found
postMsg( sprintf( '%s\n' , flwcs('-') ) )
postMsg( sprintf( ' >> Roots = %12.10f (real)\n' , real(neffsNR) ) );
postMsg( sprintf( '%s\n' , flwcs('-') ) )
    
% Which neff set to out?
switch whichNeffsToOutput
    case 0, neffs = neffsNR(:); % my NewtonRaphson
    case 1, neffs = neffsEI(:); % FIDODIES (eigenmode)
    case 2, neffs = neffsIP(:); % InterPinv "graphical" method
    case 3, neffs = neffsFZ(:); % MATLAB's fzero
end

% -------------------------------------------------------------------------
% Plot: XE & Derivative-of-XE (XE=CharEq) with Roots found
% -------------------------------------------------------------------------
if DoPlotCharEq == 1
    
    maxXE = 1e3;
    ne = linspace( neSL(1) , neSL(2) , 1e5 );
    [ XE, dXE ] = MLSWG_CharEq( ne , nL , nR , ns , wl , ts , ModePol );
    
    Xfoc =  neSL; % "focus" plot on  neff within the search-limits
    
    subplot(2,1,1)
    plot( ne , real(XE) , 'bo' ); hold on;
    plot( ne , imag(XE) , 'b:' );
    plot( ne , 0*ne , 'k:' );
    set(gca,'YLim',maxXE*[-1 +1]);
    set(gca,'XLim',Xfoc);
    for kr = 1:length(xR);
        plot( real(xR(kr))*[1 1] , get(gca,'YLim') , 'b-' )
    end
    xlabel( 'n_{eff}' ); ylabel( 'Char. Eq. (XE)' )
    legend( 'Real' , 'Imag' )
    
    subplot(2,1,2)
    plot( ne , real(dXE) , 'ro' ); hold on;
    plot( ne , imag(dXE) , 'r:' );
    plot( ne , 0*ne , 'k:' );
    set(gca,'YLim',maxXE*[-1 0]);
    set(gca,'XLim',Xfoc);
    for kr = 1:length(xR);
        plot( real(xR(kr))*[1 1] , get(gca,'YLim') , 'r-' )
    end
    xlabel( 'n_{eff}' ); ylabel( 'dXE/dn_{eff}' )
    legend( 'Real' , 'Imag' )
    
end

% -------------------------------------------------------------------------
% Calc: Mode Profiles (knowing the neffs)
% -------------------------------------------------------------------------
if nargout >= 2 ||  DoPlotModeProfiles == 1
    
    % Initialize matrix to hold profiles
    Fy = NaN*zeros( length(neffs) , length(x) ); % main field / transverse / parallel to interfaces
    Gx = NaN*zeros( length(neffs) , length(x) ); % secondary field / trans / normal to interfaces
    Gz = NaN*zeros( length(neffs) , length(x) ); % secondary field / longitudinal
    Pz = NaN*zeros( length(neffs) , length(x) ); % power-flux
    
    % Polarization-dependent Factors
    if strcmp( ModePol , 'TM' )
        pdFL = ns( 1 )^2 / nL^2;
        %pdFR = ns(end)^2 / nR^2;
        pdFn = ones(size(ns));
        for jj = 2 : length(ns)
            pdFn(jj) = ns(jj)^2 / ns(jj-1)^2;
        end
    else
        pdFL = 1;
        %pdFR = 1;
        pdFn = ones(size(ns));
    end
    hn = cumsum(ts);
    
    % Dielectric constant er=n^2 at each x-point
    er = nR^2 * ones( size(x) );
    for jj = length(ns):-1:1
        er( x<hn(jj) ) = ns(jj)^2;
    end
    er( x<0 ) = nL^2;
    
    % Scan the modes found
    for m = 1:length(neffs)
        
        % Take it's n-effective
        neff = neffs(m);
        
        % L/R layer variables (now they are scalars)
        gL = 2*pi/wl * sqrt( neff.^2 - nL^2 );
        gR = 2*pi/wl * sqrt( neff.^2 - nR^2 );
        
        % Intermediate layer variables (now they are scalars)
        kn = NaN*zeros(size(ns)); % initialize!
        for ii = 1 : length(ns)
            kn(ii) = 2*pi/wl * sqrt( ns(ii)^2 - neff^2 );
        end
        
        % Layer-amplitudes (A) and auxiliary "angles" (the)
        A   = NaN*zeros(size(ns)); % Amplitudes
        the = NaN*zeros(size(ns)); % Angles
        A(1) = 1; % Scaling factor: Field in 1st intermedia layer
        the(1) = atan( gL ./ kn(1) * pdFL ); % This one is easy!
        AL = A(1)*cos(-the(1)); % Left-side substrate
        if length(ns) >= 2
            for jj = 2 : length(ns)
                auxi = tan( kn(jj-1)*hn(jj-1) - the(jj-1) );
                the(jj) = + kn(jj)*hn(jj-1) - atan( pdFn(jj) * kn(jj-1)./kn(jj) .* auxi );
                A(jj) = A(jj-1)*cos( kn(jj-1)*hn(jj-1) - the(jj-1) )/cos( kn(jj)*hn(jj-1) - the(jj) );
            end
        end
        AR = A(end)*cos( kn(end)*hn(end) - the(end) ); % Right-side substrate
        
        % Now, calc field-profiles (for each layer) along the x-vector
        hs = [0 hn]; % locations of all the interfaces
        iL = x<0;        Fy( m , iL ) = AL * exp( +gL*( x(iL)-hs(1)   ) ); % Left substrate
        iR = x>=hn(end); Fy( m , iR ) = AR * exp( -gR*( x(iR)-hs(end) ) ); % Right substrate
        for jj = 1:length(ns)
            ix = x>=hs(jj) & x<hs(jj+1);
            Fy( m , ix ) = A(jj) * cos( kn(jj)*x(ix) - the(jj) );
        end
        
        % Calc the "other" field (e.g. H/E if you're TE and TM)
        k0 = 2*pi/wl;
        c0 = 2.9979e8;
        omega = k0 * c0;
        mu0 = 4*pi * 1e-7;
        eps0 = 1/c0^2/mu0;
        Gx(m,:) = 1j*k0*neff * Fy(m,:);
        Gz( m , iL ) = AL * exp( +gL*( x(iL)-hs(1)   ) ); % Left substrate
        Gz( m , iR ) = AR * exp( -gR*( x(iR)-hs(end) ) ); % Right substrate
        for jj = 1:length(ns)
            ix = x>=hs(jj) & x<hs(jj+1);
            Gz( m , ix ) = A(jj) * -sin( kn(jj)*x(ix) - the(jj) ) * kn(jj);
        end
        if strcmp( ModePol , 'TM' )
            Gx(m,:) = Gx(m,:) ./ er / ( +1j*omega*eps0 );
            Gz(m,:) = Gz(m,:) ./ er / ( +1j*omega*eps0 );
        else
            Gx(m,:) = Gx(m,:) ./ 1 / ( -1j*omega*mu0 );
            Gz(m,:) = Gz(m,:) ./ 1 / ( -1j*omega*mu0 );
        end
        
        % Power-flux along the z-axis
        if strcmp( ModePol , 'TE' )
            Pz(m,:) = 0.5*real( -Fy(m,:) .* conj( Gx(m,:) ) );
        elseif strcmp( ModePol , 'TM' )
            Pz(m,:) = 0.5*real( +Gx(m,:) .* conj( Fy(m,:) ) );
        end
        
        % normalize for max-abs == 1
        NormCost = max(abs(Fy(m,:)));
        Fy(m,:) = Fy(m,:) / NormCost ;
        Gx(m,:) = Gx(m,:) / NormCost ;
        Gz(m,:) = Gz(m,:) / NormCost ;
        Pz(m,:) = Pz(m,:) / max(Pz(m,:)) ;
        
    end
    
end

% -------------------------------------------------------------------------
% Plot: Mode Profiles
% -------------------------------------------------------------------------
if DoPlotModeProfiles == 1
    figure;
    nr = ceil( sqrt( length(neffs) / (1920/1080) ) ); % num-rows in subplot
    nc = ceil( length(neffs) / nr ); % num-cols in subplot
    for m = 1:length(neffs)
        subplot(nr,nc,m)
        plot( x , real(Fy(m,:)) , 'b' ); hold on;
        plot( x , imag(Fy(m,:)) , 'r' );
        plot( x , Pz(m,:) , 'ko' , 'markersize' , 1.5 );
        legend( 'Re(Fy)' , 'Im(Fy)' , 'Pz' );
        title( sprintf( '\\bf[%s.%d] \\rm neff = %6.4f' , ModePol , m , neffs(m) ) )
        for jj = 1:length(hs)
            plot( hs(jj)*[1 1] , [-1 1] , 'k-' )
        end
        plot( x ,0*x , 'k:' );
        axis( [x([1 end]) -1 1] )
        
    end
    
    % if ~isempty( Pz )
    %     figure;
    %     plot( x , Pz , 'ko' , 'markersize' , 2 ); hold on;
    %     for jj = 1:length(hs)
    %         plot( hs(jj)*[1 1] , [-1 1] , 'k-' )
    %     end
    %     plot( x ,0*x , 'k:' );
    %     axis( [x([1 end]) -1 1] )
    % end
    
    FixMultiFigsPos
    
end

end




