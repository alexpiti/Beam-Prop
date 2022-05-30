function [ XE, dXE ] = MLSWG_CharEq( ne, nL, nR, ns, wl, ts, ModePol, ForceReal )

% MultiLayer Slab WaveGuide (MLSWG) - Characteristic Equation
%
% This function is meant for the numerical solution of the transcedental
% characteristic equation (CharEq or XE) of a multilayer 1D/slab dielectric
% waveguide. It returns the XE and dXE/dneff value for a given n_eff value.
% The dXE/dneff is required for the Newton-Raphson (NR) method.
%
% - ne     : n_effective (scalar)
% - nL/nR  : refr.indx of Left- and Right-most semi-inf layers (scalars)
% - ns     : refr.indx of intermediate layers (1D vector)
% - wl     : wavelength (real, same length-units as in "ts")
% - ts     : thicknesses of intermediate layers "ns" (1D vector)
% - ModePol   : Mode-Polarization, 'TE' or 'TM'
% - ForceReal : If ==1, then remove Imag-part of XE (default=0)
%
% Refractive indices can be complex-valued, but that typically makes the
% numerical solutions of the XE of multi-layer/desync'd WGs even more
% difficult. Please note that only the NR and Secant methods can find roots
% of complex-valued nonlinear equations.
%
% Alexandros Pitilakis / Thessaloniki, Greece, 
%  2015 Sept-Oct: Original version

% Test inputs
if nargin == 0
   nL = 1.45;
   nR = 1.00;
   ns = [3.20 2.00 3.20-0.000j];
   wl = 1.55;
   ts = [0.40 1.00 0.40];
   ModePol = 'TE';
   ne = linspace( max( [nR,nL] ) , max(ns) , 1e4 );
   ForceReal = 0;
end

if nargin < 8, ForceReal = 0; end % Force CharEq to Real, removing Imag

% ------------------------------------------------------------------------
% Some Aux Variables
% ------------------------------------------------------------------------

% Free-space wavenumber 
k0 = 2*pi / wl; % [um^-1]

% Number of intermediate Layers (minimum:1)
N = length(ns);

% L/R layer variables (same size as ne space). Field in those layers is of
% the form: A*exp(-g*(x-d)), to exponentially decline inside them.
gL = k0 * sqrt( ne(:).^2 - nL^2 ); % always Real/positive => ne>nL
gR = k0 * sqrt( ne(:).^2 - nR^2 ); % always Real/positive => ne>nR

% pdFL/RL: polarization-dependent Factor used at the interfaces between
% [Left|#1] and [#N|Right] layers
if strcmp( ModePol , 'TM' )
    pdFL = ns( 1 )^2 / nL^2; 
    pdFR = ns(end)^2 / nR^2; 
else
    pdFL = 1;
    pdFR = 1;
end

% Intermediate layer variables (each has same size as ne space)
kn = NaN*zeros( length(ne), length(ns) );
pdFn = ones(size(ns)); % polarization-dependent Factor
for ii = 1 : length(ns)
    
    kn(:,ii) = k0 * sqrt( ns(ii)^2 - ne.^2 ); % Real or Imag (not complex!)
    
    % polarization-dependent Factor  for interfaces between intermediate layers 
    if ii>=2 && strcmp( ModePol , 'TM' )
        pdFn(ii) = ns(ii)^2 / ns(ii-1)^2; 
    end   
    
end

% Assuming the interface [L|#1] is at x=0, we calc the locations of the 
% rest of the interfaces from the layer-thicknesses
hn = cumsum( ts );

% ------------------------------------------------------------------------
% CharEq (XE) --> Used for Interpolation method
% ------------------------------------------------------------------------

% In order to form the CharEq, we scan the boundary conditions (BC) at all
% the layer-interfaces. We move from the [L|#1] interface, to [#1|#2], and
% finally to [#N|R]. The parameter involved is here called "Theta".
Th = NaN*zeros( length(ne), length(ns) );
Th(:,1) = atan( gL ./ kn(:,1) * pdFL ); % from the [L|#1] interface BC
for jj = 2 : N   
    auxi = tan( kn(:,jj-1)*hn(jj-1) - Th(:,jj-1) );
    Th(:,jj) = kn(:,jj)*hn(jj-1) - atan( pdFn(jj) * kn(:,jj-1)./kn(:,jj) .* auxi );
end

% Having calculated the "Theta" angle of the last interface [#N|R], we form
% the final CharEq (XE):
XE = kn(:,N) .* tan( kn(:,N)*hn(N) - Th(:,N)  ) - gR * pdFR;

if ForceReal == 1
    XE = real(XE);
    if imag(mean(XE))/real(mean(XE)) > 1e-12
       disp( ' ## MLSPWG_CharEq: Warning -- Imag(XE)/Real(XE) is not negligible!' ); 
    end
end

% Note: Re-casting of the CharEq as
%        XE = tan( kn(:,N)*hn(N) - Th(:,N) ) - gR * pdFR ./ kn(:,end) ;
% leads to a more "steep" curve, so that its roots are harded-to-find. 
% So, we stick to the previous casting: XE=k*tan(...)-gR

% ------------------------------------------------------------------------
% Derivatives of CharEq (XE) --> Used for a Newton-Raphson Method
% ------------------------------------------------------------------------

% Derivatives of k,gL/gR with respect to neff (ne)
dkn = NaN * zeros( length(ne) , 1 );
for jj = 1 : N
    dkn(:,jj) = -k0^2*ne(:)./kn(:,jj);
end
dgL = +k0^2*ne(:)./gL; % TE
dgR = +k0^2*ne(:)./gR; % TE


% Derivative of "Theta_N" wrt neff, calculated iteratively. We start from
% Theta_1 [L|#1] interface, and move onward until [#N|R] interface.
dTh = NaN*zeros( length(ne), length(ns) );
dTh(:,1) = pdFL ./( 1 + (gL./kn(:,1)*pdFL).^2 ).* (  dgL.*kn(:,1) - dkn(:,1).* gL ) ./ kn(:,1).^2 ;
for jj = 2 : N
   
    Phi = kn(:,jj-1)*hn(jj-1) - Th(:,jj-1); % aux var, here defined at jj-1
    Om = pdFn(jj) * kn(:,jj-1)./kn(:,jj) .* tan( Phi ); % another aux var (check my notes)
    
    % Parts and full expression of derivative of this aux Om(ega) var
    dOm1 = ( dkn(:,jj-1).*kn(:,jj) - dkn(:,jj).* kn(:,jj-1) ) ./ kn(:,jj).^2 .* tan(Phi) ;
    dOm2 = kn(:,jj-1) ./ kn(:,jj) .* ( 1 + tan(Phi).^2 );
    dOm3 = dkn(:,jj-1)*hn(jj-1) - dTh(:,jj-1);
    dOm = pdFn(jj) * ( dOm1 + dOm2 .* dOm3 );
    
    dTh(:,jj) = dkn(:,jj)*hn(jj-1) - dOm ./ ( 1 + Om.^2 );
    
end

% Final expression of the dXE/dneff
Phi = kn(:,N)*hn(N) - Th(:,N); % aux var, defined at jj==N
if strcmp( ModePol , 'TM' )
    dXE = ns(end)^-2 * ( dkn(:,N).*tan(Phi) + kn(:,N).*(1+tan(Phi).^2).*( dkn(:,N)*hn(N) - dTh(:,N) ) ) ...
        - nR^-2 * dgR;
else
    dXE = dkn(:,N).*tan(Phi) + kn(:,N).*(1+tan(Phi).^2).*( dkn(:,N)*hn(N) - dTh(:,N) ) ...
        - dgR;
end
%dXE = real(dXE);

% ------------------------------------------------------------------------
% Test plots
% ------------------------------------------------------------------------
if nargin == 0
    
    figure;
    
    Xfoc = ne([1,end]);
    
    subplot(2,1,1)
    plot( ne , real(XE) , 'bo' ); hold on;
    plot( ne , imag(XE) , 'co' , 'markersize' , 2 );
    set(gca,'YLim',1e2*[-1 +1]);
    set(gca,'XLim',Xfoc);
    legend( 'Real(XE)' , 'Imag(XE)' );

    subplot(2,1,2)
    plot( ne , real(dXE) , 'ro' ); hold on;
    plot( ne , imag(dXE) , 'yo' , 'markersize' , 2 );
    set(gca,'YLim',1e3*[-1 0]);
    set(gca,'XLim',Xfoc);
    legend( 'Real(dXE)' , 'Imag(dXE)' );
    
end

