function [gain,beta,Rspon_total] = gain_numerical_cgp(E0,Nc)
% computes the gain and beta factor for a given range of E0
% Nc is the carrier density in SI units.
% E0 should be given in eV (not SI, but this is ok for the moment)
% 
% Based on the notes of Chris/Arti, early 2023
%

global qe hbar cspeed eps0 kB   % fundamental constants
global T na nc Eg me mh mr      % material constants
global NBANDS

DEBUG = 0;

% fundamental constants. Everything is in SI units.

qe=1.60217662e-19;           % charge on an electron /C
h = 6.62607015e-34;          % planck's constant /Js
hbar = h/(2*pi);             % reduced planck's constant
c = 2.99792458e8;            % speed of light /(m/s)
m0 = 9.109e-31;              % electron rest mass / kg
eps0 = 8.854e-12;            % permittivity of free space /SI 
kB = 1.3807e-23;             % Boltzmann constant

% material properties

T = 5;                  % temperature /K
na = 3.6;               % refractive index of the gain material
nc = 3.6;               % refractive index of cavity

Eg=1.2380*qe;           % band gap energy /J
Ev = 0;                 % top of valence band
Ec = Eg;                % botton of conduction band

tin = 40e-15;           % broadening lifetime /s
Gamma_in = 2*hbar/tin;  % Lorentzian linewidth /J

Ep=25.18*qe;
D2 = m0*Ep/3;           % square of the dipole interaction strength

PInAs = 0.15;           % proportion of InAs in the QW 
PGaAs = 1-PInAs;        % proportion of GaAs in the QW

me=(PInAs*0.023 + PGaAs*0.067)*m0;           % effective mass of electrons in conduction band (interpolation using 85% GaAs)
mh100=(PInAs*0.35 + PGaAs*0.34)*m0;             % effective mass of heavy holes in valence band (geometric mean, from Mollenkamp measurements, via Nawaski)
mh111=(PInAs*0.43 + PGaAs*0.7)*m0;              % InAs is room termperature from [169] in Nawaski
mh= sqrt(mh100*mh111);      %take the geometric mean

mr=(me*mh)/(me+mh);     % reduced effective mass

NBANDS = 5;             % number of energy bands considered in the Quantum well

%Resonator and QW geometry

d = 200e-9;             % cavity width /m
L=2.2e-6;               % cavity length /m
Vcav = (sqrt(3)/2)*d^2*L; % cavity volume /m^3

Lz = 19e-9;             % quantum well width /m

% computation parameters

NPTS = 5000;            % number of points of integration over Ecv
maxE = 5*Eg;            % maximum energy to integrate to
Erange0 = reshape(E0,1,length(E0)).*qe; 
Erange1 = Erange0;
Ecvrange = linspace(Eg,maxE,NPTS+1).';

% preliminaries:

kBT = kB.*T;
pr = mr/(pi*hbar^2*Lz)                 % compute electronic joint DOS (not used)

% optical cavity resonances

Fcav = [281.7788 297.3256 298.3676... %List of cavity frequencies / Hz
        300.5511 308.0342 312.2962...
        320.6671 326.2252 333.7105...
        339.8508].*1e12;      
Qcav = [9.7742 20.7069 5.2981 ... %List of quality factors
        5.1025 10.0198 33.3025 ...
        19.0478 44.3064 32.7116 ...
        54.3495];                         

NMODE = 6; % pick out the lasing mode - this will be the mode for computing beta.

Gcav = Fcav./Qcav;                                                                                             % frequency linewidth of cavities (1/s)
GammaCav = h.*Gcav;                                                                                            %energy linewidth of the cavity (J)
Ecav = h*Fcav;                                                                                                 %energies of the cavity (J)

% compute Quasi-Fermi-levels:

[Fc, Fv] = getQuasiFermiLevels(Nc,Lz);

% define functions that have to be integrated over

C0 = @(E) hbar*pi*qe^2./(na*c*eps0*m0^2.*E); %C0 function from notes
L = @(E,Ecv,Gamma) Gamma/2./((Ecv-E).^2+(Gamma/2).^2)/pi;     % Lorentzian
p2D = @(Ecv) rho2Dr(mr,Ecv,Lz); %2D Density of States

fc = @(Ecv) 1./(1+exp((Eg+mr/me*(Ecv-Eg)-Fc)./kBT)); % carrier distribution
fv = @(Ecv) 1./(1+exp((  -mr/mh*(Ecv-Eg)-Fv)./kBT)); % vacancy distribution

if DEBUG
    figure
    plot(Ecvrange,p2D(Ecvrange))
end

% compute the gain:
E = Erange0; 
integrand_over_Ecv = @(Ecv) p2D(Ecv).*(fc(Ecv)-fv(Ecv)).*D2.*L(E,Ecv,Gamma_in);

gain = C0(E).*trapz(Ecvrange,integrand_over_Ecv(Ecvrange),1);

% compute photon DOS fro the cavity modes
E = Erange1;
Nph = zeros(size(E));
for j = 1:length(Fcav)
    Nph = Nph + (1/Vcav)*L(E,Ecav(j),GammaCav(j)); 
end
Nph0=(8.*pi.*(nc.^3).*(E.^2))./((h.^3).*(c.^3));                             
Nph = Nph+Nph0;

Nphj=(1/Vcav)*L(E,Ecav(NMODE),GammaCav(NMODE));

if DEBUG
    figure
    plot(E,Nph)
end

% compute spontaneous emission spectrum:

integrand_over_Ecv2 = @(Ecv) p2D(Ecv).*fc(Ecv).*(1-fv(Ecv)).*D2.*L(E,Ecv,Gamma_in);

rs_total_E = Nph.*C0(E).*(c/nc).*trapz(Ecvrange,integrand_over_Ecv2(Ecvrange),1);

if DEBUG
    figure
    plot(E,rs_total_E)
end

% compute spontaneous emission rate into the j-th mode:
rs_j_E = Nphj.*C0(E).*(c/nc).*trapz(Ecvrange,integrand_over_Ecv2(Ecvrange),1);


if DEBUG
    hold on
    plot(E,rs_j_E)
    pause
end

Rspon_total = trapz(E,rs_total_E);
Rspon_j = trapz(E,rs_j_E);

beta = Rspon_j/Rspon_total;


