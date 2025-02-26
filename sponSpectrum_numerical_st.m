function [sponSpec] = sponSpectrum_numerical_st(E0,Nc)

qe=1.60217662e-19;
tin = 1e-14;
Eg=1.306*qe;           % band gap energy /J

% computation parameters

NPTS = 5000;            % number of points of integration over Ecv
maxE = 5*Eg;            % maximum energy to integrate to

c = 2.99792458e8;            % speed of light /(m/s)
m0 = 9.109e-31;              % electron rest mass / kg
eps0 = 8.854e-12;            % permittivity of free space /SI 
kB = 1.3807e-23;             % Boltzmann constant
h = 6.62607015e-34;          % planck's constant /Js
hbar = h/(2*pi);             % reduced planck's constant


T = 4;                          % temperature /K
Lz = 19e-9;                     % quantum well width /m
kBT = kB.*T;

Ep=25.18*qe;
D2 = m0*Ep/6;           % square of the dipole interaction strength

me = m0 * 0.062;   % Effective electron mass
mh = m0 * 0.5;     % Effective hole mass
mr = m0 * 0.055;   % reduced effective mass
na = 3.7;               % refractive index of the gain material


Erange0 = reshape(E0,1, length(E0)) .* qe;
%Erange0 = [1 3 5; 2 4 6; 7 8 10];
%Erange0 = E0;
%Erange0 = Erange0.*qe;
%Erange0 = E0*qe;
Ecvrange = linspace(Eg,maxE,NPTS+1).';

[Fc, Fv] = getQuasiFermiLevels(Nc,Lz);  

fc = @(Ecv) 1./(1+exp((Eg+mr/me*(Ecv-Eg)-Fc)./kBT)); 
fv = @(Ecv) 1./(1+exp((-mr/mh*(Ecv-Eg)-Fv)./kBT));
C0 = @(Ecv) hbar*pi*qe^2./(eps0 * na^2* m0^2.*Ecv);
freespaceDOS = @(E) (8*pi*na^3.*E.^2)./(h^3*c^3);
p2d = @(Ecv) rho2Dr(mr,Ecv,Lz); 
L = @(Ecv, E, t) hbar/(pi*t)./((Ecv - E).^2 + (hbar/t).^2);

sponIntegral = @(Ecv) p2d(Ecv).*fc(Ecv).*(1-fv(Ecv)).*D2.*L(Ecv, Erange0, tin);
%Ecv1 = @(n1, n2) pi^2*h^2/(2*mr*Lz)*(1);
%e_gap = Ecv1(1,2);
%sponSpec = freespaceDOS(E0) .* C0(E0).*trapz(e_gap, sponIntegral(e_gap),1);
sponSpec = freespaceDOS(Erange0) .* C0(Erange0).*trapz(Ecvrange, sponIntegral(Ecvrange),1);




