function [gain] = gain_numerical_st(E0,Nc)

qe=1.60217662e-19;
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
tin = 40e-15;           % broadening lifetime /s
Gamma_in = 2*hbar/tin;  % Lorentzian linewidth /J
na = 3.7;               % refractive index of the gain material

Eg=1.306*qe;           % band gap energy /J

NPTS = 5000;            % number of points of integration over Ecv
maxE = 5*Eg;            % maximum energy to integrate to
Erange0 = reshape(E0,1,length(E0)).*qe; 
Ecvrange = linspace(Eg,maxE,NPTS+1).';


[Fc, Fv] = getQuasiFermiLevelsBulk(Nc);  

fc = @(Ecv) 1./(1+exp((Eg+mr/me*(Ecv-Eg)-Fc)./kBT)); 
fv = @(Ecv) 1./(1+exp((-mr/mh*(Ecv-Eg)-Fv)./kBT));

L = @(E,Ecv, gamma) gamma/(2*pi)./((Ecv-E).^2 + (gamma/2).^2);
C0 = @(Ecv) hbar*pi*qe^2./(eps0 * na* c* m0^2.*Ecv);
%p2d = mr/(pi * hbar^2 * Lz);           
%p2d = @(Ecv) rho2Dr(mr,Ecv,Lz);    % Quantum Well DOS
pd = @(E) 1/(2*pi^2) * (2*mr/hbar^2)^(3/2) * sqrt(E - Eg); % Bulk DOS

gainIntegral = @(Ecv) C0(Erange0).*pd(Erange0).*(fc(Ecv) - fv(Ecv)).*D2.*L(Erange0, Ecv ,Gamma_in);
%gainIntegral = @(Ecv) C0(Erange0).*p2d.*(fc(Ecv) - fv(Ecv)).*D2.*L(Erange0, Ecv ,Gamma_in);

gain = trapz(Ecvrange, gainIntegral(Ecvrange), 1);




