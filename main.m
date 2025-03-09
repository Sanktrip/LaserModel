
close all
clear
global qe c hbar kB kBT maxE Eg NPTS tin Gamma_in na m0 eps0 h Lz T Ep D2 
global me mh mr h 

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
D2 = m0*Ep/3;           % square of the dipole interaction strength

me = m0 * 0.062;   % Effective electron mass
mh = m0 * 0.5;     % Effective hole mass
mr = m0 * 0.055;   % reduced effective mass
tin = 40e-15;           % broadening lifetime /s
Gamma_in = 2*hbar/tin;  % Lorentzian linewidth /J
na = 3.7;               % refractive index of the gain material

Eg=1.306*qe;           % band gap energy /J

NPTS = 5000;            % number of points of integration over Ecv
maxE = 5*Eg;            % maximum energy to integrate to


% Parameters 
%E0 = linspace(-0.05, 0.25,100);   % In eV
E0 = linspace(0.4,2.0,1000);
NcRange = logspace(20,25,20);
%NcRange = [4.7e18, 3.6e18, 2.6e18, 1.7e18] .* 1e6; % Carrier densities

% Test Bulk

ib = 0;
for Nc = NcRange
    ib = ib+1;
    [gain, sponSpec] = gain_numerical_st(E0, Nc);
    gain_vals(ib,:) = gain;
    sponSpec_vals(ib, :) = sponSpec;
 
 end
 
plot(E0, gain_vals./1e2);

hold on
%plot(E0.*qe, sponSpec_vals);
% figure 

% Test 2D
%{
ib = 0;
for Nc = NcRange
    ib = ib+1;
    [gain2D, sponSpec2D] = gain2D_numerical_st(E0, Nc);
    gain2D_vals(ib,:) = gain2D;
    sponSpec2D_vals(ib, :) = sponSpec2D;

end

close all
%plot(E0.*qe, gain2D_vals);
%hold on
lambda = (h.*c)./(E0.*qe) * 1e9;
plot(lambda, sponSpec2D_vals);
hold on
%figure 
%}



