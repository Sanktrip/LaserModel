function [sponSpec] = sponSpectrum_numerical_st(E0,Nc)

global hbar h qe kBT Lz eps0 na m0 mr Eg me D2 c mh

tin = 1e-14;

% computation parameters

NPTS = 5000;            % number of points of integration over Ecv
maxE = 5*Eg;            % maximum energy to integrate to
Erange0 = reshape(E0,1,length(E0)).*qe; 
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




