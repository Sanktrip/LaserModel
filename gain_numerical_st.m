function [gain, sponSpec] = gain_numerical_st(E0,Nc)

% Bulk Gain and Spontaneous Emission

global Eg kBT hbar eps0 na m0 c Lz mr D2 Gamma_in qe maxE NPTS me mh h tin


Erange0 = reshape(E0,1,length(E0)).*qe; 
Ecvrange = linspace(Eg,maxE,NPTS+1).';


%[Fc, Fv] = getQuasiFermiLevelsBulk(Nc);  
[Fc, Fv] = getQuasiFermiLevels(Nc, Lz);

fc = @(Ecv) 1./(1+exp((Eg+mr/me*(Ecv-Eg)-Fc)./kBT)); 
fv = @(Ecv) 1./(1+exp((-mr/mh*(Ecv-Eg)-Fv)./kBT));

L = @(E,Ecv, gamma) gamma/(2*pi)./((Ecv-E).^2 + (gamma/2).^2);
C0 = @(Ecv) hbar*pi*qe^2./(eps0 * na* c* m0^2.*Ecv);

pd = @(E) 1/(2*pi^2) * (2*mr/hbar^2)^(3/2) * sqrt(E - Eg); % Bulk DOS

gainIntegral = @(Ecv) C0(Erange0).*pd(Ecv).*(fc(Ecv) - fv(Ecv)).*D2.*L(Erange0, Ecv ,Gamma_in);
gain = trapz(Ecvrange, gainIntegral(Ecvrange), 1);


% Spontaneous Emission Spectrum
freespaceDOS = @(E) (8*pi*na^3.*E.^2)./(h^3*c^3);

sponIntegral = @(Ecv) pd(Ecv).*fc(Ecv).*(1-fv(Ecv)).*D2.*L(Ecv, Erange0, Gamma_in);
sponSpec = freespaceDOS(Erange0) .* C0(Erange0).*  c/na .* trapz(Ecvrange, sponIntegral(Ecvrange),1);










