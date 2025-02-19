function Nd = get2Dchargedensity(F,Eoffset,m,Lz)
% Computes the charge density for a 2D quantum well, given the Quasi-fermi
% level F which is given in Joules measured from the top of the valence
% band. Eoffset is the offset energy from which the bands are measured: this is
% Eg for the conduction band and 0 for the valence band.

global hbar qe NBANDS kB T

nb = 1:NBANDS;

%compute the energy bands
En = 2*(nb.*pi./Lz).^2.*hbar^2./(2*m);

nsum = 0;
for n = 1:NBANDS
    nsum = nsum + log(1+exp((F-Eoffset+En(n))./(kB*T)));
end
Nd = nsum.*(m.*kB.*T)./(pi*hbar^2*Lz);

% hold on
% plot(F,Nd,'ko')
% pause

end

