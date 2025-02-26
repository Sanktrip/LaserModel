function Ne = get2Delectrondensity(Fc,Lz)
% Computes the excited electron number density for a 2D quantum well, given the Quasi-fermi
% level F which is given in Joules measured from the top of the valence
% band.

global hbar me Eg qe NBANDS kB T

nb = 1:NBANDS;

%compute the energy bands
En = (nb.*pi./Lz).^2.*hbar^2./(2*me);

nsum = 0;
for n = 1:NBANDS
    if exp((Fc-Eg-En(n))./(kB*T))<10000
        nsum = nsum + log(1+exp((Fc-Eg-En(n))./(kB*T)));
    else
        nsum = nsum + (Fc-Eg-En(n))./(kB*T);
    end
    
end
Ne = nsum.*(me.*kB.*T)./(pi*hbar^2*Lz);

% hold on
%plot(F,Nd,'ko')
% pause

end

