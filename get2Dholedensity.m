function Nh = get2Dholedensity(Fv,Lz)
% Computes the excited electron number density for a 2D quantum well, given the Quasi-fermi
% level F which is given in Joules measured from the top of the valence
% band.

global hbar mh qe NBANDS kB T

nb = 1:NBANDS;

%compute the energy bands
En = -(nb.*pi./Lz).^2.*hbar^2./(2*mh);

nsum = 0;
for n = 1:NBANDS
    if exp((En(n)-Fv)./(kB*T))<10000
        nsum = nsum + log(1+exp((En(n)-Fv)./(kB*T)));
    else
        nsum = nsum + (En(n)-Fv)./(kB*T);
    end
    
end
Nh = nsum.*(mh.*kB.*T)./(pi*hbar^2*Lz);

% hold on
% plot(Fv,Nh,'ko')
% pause

end

