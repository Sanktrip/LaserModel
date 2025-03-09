function Ne = get3Dholedensity(Fv)
% Computes the excited electron number density for a 2D quantum well, given the Quasi-fermi
% level F which is given in Joules measured from the top of the valence
% band.

global hbar mh Eg kBT mr

integrand = @(E) 1./(1+exp( (-mr/mh*(E-Eg)-Fv)./kBT)) .* sqrt(E - Eg);

nsum = integral(@(E) integrand(E), Eg, inf);

Ne = nsum .* 1/(2*pi^2) .* (2*mh/hbar^2)^(3/2);

% hold on
% plot(F,Nd,'ko')
% pause

end

