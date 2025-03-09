function Ne = get3Delectrondensity(Fc)
% Computes the excited electron number density for a 2D quantum well, given the Quasi-fermi
% level F which is given in Joules measured from the top of the valence
% band.

global hbar me Eg kBT mr
integrand = @(E) 1./(1+exp( (Eg+mr/me*(E-Eg) - Fc)./kBT)) .* sqrt(E - Eg);

nsum = integral(@(E) integrand(E), Eg, inf);

Ne = nsum .* 1/(2*pi^2) .* (2*me/hbar^2)^(3/2);

%hold on
%plot(F,Ne,'ko')
%pause

end

