function [Fc, Fv] = getQuasiFermiLevels(Nc,Lz)
% function to compute the Quasi-fermi levels for a given carrier density
% Nc.

global Eg me mh qe


Nfunc = @(Fc) (get2Delectrondensity(Fc,Lz)-Nc);

OPTS = optimset('TolX',qe*1e-6);
Fc = fzero(Nfunc,Eg,OPTS);

Nfunc2 = @(Fv) (get2Dholedensity(Fv,Lz)-Nc);

OPTS = optimset('TolX',qe*1e-8);
Fv = fzero(Nfunc2,-1e-20,OPTS);


return


