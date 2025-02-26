function [Fc, Fv] = getQuasiFermiLevelsBulk(Nc)


global Eg qe

Nfunc = @(Fc) (get3Delectrondensity(Fc)-Nc);

OPTS = optimset('TolX',qe*1e-6);
Fc = fzero(Nfunc,Eg,OPTS);

Nfunc2 = @(Fv) (get3Dholedensity(Fv)-Nc);

OPTS = optimset('TolX',qe*1e-8);
Fv = fzero(Nfunc2,-1e-20,OPTS);
