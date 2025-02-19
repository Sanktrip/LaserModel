close all
clear
global c
qe=1.60217662e-19;

E0 = linspace(0.4,2.0,1000);
Nc = 1.0e20;

Ncrange = logspace(20,25,20);
ib = 0;
for Nc = Ncrange
    ib = ib+1;
    [sponSpec] = sponSpectrum_numerical_st(E0, Nc);
    [gain,beta,Rspon_total, em] = gain_numerical_cgp(E0, Nc);
    
    expected_vals(ib, :) = em;
    sponSpec_vals(ib, :) = sponSpec;
end

close all
plot(E0, sponSpec_vals);
figure
hold on 
plot(E0,expected_vals);
