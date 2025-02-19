close all
clear

qe=1.60217662e-19;

E0 = linspace(-0.05, 0.25,1000);   % In eV
NcRange = [4.7e18, 3.6e18, 2.6e18, 1.7e18];

ib = 0;
for Nc = NcRange
    ib = ib+1;
    [gain] = gain_numerical_st(E0, Nc);
    gain_vals(ib,:) = gain;
end

close all
plot(E0, gain_vals);
figure 