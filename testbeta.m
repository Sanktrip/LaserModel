close all
clear

qe=1.60217662e-19;

E0 = linspace(0.4,2.0,1000);
Nc = 1.0e20;
[gain,beta,Rspon] = gain_numerical_cgp(E0,Nc);

%[gain2,E2,beta2] = gain_numerical_parya(Nc);

%disp(beta)
%disp(beta2)

% 
% close all
% figure
% plot(E2./qe,gain2)
% hold on
% plot(E0,gain)


Ncrange = logspace(20,25,20);
ib = 0;
for Nc = Ncrange
    ib = ib+1;
    [gain,beta,Rspon] = gain_numerical_cgp(E0,Nc);
    %[gain2,E2,beta2] = gain_numerical_parya(Nc);

    beta_vals(ib) = beta;
    %beta2_vals(ib) = beta2;

    gain_vals(ib,:) = gain;
    %gain_vals2(ib,:) = gain2;
    rspon_vals(ib,:) = Rspon;
end

close all
%plot(E0.*qe,gain_vals)
%hold on
%plot(E2,gain_vals2,'k')
%figure
%semilogx(Ncrange,beta2_vals)
%hold on
%plot(Ncrange,beta_vals)
plot(E0.*qe,rspon_vals);
figure;
hold on


