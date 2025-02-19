function p2D = rho2D(m,E,Lz)
% computes the 2D DOS. m is the relevant mass. E is the Energy in J

global hbar Eg qe NBANDS

N = NBANDS;      %number of energy bands in the quantum well
nb = 1:N;

En = (nb.*pi./Lz).^2.*hbar^2./(2*m)

Enqe = En./qe;

p0 = m/(pi*hbar^2*Lz)

psum = 0;
for n = nb
    psum = psum + real((E-Eg-En(n))>0);
end

p2D = psum.*p0;

end

% function HV(x)
% %heaviside function
%     HV(x) = real(x>0);
% end
