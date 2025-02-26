function p2D = rho2Dr(m,E,Lz)
% computes the 2D DOS. m is the reduced mass. E is the Energy in J

global hbar Eg qe NBANDS

N = NBANDS;      %number of energy bands in the quantum well
nb = 1:N;

En = 2*(nb.*pi./Lz).^2.*hbar^2./(2*m);  % note that the bands are spaced twice as far 
                                        % apart in the reduced coordinate system

Enqe = En./qe;

p0 = m/(pi*hbar^2*Lz);

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
