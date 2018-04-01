function [Npre,slope1,slope2] = Nt(I0,tp,t,lambda1,Ip)
%calculate the pre-exponential factor of the Yudin ionization rate
%Npre = pre-exponential factor
%slope1 and slope2 = asymptote

%cst
etha0	= 377;
me		= 9.109e-31;
hbar	= 1.055e-34; 
c0		= 299792458;
w1		= c0*2*pi/lambda1;
e		= 1.602e-19;



sigma	= tp/(2*sqrt(2*log(2)));
I		= I0*exp(-t.^2/(2*sigma^2));



%prexeponential factor calculation
Up      = 0.0933*lambda1^2.*I*1e-4; %eV
gamma0   = sqrt(Ip./(2.*Up));

nstar   = 1/sqrt(2*(Ip/27.211)); %effective principal quantum number
lstar   = 0; %effective angular momentum
l       = 0; %actual angular momentum
m       = 0; %projection of l on the laser polarisation vector


Anl = 2^(2*nstar)/(nstar*gamma(nstar+lstar+1)*gamma(nstar-lstar));
Blm = (2*l+1)*factorial(l+m)/(2^m*factorial(m)*factorial(l-m));
k = log(gamma0 + sqrt(gamma0.^2 + 1)) - gamma0./sqrt(gamma0.^2 + 1);
C = 1*(gamma0 < 1.2) + (1.2)./sqrt(gamma0).*(gamma0 > 1.2);

Npre =Anl*Blm*sqrt(3*k./(gamma0.^3)).*C.*(Ip*e/hbar).*(4*Ip*e/hbar*sqrt(2*me*Ip*e)./(e*sqrt(2*etha0*I))).^(2*nstar-m-1);
slope1 = 10^(-152.35)*I.^(Ip/(hbar*w1)/6.242e18+(-1/2*(2*nstar-m-3/2))+0.6);
slope2 = 10^(27.2)*I.^(-1/2*(2*nstar-m-1));
end

