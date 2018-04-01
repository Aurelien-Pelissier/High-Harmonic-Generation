function [Phi,K] = phase(z,r,t,I0,tp,Pm,lp,profile,znozzle,q,alpha,lambda1,R0,gas,zmax,rmax,nres,ethaft,xHe)
%main function to calculate phase matching DeltaPhi and DeltaK%
%Phi = phase mismatch array (r,z)
%K = wavevector mismatch array (r,z)

%universe constants
c		= 299792458;	
hbar	= 1.055e-34; 	
e		= 1.602e-19;		
me		= 9.109e-31;	
eps0	= 8.854e-12; 

%experiment constants
w1		= c*2*pi/lambda1;       %fundamental pulsation
k1		= w1/c;                 %fundamental wave vector	
wq		= q*w1;                 %harmonic pulsation							
kq		= wq/c;                 %harmonic wave vector
Ne		= Ndens(gas);           %electron density at STP for eta = 1
wp		= sqrt(Ne*e*e/me/eps0); %plasma frequency
P0		= 1013;                 %pressure at STP = 1atm
deltan	= refractive(lambda1,gas) - refractive(lambda1/q,gas);  %difference between the fundamental and high order harmonic refractive index%								
deltaHe =  refractive(lambda1,'He') - refractive(lambda1/q,'He'); %only important for non zero helium fraction
Zr		= k1*R0^2/2;            %Rayleigh length
b       = Zr*2;                 %confocal parameter



%gaussian beam
[I,dIz,dIr] = Igauss(I0,tp,R0,lambda1,z,r,t);


%pressure
[P,dPz,dPr] = Press(z,r,Pm,lp,profile,znozzle);


%ionization fraction
[etha,dethaz,dethar] = ioniz(z,r,zmax,rmax,nres,ethaft);


%refractive index
n1		= 1+P/P0*((1-etha)*deltan-wp^2/(2*w1^2)*etha) + P/P0*xHe/(1-xHe)*deltaHe; %last part is helium index modification
nq		= 1-P/P0*(wp^2/(2*wq^2)*etha);
dn1z	= (wp^2/(2*P0*w1^2))*(dPz*etha+P*dethaz) + dPz/P0*(1-etha)*deltan - P/P0*dethaz*deltan;
dn1r	= (wp^2/(2*P0*w1^2))*(dPr*etha+P*dethar) + dPr/P0*(1-etha)*deltan - P/P0*dethar*deltan;
dnqz	= (wp^2/(2*P0*wq^2))*(dPz*etha+P*dethaz);
dnqr	= (wp^2/(2*P0*wq^2))*(dPr*etha+P*dethar);




%-------------------------- phase --------------------------%

%dipole phase
Phid = -alpha*I;


%focusing phase
Phi1 = k1*n1*z + angle((1/(n1*b+2i*z))*exp(-k1*n1*r^2/(n1*b+2i*z)));
Phiq = kq*nq*z; %+ angle((1/(nq*bq+2i*z))*exp(-kq*nq*r^2/(nq*bq+2i*z))); %commented because negligible anyway


%Phi modulo 2Pi
Phit = Phiq - (q*Phi1 + Phid);
Phi = mod(Phit,2*pi);

%--------------------------- K ---------------------------%
%K=grad(Phi)

%dipole wavevector
Kzd = -alpha*dIz;
Krd =  -alpha*dIr;

%plane wave wavevector
Kzq = kq*nq + kq*z*dnqz;
Krq = dnqr*kq*z;
k1z = k1*n1 + k1*z*dn1z;
k1r = dn1r*k1*z;

%radial dependance wavevector
Krz1 = (k1*r^2*dn1z)/(2*z*((Zr^2*n1^2)/z^2 + 1)) - (k1*r^2*n1)/(2*z^2*((Zr^2*n1^2)/z^2 + 1)) + (k1*r^2*n1*((2*Zr^2*n1^2)/z^3 - (2*Zr^2*dn1z*n1)/z^2))/(2*z*((Zr^2*n1^2)/z^2 + 1)^2);
Krr1 = (k1*r^2*dn1r)/(2*z*((Zr^2*n1^2)/z^2 + 1)) + (k1*r*n1)/(z*((Zr^2*n1^2)/z^2 + 1)) - (Zr^2*k1*r^2*dn1r*n1^2)/(z^3*((Zr^2*n1^2)/z^2 + 1)^2);

%gouy wavevector
Kgz1 = (-1/(Zr*n1)+z*dn1z/(Zr*n1^2))/(1+(z/(n1*Zr))^2);
Kgr1 = -z/(Zr*n1^2)/(1+(z/(n1*Zr))^2)*dn1r;

%final wavevector
Kz1 = k1z + Kgz1 + Krz1;
Kr1 = k1r + Kgr1 + Krr1;


K = sqrt(Kzq^2+Krq^2) - sqrt((q*Kz1+Kzd)^2 + (q*Kr1+Krd)^2);




end
