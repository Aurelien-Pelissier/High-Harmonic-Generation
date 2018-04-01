function [I,dIz,dIr] = Igauss(I0,tp,R0,lambda1,z,r,t)
%return gaussian beam Intesity and its derivative at position (z,r) and time t%

%experiment cst
c		= 299792458;
w1		= c*2*pi/lambda1;
k1		= w1/c;																	
Zr		= k1*R0^2/2;

%function
sigma	= tp/(2*sqrt(2*log(2)));
Rzz		= R0*sqrt(1+(z/Zr)^2);
I		= I0*(R0/Rzz)^2*exp((-2*r^2)/(Rzz^2))*exp(-t^2/(2*sigma^2));

%gradient
dIz		= (4*I0*r^2*z*exp(-t^2/(2*sigma^2))*exp(-(2*r^2)/(R0^2*(z^2/Zr^2 + 1))))/(R0^2*Zr^2*(z^2/Zr^2 + 1)^3) - (2*I0*z*exp(-t^2/(2*sigma^2))*exp(-(2*r^2)/(R0^2*(z^2/Zr^2 + 1))))/(Zr^2*(z^2/Zr^2 + 1)^2);
dIr		= -(4*I0*r*exp(-t^2/(2*sigma^2))*exp(-(2*r^2)/(R0^2*(z^2/Zr^2 + 1))))/(R0^2*(z^2/Zr^2 + 1)^2);

end


