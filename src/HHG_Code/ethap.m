function etha = ethap(I0,tp,lambda1,gas,t1)
% calculate the ionization fraction at time t1 of the pulse,  intensity I0, and FWHM tp

if t1 <=-tp %bigger than tp is equivalent to inf
    t = -tp;
elseif t1 >=tp
    t = tp;
else
    t = t1;
end


%univers cst
c0		= 299792458;
me		= 9.109e-31;
hbar	= 1.055e-34; 	
e		= 1.602e-19;
etha0	= 377; %impedence of free space

%experiment cst
w1		= c0*2*pi/lambda1;
k1		= w1/c0;		


if I0 <= 1e13 && numel(gas) == 2 %we know its 0 (<3e-5), its to avoide useless calculations
    etha = 0;
else

%gaussian beam
sigma	= tp/(2*sqrt(2*log(2)));
I		= @(t)I0*exp(-t.^2/(2*sigma^2));


%function
Ip      = Ipot(gas); %eV
%Up      = (etha0*e^2)/(2*me*w1^2)*I*6.242e18*1e4
Up      = @(t)0.0933*lambda1^2.*I(t);
gamma   = @(t)sqrt(Ip./(2.*Up(t)));


%exponential factor
a = @(t)1 + gamma(t).^2 - sin(w1.*t).^2;
b = @(t)sqrt(a(t).^2 + 4.*gamma(t).^2.*sin(w1.*t).^2);
c = @(t)sqrt((sqrt((b(t)+a(t))/2)+gamma(t)).^2+(sqrt((b(t)-a(t))./2)+sin(abs(w1.*t))).^2);
d = @(t)(gamma(t).^2 + sin(w1.*t).^2 + 1/2).*log(c(t)) - 3.*sqrt(b(t)-a(t))./2./sqrt(2).*sin(abs(w1.*t))-gamma(t).*sqrt(b(t)+a(t))./2./sqrt(2);
kappa = 1e4*(2*etha0*e^2/hbar/me/(w1)^3);
W = @(t)(exp(-kappa.*d(t).*I(t)));

%pre-exponential factor
Npre = @(t)Nt(I0*1e4,tp,t,lambda1,Ip);


Wf = @(t)Npre(t).*W(t);



%integration
etha = 1 - exp(-integral(Wf,-tp,t));
end

end

