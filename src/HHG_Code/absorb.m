function sigma = absorb(lambda1,q,gas)
%lambda = fundamental wavelength
%q = harmonic order
%gas = rare gas used for HHG
%return the absorption cross section sigma in m2

h = 6.62607004e-34; %Js
c = 299792458; %m/s
e = 1.60217662e-19; %C
re = 2.8179e-15; %classical radius of electron
Natm = 2.7e25; %gas density at STP

%%reading the data
if gas == 'Ar'
    myfile = fullfile('absorption','Argon_Xray.txt');
    Cst = 10^(18.0); %arbitrary cst to calculate absorption after measurement
elseif gas == 'Kr'
    myfile = fullfile('absorption','Krypton_Xray.txt');
    Cst = 10^(18.1);
elseif gas == 'Xe'
    myfile = fullfile('absorption','Xenon_Xray.txt');
    Cst = 10^(18.3);   
elseif gas == 'He'
    myfile = fullfile('absorption','Helium_Xray.txt');
    Cst = 10^(18);
end

A = dlmread(myfile);
ev = A(:,1);
lambda = h*c/e./ev;
f1 = A(:,2);
f2 = A(:,3);

refra = 1/2/pi*re*Natm.*lambda.^2.*(f1+i*f2); %refractive index at STP

a = imag(refra);

%%indexing the absorption
lambdaq = lambda1/q;
if lambdaq <= 80e-9
    [~,index] = min(abs(lambdaq-lambda));
    kappa = a(index);
else
    kappa =  lambdaq.^(3)*Cst;
end

alpha = 4*pi*kappa/lambdaq; %absorption coefficient at STP
sigma = alpha/Natm;

end

