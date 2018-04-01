%Studying the Harmonic amplitude output as a function of pressure for different helium fraction


clear all
clc

q = 21;
I0 = 7e13;

Vgas  = 194;          %gas sound velocity of Krypton
VHe = 880;           %gas sound velocity of Helium
xHe = linspace(0,90e-2,50);  %helium fraction

lp = 150e-6;        %interaction length
profile = 'gauss';  %density profile
t = 0;
alpha	= 2e-14;    %phase coefficient
tp		= 130e-15;  %pulse duration
lambda1 = 1050e-9;  %fundamental wavelength
R0		= 19.6e-6;  %beam waist
gas     = 'Kr';     %Kr for Krypton and Xe for Xenon
Te = 3;             %freed electron temperature
f = 60e6;           %pulse frequency
znozzle = 0;        %nozzle position

zmax	= 1e-3;     %graph parameter zmax
rmax	= 50e-6;    %graph parameter rmax
nres    = 200;      %graph parameter resolution
nresP = 50;

V = Vgas.*(1-xHe) + VHe.*xHe;
Ptot = linspace(50,3000,nresP);

for mx = 1:length(xHe)
    for mP = 1:length(Ptot)
        dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q,t,zmax,nres,figure(1));
        close(1)
        [~,Dp,~,~,ethaft] = phase_matching(I0,V(mx),Ptot(mP)*(1-xHe(mx)),lp,profile,znozzle,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,0);
        [~,Iq] = amplitude_He(dipqz,ethaft,Dp,Ptot(mP),lp,profile,q,lambda1,znozzle,zmax,nres,gas,xHe(mx));
        Iqz(mx,mP) = Iq(end);
    end
    [Iqz_m(mx),indexP] = max(Iqz(mx,:));
    P_opt(mx) = Ptot(indexP);
end

figure(2)
plot(xHe*1e2,abs(Iqz_m),'linewidth',2)
xlabel('Helium fraction [%]')
ylabel('Harmonic amplitude at P_{opt}[arb. unit]')
set(gca,'fontsize',15)

figure(3)
plot(xHe*1e2,P_opt,'linewidth',2)
xlabel('Helium fraction [%]')
ylabel('Optimum total pressure [mbar]')
set(gca,'fontsize',15)

figure(4)
plot(Ptot,abs(Iqz(1,:)),'linewidth',2)
hold on
for mx = [10,20,30,40,50]
    plot(Ptot,abs(Iqz(mx,:)),'linewidth',2)
end
hold off
xlabel('Total pressure [mbar]')
ylabel('Harmonic amplitude')
legend('0%' ,sprintf('%.0f%%',xHe(10)*100) ,sprintf('%.0f%%',xHe(20)*100) ,sprintf('%.0f%%',xHe(30)*100) ,sprintf('%.0f%%',xHe(40)*100) ,sprintf('%.0f%%',xHe(50)*100))
set(gca,'fontsize',15)