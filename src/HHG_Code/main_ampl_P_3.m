%Studying the Harmonic amplitude output as a function of pressure for different interaction length


clear all
clc

q = 21;
I0 = 6e13;

V  = 200;           %gas velocity
P = linspace(50,3000,100);  %maximum pressure
lp = [50e-6,125e-6,200e-6,300e-6];        %interaction length
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
xHe = 0;

zmax	= 1e-3;     %graph parameter zmax
rmax	= 50e-6;    %graph parameter rmax
nres    = 200;       %graph parameter resolution

for mV = 1:length(V)
for ml = 1:length(lp)



for j = 1:length(q)
    dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q(j),t,zmax,nres,figure(1));
    close(1)
    for i = 1:length(P)
        [~,Dp,~,~,ethaft] = phase_matching(I0,V(mV),P(i),lp(ml),profile,znozzle,gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqz(i,j,ml) = Iq(end);
        
        [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas); 
        Iqzp(i,j,ml) = Iq2(end);
        
        [~,Iq3] = amplitude_ABSfree(dipqz,ethaft,Dp,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqza(i,j,ml) = Iq3(end);
        
        [~,Iq4] = amplitude_ABSPMfree(dipqz,ethaft,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqzap(i,j,ml) = Iq4(end);
    end
end
end

figure(2)
h = plot(1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10);
c = get(h,'Color');

figure(2)
semilogy(P,abs(Iqz(:,1,1)),'linewidth',2)
hold on
semilogy(P,abs(Iqz(:,1,2)),'linewidth',2)
semilogy(P,abs(Iqz(:,1,3)),'linewidth',2)
semilogy(P,abs(Iqz(:,1,4)),'linewidth',2)
semilogy(P,abs(Iqzp(:,1,1)),'--','linewidth',2,'color',c{1})
semilogy(P,abs(Iqzp(:,1,2)),'--','linewidth',2,'color',c{2})
semilogy(P,abs(Iqzp(:,1,3)),'--','linewidth',2,'color',c{3})
semilogy(P,abs(Iqzp(:,1,4)),'--','linewidth',2,'color',c{4})
hold off
legend('50um','125um','200um','300um','location','southwest')
legend boxoff
ylabel('Harmonic amplitude [arb. unit]')
ylim([6e11 1e16])
xlabel('Pressure [mbar]')
set(gca,'fontsize',15)
print(sprintf('Pressure_Intensity_%.0fms.png',V(mV)),'-dpng')
end