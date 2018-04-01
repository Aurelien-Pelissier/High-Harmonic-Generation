%Studying the Harmonic amplitude output as a function of pressure for different harmonic order


clear all
close all

q = [21,23,25,27,29,31];
I0 = 6e13;

V  = 250;           %gas velocity
P = linspace(50,4000,100);  %maximum pressure
lp = 150e-6;        %interaction length
profile = 'gauss';  %density profile
t = 0;
alpha	= 2e-14;    %phase coefficient
tp		= 130e-15;  %pulse duration
lambda1 = 1050e-9;  %fundamental wavelength
R0		= 19.6e-6;  %beam waist
gas     = 'Kr';     %Kr for Krypton,Xe for Xenon and Ar for Argon
Te = 3;             %freed electron temperature
f = 60e6;           %pulse frequency
znozzle = 0;        %nozzle position
xHe = 0;

zmax	= 1e-3;     %graph parameter zmax
rmax	= 50e-6;    %graph parameter rmax
nres    = 200;       %graph parameter resolution

for mV = 1:length(V)
for mI = 1:length(I0)



for j = 1:length(q)
    dipqz = detectdipole(I0(mI),lambda1,znozzle,tp,R0,gas,q(j),t,zmax,nres,figure(1));
    close(1)
    for i = 1:length(P)
        [~,Dp,~,~,ethaft] = phase_matching(I0(mI),V(mV),P(i),lp,profile,znozzle,gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P(i),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqz(i,j,mI) = Iq(end);
        
        [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P(i),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas); 
        Iqzp(i,j,mI) = Iq2(end);
        
        [~,Iq3] = amplitude_ABSfree(dipqz,ethaft,Dp,P(i),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqza(i,j,mI) = Iq3(end);
        
        [~,Iq4] = amplitude_ABSPMfree(dipqz,ethaft,P(i),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqzap(i,j,mI) = Iq4(end);
    end
end
end


figure(2)
h = plot(1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10);
c = get(h,'Color');


figure(2)
semilogy(P,abs(Iqz(:,1,1)),'linewidth',2,'color',c{1})
hold on
semilogy(P,abs(Iqz(:,2,1)),'linewidth',2,'color',c{2})
semilogy(P,abs(Iqz(:,3,1)),'linewidth',2,'color',c{3})
semilogy(P,abs(Iqz(:,4,1)),'linewidth',2,'color',c{4})
semilogy(P,abs(Iqz(:,5,1)),'linewidth',2,'color',c{5})
semilogy(P,abs(Iqz(:,6,1)),'linewidth',2,'color',c{6})
% semilogy(P,abs(Iqzp(:,1,1)),'--','linewidth',2,'color',c{1})
% semilogy(P,abs(Iqzp(:,2,1)),'--','linewidth',2,'color',c{2})
% semilogy(P,abs(Iqzp(:,3,1)),'--','linewidth',2,'color',c{3})
% semilogy(P,abs(Iqzp(:,4,1)),'--','linewidth',2,'color',c{4})
% semilogy(P,abs(Iqzp(:,5,1)),'--','linewidth',2,'color',c{5})
% semilogy(P,abs(Iqzp(:,6,1)),'--','linewidth',2,'color',c{6})
hold off
legend('21th','23th','25th','27th','29th','31th','location','southeast')
ylabel('Harmonic amplitude [arb. unit]')
xlabel('Pressure [mbar]')
set(gca,'fontsize',15)
ylim([1e13 2e16])
print(sprintf('Pressure_Intensity_%.0fms.png',V(mV)),'-dpng')
end