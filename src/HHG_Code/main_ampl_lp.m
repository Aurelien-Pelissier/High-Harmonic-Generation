%Studying the Harmonic amplitude output as a function of interaction length for different pressure

clear all
clc

q = 21;
I0 = 6e13;

V  = 200;           %gas velocity
P = [100,250,500,1000];
lp = linspace(22.7*1e-6,500*1e-6,106);        %interaction length
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





for j = 1:length(q)
    dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q(j),t,zmax,nres,figure(1));
    close(1)
    for i = 1:length(P)
    for ml = 1:length(lp)
        [~,Dp,~,~,ethaft] = phase_matching(I0,V,P(i),lp(ml),profile,znozzle,gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqz(i,ml,j) = Iq(end);
        
        [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas); 
        Iqzp(i,ml,j) = Iq2(end);
        
        [~,Iq3] = amplitude_ABSfree(dipqz,ethaft,Dp,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqza(i,ml,j) = Iq3(end);
        
        [~,Iq4] = amplitude_ABSPMfree(dipqz,ethaft,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqzap(i,ml,j) = Iq4(end);
    end
    end
end


figure(8)
h = plot(1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10);
c = get(h,'Color');

figure(8)
%Studying the Harmonic amplitude output as a function of interaction length
%with perfect phaseatching
semilogy(lp*1e6,abs(Iqz(1,:,1)),'linewidth',2,'color',c{1})
hold on
semilogy(lp*1e6,abs(Iqz(2,:,1)),'linewidth',2,'color',c{2})
semilogy(lp*1e6,abs(Iqz(3,:,1)),'linewidth',2,'color',c{3})
semilogy(lp*1e6,abs(Iqz(4,:,1)),'linewidth',2,'color',c{4})
semilogy(lp*1e6,abs(Iqzp(1,:,1)),'--','linewidth',2,'color',c{1})
semilogy(lp*1e6,abs(Iqzp(2,:,1)),'--','linewidth',2,'color',c{2})
semilogy(lp*1e6,abs(Iqzp(3,:,1)),'--','linewidth',2,'color',c{3})
semilogy(lp*1e6,abs(Iqzp(4,:,1)),'--','linewidth',2,'color',c{4})
hold off
legend('100mbar','250mbar','500mbar','1000mbar')
legend boxoff
ylabel('Harmonic amplitude [arb. unit]')
xlabel('Interaction length [um]')
set(gca,'fontsize',15)
% ylim([1e13 1e19])
print(sprintf('lp_Pressure_%.0fms.png',V),'-dpng')



figure(9)
%Studying the Harmonic amplitude output as a function of interaction length
%with negligible absorption
semilogy(lp*1e6,abs(Iqz(1,:,1)),'linewidth',2,'color',c{1})
hold on
semilogy(lp*1e6,abs(Iqz(2,:,1)),'linewidth',2,'color',c{2})
semilogy(lp*1e6,abs(Iqz(3,:,1)),'linewidth',2,'color',c{3})
semilogy(lp*1e6,abs(Iqza(1,:,1)),'--','linewidth',2,'color',c{1})
semilogy(lp*1e6,abs(Iqza(2,:,1)),'--','linewidth',2,'color',c{2})
semilogy(lp*1e6,abs(Iqza(3,:,1)),'--','linewidth',2,'color',c{3})
hold off
legend('100mbar','250mbar','500mbar','1000mbar')
ylabel('Harmonic amplitude [arb. unit]')
xlabel('Interaction length [um]')
ylim([1e13 1e19])
print(sprintf('lp_Pressure_%.0fms.png',V),'-dpng')



figure(10)
%Studying the Harmonic amplitude output as a function of interaction length
%with perfect phaseatching and negligible absorption
semilogy(lp*1e6,abs(Iqz(1,:,1)),'linewidth',2,'color',c{1})
hold on
semilogy(lp*1e6,abs(Iqz(2,:,1)),'linewidth',2,'color',c{2})
semilogy(lp*1e6,abs(Iqz(3,:,1)),'linewidth',2,'color',c{3})
semilogy(lp*1e6,abs(Iqzap(1,:,1)),'--','linewidth',2,'color',c{1})
semilogy(lp*1e6,abs(Iqzap(2,:,1)),'--','linewidth',2,'color',c{2})
semilogy(lp*1e6,abs(Iqzap(3,:,1)),'--','linewidth',2,'color',c{3})
hold off
legend('100mbar','250mbar','500mbar','1000mbar')
ylabel('Harmonic amplitude [arb. unit]')
xlabel('Interaction length [um]')
ylim([1e13 1e19])
print(sprintf('lp_Pressure_%.0fms.png',V),'-dpng')
