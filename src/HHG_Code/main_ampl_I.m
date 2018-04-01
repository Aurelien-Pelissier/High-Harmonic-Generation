%Studying the Harmonic amplitude output as a function of intensity for different pressure


q = 21;
I0 = linspace(1e13,1.35e14,100);

V  = 250;           %gas velocity
% V = [200,250,300]
P = [250,500,750];  %maximum pressure
lp = 150e-6;        %interaction length
profile = 'squar';  %density profile
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
nres    = 200;      %graph parameter resolution

for mV = 1:length(V)
for mP = 1:length(P)


for i = 1:length(I0)
    for j = 1:length(q)
        dipqz = detectdipole(I0(i),lambda1,znozzle,tp,R0,gas,q(j),t,zmax,nres,figure(1));
        close(1)
        [~,Dp,~,~,ethaft] = phase_matching(I0(i),V,P(mP),lp,profile,znozzle,gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P(mP),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqz(i,mP,j) = Iq(end);
        
        [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P(mP),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas); 
        Iqzp(i,mP,j) = Iq2(end);
        
        [~,Iq3] = amplitude_ABSfree(dipqz,ethaft,Dp,P(mP),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqza(i,mP,j) = Iq3(end);
        
        [~,Iq4] = amplitude_ABSPMfree(dipqz,ethaft,P(mP),lp,profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqzap(i,mP,j) = Iq4(end);
    end
end
end


figure(9)
%Harmonic amplitude as a function of intensity for different pressure
%With free absorption
h = plot(1:10,1:10,1:10,1:10,1:10,1:10);
c = get(h,'Color');


figure(9)
plot(I0(1:end-8),abs(Iqz(1:end-8,1,1)),'linewidth',2,'color',c{1})
hold on

plot(I0(1:end-8),abs(Iqz(1:end-8,2,1)),'linewidth',2,'color',c{2})
plot(I0(1:end-8),abs(Iqz(1:end-8,3,1)),'linewidth',2,'color',c{3})
plot(I0(1:end-8),abs(Iqzp(1:end-8,1,1)),'--','linewidth',1,'color',c{1})
plot(I0(1:end-8),abs(Iqzp(1:end-8,2,1)),'--','linewidth',1,'color',c{2})
plot(I0(1:end-8),abs(Iqzp(1:end-8,3,1)),'--','linewidth',1,'color',c{3})

hold off
legend('250mbar','500mbar','750mbar')
xlabel('Peak Intensity [W/cm2]')
ylabel('Harmonic amplitude [arb. unit]')
set(gca,'fontsize',15)
print(sprintf('Intensity_Pressure_%.0fms.png',V(mV)),'-dpng')



figure(10)
%Harmonic amplitude as a function of intensity for different pressure, semilog
%With free absorption
semilogy(I0(1:end-8),abs(Iqz(1:end-8,1,1)),'linewidth',2,'color',c{1})
hold on

semilogy(I0(1:end-8),abs(Iqz(1:end-8,2,1)),'linewidth',2,'color',c{2})
semilogy(I0(1:end-8),abs(Iqz(1:end-8,3,1)),'linewidth',2,'color',c{3})
semilogy(I0(1:end-8),abs(Iqzp(1:end-8,1,1)),'--','linewidth',1,'color',c{1})
semilogy(I0(1:end-8),abs(Iqzp(1:end-8,2,1)),'--','linewidth',1,'color',c{2})
semilogy(I0(1:end-8),abs(Iqzp(1:end-8,3,1)),'--','linewidth',1,'color',c{3})

hold off
legend('250mbar','500mbar','750mbar','location','southeast')
xlabel('Peak Intensity [W/cm2]')
ylabel('Harmonic amplitude [arb. unit]')
set(gca,'fontsize',15)
print(sprintf('Intensity_Pressure_log_%.0fms.png',V(mV)),'-dpng')



figure(11)
%Harmonic amplitude as a function of intensity for different pressure
%With perfect phase matching
semilogy(I0(1:end-8),abs(Iqz(1:end-8,1,1)),'linewidth',2,'color',c{1})
hold on

semilogy(I0(1:end-8),abs(Iqz(1:end-8,2,1)),'linewidth',2,'color',c{2})
semilogy(I0(1:end-8),abs(Iqz(1:end-8,3,1)),'linewidth',2,'color',c{3})
semilogy(I0(1:end-8),abs(Iqza(1:end-8,1,1)),'--','linewidth',1,'color',c{1})
semilogy(I0(1:end-8),abs(Iqza(1:end-8,2,1)),'--','linewidth',1,'color',c{2})
semilogy(I0(1:end-8),abs(Iqza(1:end-8,3,1)),'--','linewidth',1,'color',c{3})

hold off
legend('250mbar','500mbar','750mbar','location','southeast')
xlabel('Peak Intensity [W/cm2]')
ylabel('Harmonic amplitude [arb. unit]')
set(gca,'fontsize',15)
print(sprintf('Intensity_Pressure_log_%.0fms.png',V(mV)),'-dpng')




figure(12)
%Harmonic amplitude as a function of intensity for different pressure, semilog
%With perfect phase matching
semilogy(I0(1:end-8),abs(Iqz(1:end-8,1,1)),'linewidth',2,'color',c{1})
hold on

semilogy(I0(1:end-8),abs(Iqz(1:end-8,2,1)),'linewidth',2,'color',c{2})
semilogy(I0(1:end-8),abs(Iqz(1:end-8,3,1)),'linewidth',2,'color',c{3})
semilogy(I0(1:end-8),abs(Iqzap(1:end-8,1,1)),'--','linewidth',1,'color',c{1})
semilogy(I0(1:end-8),abs(Iqzap(1:end-8,2,1)),'--','linewidth',1,'color',c{2})
semilogy(I0(1:end-8),abs(Iqzap(1:end-8,3,1)),'--','linewidth',1,'color',c{3})

hold off
legend('250mbar','500mbar','750mbar','location','southeast')
xlabel('Peak Intensity [W/cm2]')
ylabel('Harmonic amplitude [arb. unit]')
set(gca,'fontsize',15)
print(sprintf('Intensity_Pressure_log_%.0fms.png',V(mV)),'-dpng')
end