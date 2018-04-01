%Studying the Harmonic amplitude output as a function of intensity for different harmonic order

q = [21,23,25,27,29,31];
I0 = linspace(1e13,1.35e14,100);

V  = 250;           %gas velocity
% V = [200,250,300]
P = 500;
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
    end
end
end

figure(15)
h = plot(1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10,1:10);
c = get(h,'Color');


figure(15)
%Harmonic amplitude as a function of intensity for different harmonic order
%With perfect phasematching
semilogy(I0,abs(Iqz(:,1,1)),'linewidth',2)
hold on

semilogy(I0,abs(Iqz(:,1,2)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,3)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,4)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,5)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,6)),'linewidth',2)

semilogy(I0,abs(Iqzp(:,1,1)),'--','linewidth',2,'color',c{1})
semilogy(I0,abs(Iqzp(:,1,2)),'--','linewidth',2,'color',c{2})
semilogy(I0,abs(Iqzp(:,1,3)),'--','linewidth',2,'color',c{3})
semilogy(I0,abs(Iqzp(:,1,4)),'--','linewidth',2,'color',c{4})
semilogy(I0,abs(Iqzp(:,1,5)),'--','linewidth',2,'color',c{5})
semilogy(I0,abs(Iqzp(:,1,6)),'--','linewidth',2,'color',c{6})

hold off
legend('21th','23th','25th','27th','29th','31th','location','southeast')
xlabel('Peak Intensity [W/cm2]')
ylabel('Harmonic amplitude [arb. unit]')
set(gca,'fontsize',15)
print(sprintf('Intensity_Pressure_log_%.0fms.png',V(mV)),'-dpng')



figure(16)
%Harmonic amplitude as a function of intensity for different harmonic order, semilog
semilogy(I0,abs(Iqz(:,1,1)),'linewidth',2)
hold on

semilogy(I0,abs(Iqz(:,1,2)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,3)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,4)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,5)),'linewidth',2)
semilogy(I0,abs(Iqz(:,1,6)),'linewidth',2)

hold off
legend('21th','23th','25th','27th','29th','31th','location','southeast')
xlabel('Peak Intensity [W/cm2]')
ylabel('Harmonic amplitude [arb. unit]')
set(gca,'fontsize',15)
print(sprintf('Intensity_Pressure_log_%.0fms.png',V(mV)),'-dpng')



figure(17)
%Harmonic amplitude as a function of intensity for different harmonic order
plot(I0,abs(Iqz(:,1,1)),'linewidth',2)
hold on

plot(I0,abs(Iqz(:,1,2)),'linewidth',2)
plot(I0,abs(Iqz(:,1,3)),'linewidth',2)
plot(I0,abs(Iqz(:,1,4)),'linewidth',2)
plot(I0,abs(Iqz(:,1,5)),'linewidth',2)
plot(I0,abs(Iqz(:,1,6)),'linewidth',2)

hold off
legend('21th','23th','25th','27th','29th','31th','location','southeast')
xlabel('Peak Intensity [W/cm2]')
ylabel('Harmonic amplitude [arb. unit]')
set(gca,'fontsize',15)
print(sprintf('Intensity_Pressure_%.0fms.png',V(mV)),'-dpng')
end