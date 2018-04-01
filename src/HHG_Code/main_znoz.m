%Studying the harmonic amplitude output as a function of nozzle position for different interaction length


clear all
clc

q = 21;
I0 = 7e13;

V  = 250;           %gas velocity
P = 100;            %gas pressure
profile = 'gauss';  %density profile
t = 0;
lp = [50,100,150,200].*1e-6;        %interaction length
alpha	= 2e-14;    %phase coefficient
tp		= 130e-15;  %pulse duration
lambda1 = 1050e-9;  %fundamental wavelength
R0		= 19.6e-6;  %beam waist
gas     = 'Kr';     %Kr for Krypton and Xe for Xenon
Te = 3;             %freed electron temperature
f = 60e6;           %pulse frequency
xHe = 0;

zmax	= 1e-3;     %graph parameter zmax
rmax	= 50e-6;    %graph parameter rmax
nres    = 200;      %graph parameter resolution

znozzle = linspace(-zmax, zmax, 100);        %nozzle position

for ml = 1:length(lp)
for i = 1:length(znozzle)
    dipqz = detectdipole(I0,lambda1,znozzle(i),tp,R0,gas,q,t,zmax,nres,figure(1));
    close(1)
        [~,Dp,~,~,ethaft] = phase_matching(I0,V,P,lp(ml),profile,znozzle(i),gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P,lp(ml),profile,q,lambda1,znozzle(i),zmax,nres,gas);
        Iqz(i,ml) = Iq(end);
        
        [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P,lp(ml),profile,q,lambda1,znozzle(i),zmax,nres,gas);
        Iqzp(i,ml) = Iq2(end);
        
        [~,Iq3] = amplitude_ABSfree(dipqz,ethaft,Dp,P,lp(ml),profile,q,lambda1,znozzle(i),zmax,nres,gas);
        Iqza(i,ml) = Iq3(end);
        
        [~,Iq4] = amplitude_ABSPMfree(dipqz,ethaft,P,lp(ml),profile,q,lambda1,znozzle(i),zmax,nres,gas);
        Iqzap(i,ml) = Iq4(end);
end
end

for ml = 1:length(lp)
    figure(ml)
    plot(znozzle*1e3,abs(Iqz(:,ml)),'linewidth',3,'color','green')
    hold on
    plot(znozzle*1e3,abs(Iqzp(:,ml)),'--','linewidth',3,'color','green')
    hold off
    ylabel('Harmonic amplitude')
    xlabel('nozzle position [mm]')
    ylim([0 6.4e14])
    set(gca,'fontsize',25)
end
