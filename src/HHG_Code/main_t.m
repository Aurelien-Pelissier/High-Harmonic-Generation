%Studying the generation of high Harmonic through a ultrafast pulse


clear all
clc

q = 21;
I0 = 7e13;

V  = 250;           %gas velocity
P = 500;            %gas pressure
profile = 'gauss';  %density profile
lp = 150*1e-6;      %interaction length
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
nres    = 50;      %graph parameter resolution
t = linspace(-tp,tp,50); %time
znozzle = 0;        %nozzle position

for ml = 1:length(lp)
for i = 1:length(t)
    dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q,t(i),zmax,nres,figure(1));
    close(1)
        [~,Dp,~,~,ethaft] = phase_matching(I0,V,P,lp(ml),profile,znozzle,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t(i),xHe);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P,lp(ml),profile,q,lambda1,znozzle,zmax,nres,gas);
        Iqz(i,ml) = Iq(end);
        
        [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P,lp(ml),profile,q,lambda1,znozzle,zmax,nres,gas);
        Iqzp(i,ml) = Iq2(end);
        
        [~,Iq3] = amplitude_ABSfree(dipqz,ethaft,Dp,P,lp(ml),profile,q,lambda1,znozzle,zmax,nres,gas);
        Iqza(i,ml) = Iq3(end);
        
        [~,Iq4] = amplitude_ABSPMfree(dipqz,ethaft,P,lp(ml),profile,q,lambda1,znozzle,zmax,nres,gas);
        Iqzap(i,ml) = Iq4(end);
end
end


for ml = 1:length(lp)  
    
    for i = 1:length(t)
        Ig(i) = max(abs(Iqz(:,ml)))*Igauss(1,tp,R0,lambda1,znozzle,0,t(i));
    end
    figure(ml)
    plot(t*1e15,abs(Iqz(:,ml)),'linewidth',2,'color','blue')
    hold on
    plot(t*1e15,abs(Iqzp(:,ml))/max(abs(Iqzp(:,ml)))*max(abs(Iqz(:,ml))),'--','linewidth',2,'color','blue')
    plot(t*1e15,Ig,'--','linewidth',1.5,'Color',[0.6,0.6,0.6])
%     plot([0,0],[1.1*min(abs(Iqz(:,ml))),1.1*max(abs(Iqz(:,ml)))],'--','linewidth',1,'Color',[0.4,0.4,0.4])
    hold off
    legend('a','Dipole response','Laser pulse')
%     legend boxoff
    ylabel('Harmonic amplitude')
    xlabel('time [fs]')
    xlim([-tp tp]*1e15)
    ylim([0 3e15])
    set(gca,'fontsize',16)
end
