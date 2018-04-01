%Studying the Harmonic amplitude output as a function of interaction length for square and gaussian density profile
% also calculating the optimum pressure

q = 21;
I0 = 6e13;

V  = 250;           %gas velocity
P = linspace(50,2000,30);  %maximum pressure
lp = linspace(20*1e-6,500*1e-6,50);        %interaction length
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

    
profile = 'squar';  %density profile   
for ml = 1:length(lp)
for j = 1:length(q)
    dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q(j),t,zmax,nres,figure(1));
    close(1)
    for i = 1:length(P)
        [~,Dp,~,~,ethaft] = phase_matching(I0,V(mV),P(i),lp(ml),profile,znozzle,gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqz(i) = Iq(end);
    end
     [Iqms(ml,j),indexP] = max(Iqz);
     Pms(ml,j) = P(indexP);
     
     [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P(indexP),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas); 
     Iqmfs(ml,j) = Iq2(end);
end
end


profile = 'gauss';  %density profile   
for ml = 1:length(lp)
for j = 1:length(q)
    dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q(j),t,zmax,nres,figure(1));
    close(1)
    for i = 1:length(P)
        [~,Dp,~,~,ethaft] = phase_matching(I0,V(mV),P(i),lp(ml),profile,znozzle,gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t);
        [~,Iq] = amplitude(dipqz,ethaft,Dp,P(i),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas);
        Iqz(i) = Iq(end);
    end
     [Iqmg(ml,j),indexP] = max(Iqz);
     Pmg(ml,j) = P(indexP);
     
     [~,Iq2] = amplitude_PMfree(dipqz,ethaft,P(indexP),lp(ml),profile,q(j),lambda1,znozzle,zmax,nres,gas); 
     Iqmfg(ml,j) = Iq2(end);
end
end

figure(2)
semilogy(lp*1e6,abs(Iqmg(:,1)),'linewidth',2,'color','red')
hold on
semilogy(lp*1e6,abs(Iqms(:,1)),'linewidth',2,'color','blue')
semilogy(lp*1e6,abs(Iqmfg(:,1)),'--','linewidth',2,'color','red')
semilogy(lp*1e6,abs(Iqmfs(:,1)),'--','linewidth',2,'color','blue')
hold off
legend('Gaussian','Square','location','southwest')
legend boxoff
ylabel('Max_P Harmonic amplitude [arb. unit]')
xlabel('Interaction length [um]')
set(gca,'fontsize',15)
print(sprintf('Pressure_lp_%.0fms.png',V(mV)),'-dpng')


figure(3)
semilogy(lp*1e6,Pmg(:,1),'linewidth',2,'color','red')
hold on
semilogy(lp*1e6,Pms(:,1),'linewidth',2,'color','blue')
hold off
legend('Gaussian','Square','location','southwest')
legend boxoff
ylabel('Optimum pressure')
xlabel('Interaction length [um]')
set(gca,'fontsize',15)
end