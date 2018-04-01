%Studying the Harmonic amplitude output as a function of pressure for different interaction length
% Looking at the output with optimum nozzle position instead of just znoz = 0;


close all
clc
q = 21;
I0 = 6e13;

V  = 200;           %gas velocity
P = linspace(50,3000,30);  %maximum pressure
lp = [50e-6,125e-6,200e-6,300e-6];        %interaction length
profile = 'squar';  %density profile
t = 0;
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
nres    = 200;       %graph parameter resolution
nresnoz = 30;

zmin = -zmax;
znoz = [zmin : (zmax-zmin)/nresnoz : zmax];

for mV = 1:length(V)
for ml = 1:length(lp)



for j = 1:length(q)
    dipqz = detectdipole(I0,lambda1,0,tp,R0,gas,q(j),t,zmax,nres,figure(1));
    close(1)
    for i = 1:length(P)
        display(i)
        
        for nz = 1:length(znoz)
            [~,Dp,~,~,ethaft] = phase_matching(I0,V(mV),P(i),lp(ml),profile,znoz(nz),gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
            [~,Iq] = amplitude(dipqz,ethaft,Dp,P(i),lp(ml),profile,q(j),lambda1,znoz(nz),zmax,nres,gas);
            Iqz(nz) = Iq(end);
        end
        Iqz0(i,ml,j) = Iqz(floor(nresnoz/2)+1);
        [MaxI,index] = max(Iqz);
        Iqz1(i,ml,j) = MaxI;
        znozm(i,ml,j) = znoz(floor(index*nresnoz/nres)+1);  
    end
end
end

figure(8)
semilogy(P,abs(Iqz0(:,1,1)),'linewidth',2)
hold on
plot(P,abs(Iqz0(:,2,1)),'linewidth',2)
plot(P,abs(Iqz0(:,3,1)),'linewidth',2)
plot(P,abs(Iqz0(:,4,1)),'linewidth',2)
hold off
legend('50um','125um','200um','300um','location','southeast')
ylabel('Harmonic amplitude at znoz = 0 [arb. unit]')
xlabel('Pressure[mbar]')
print(sprintf('lp_Pressure_%.0fms.png',V),'-dpng')



figure(9)
semilogy(P,abs(Iqz1(:,1,1)),'linewidth',2)
hold on
plot(P,abs(Iqz1(:,2,1)),'linewidth',2)
plot(P,abs(Iqz1(:,3,1)),'linewidth',2)
plot(P,abs(Iqz1(:,4,1)),'linewidth',2)
hold off
legend('50um','125um','200um','300um','location','southeast')
ylabel('Max. Harmonic amplitude [arb. unit]')
xlabel('Pressure[mbar]')
print(sprintf('lp_Pressure_%.0fms.png',V),'-dpng')
end