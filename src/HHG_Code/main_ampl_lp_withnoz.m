%Studying the Harmonic amplitude output as a function of interaction length for different pressure
% Looking at the output with optimum nozzle position instead of just znoz = 0;

close all
clc
q = 21;
I0 = 6e13;

V  = 200;           %gas velocity
P = [100,500,1000];  %maximum pressure
lp = linspace(10*1e-6,500*1e-6,20);        %interaction length
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
nresnoz = 30;

zmin = -zmax;
znoz = [zmin : (zmax-zmin)/nresnoz : zmax];

dipqz = zeros(length(znoz),length(znoz));

for j = 1:length(q) 
    dipqz = detectdipole(I0,lambda1,0,tp,R0,gas,q(j),t,zmax,nres,figure(1));
    close(1)
    for i = 1:length(P)
        for ml = 1:length(lp)
            display(ml)
            for nz = 1:length(znoz)
                [~,Dp,~,~,ethaft] = phase_matching(I0,V,P(i),lp(ml),profile,znoz(nz),gas,q(j),f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
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
semilogy(lp*1e6,abs(Iqz0(1,:,1)),'linewidth',2)
hold on
plot(lp*1e6,abs(Iqz0(2,:,1)),'linewidth',2)
plot(lp*1e6,abs(Iqz0(3,:,1)),'linewidth',2)
hold off
legend('100mbar','500mbar','1000mbar','location','southeast')
ylabel('Harmonic amplitude at znoz = 0 [arb. unit]')
xlabel('Interaction length [um]')
print(sprintf('lp_Pressure_%.0fms.png',V),'-dpng')



figure(9)
semilogy(lp*1e6,abs(Iqz1(1,:,1)),'linewidth',2)
hold on
plot(lp*1e6,abs(Iqz1(2,:,1)),'linewidth',2)
plot(lp*1e6,abs(Iqz1(3,:,1)),'linewidth',2)
hold off
legend('100mbar','500mbar','1000mbar','location','southeast')
ylabel('Max. Harmonic amplitude [arb. unit]')
xlabel('Interaction length [um]')
print(sprintf('lp_Pressure_%.0fms.png',V),'-dpng')