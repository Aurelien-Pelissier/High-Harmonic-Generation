%plot the phase matching for the square and gaussian density profile


q = 21;
I0 = 6e13;

V  = 250;           %gas velocity
P = 320;  %maximum pressure
lp = 400e-6;        %interaction length
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

zmin	= -zmax;
rmin	= -rmax;
z		= [zmin : (zmax-zmin)/nres : zmax];
r		= [rmin : (rmax-rmin)/nres : rmax];





profile = 'gauss'; %density profile     
dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q,t,zmax,nres,figure(1));
close(1)
[Dk1,Dp1,~,~,ethaft] = phase_matching(I0,V,P,lp,profile,znozzle,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
[zz1,Iq1] = amplitude(dipqz,ethaft,Dp1,P,lp,profile,q,lambda1,znozzle,zmax,nres,gas);
Cl1 = pi./Dk1;
        
        
figure(11)
subplot(3,1,[1,2])
    imagesc(z*1e3,r*1e6,log10(abs(Cl1)))
    colormap(jet)
    colormap(flipud(colormap))
    caxis([-6,0])
    set(gca,'XTick',[]);
    ylabel('r[um]')
    set(gca, 'fontsize', 19);
        
subplot(3,1,3)
    plot(zz1*1e3,abs(Iq1),'linewidth',2)
    xlabel('z [mm]'), ylabel('Amplitude [arb. unit]')
    ylim([0 4e15])
    set(gca, 'fontsize', 19);
    
    
    
    
    
profile = 'squar'; %density profile     
dipqz = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q,t,zmax,nres,figure(1));
close(1)
[Dk2,Dp2,~,~,ethaft] = phase_matching(I0,V,P,lp,profile,znozzle,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t);
[zz2,Iq2] = amplitude(dipqz,ethaft,Dp1,P,lp,profile,q,lambda1,znozzle,zmax,nres,gas);
Cl2 = pi./Dk2;
        
        
figure(12)
subplot(3,1,[1,2])
    imagesc(z*1e3,r*1e6,log10(abs(Cl2)))
    colormap(jet)
    colormap(flipud(colormap))
    caxis([-6,0])
    set(gca,'XTick',[]);
    ylabel('r[um]')
    set(gca, 'fontsize', 19);
        
subplot(3,1,3)
    plot(zz2*1e3,abs(Iq2),'linewidth',2)
    xlabel('z [mm]'), ylabel('Amplitude [arb. unit]')
    ylim([0 4e15])
    set(gca, 'fontsize', 19);

