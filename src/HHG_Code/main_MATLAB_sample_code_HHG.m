%Basic code to run the HHG simulation:


%General parameters:
	q =21;                  %Harmonic order (1:51)
	t = 0;                  %Time [s] (gaussian pulse centered in t = 0)
	alpha = 2e-14;			%Phase coefficient (taken at 2e-14 cm2/W for the whole study)
	Te = 3;	            	%Freed electron temperature [eV] (taken as 3eV for the whole study)


%Laser parameters:
	I0 = 6e13;              %Peak intensity [W/cm2]
	tp = 130e-15;			%Pulse length (FWHM)
	lambda1 = 1050e-9;		%Fundamental wavelength [m]
	R0 = 19.6e-6;			%Beam radius [m]
	f = 60e6; 	           	%Laser frequency [Hz]


%Gas parameters:
	V = 250;          		%Gas velocity [m/s]
	P = 500;                %Peak pressure [mbar]
	lp = 150e-6; 			%Interaction length (FWHM if gaussian)/ nozzle diameter [m]
	profile = 'gauss'; 	  	%Density profile ('squar' or 'gauss')
	znozzle = 0;	       	%Nozzle position [m]
	gas = 'Kr';     		%Gas ('Ar', 'Kr', 'Xe')
	xHe = 0;                %Helium fraction (0:1)


%Graph parameters:
	zmax = 1e-3; 		   	%Boundaries calculation on optical axis
	rmax = 50e-6;		   	%Boundaries calculations on nozzle axis
	nres = 50; 	  	   	%Resolution (50,100 or 200) [50 is a good approximation if ionization not to high, but 200 is needed for accurate calculation, do not go higher than 200 - or MATLAB will go out of memory]

		%The calculations are made with:
		%z	= [-zmax : (2*zmax)/nres : zmax];
		%r	= [-rmax : (2*rmax)/nres : rmax];



%Basic code to calculate the harmonic amplitude:
	%1) Calculate the dipole response amplitude along the optical axis:
	[dipqz,~,~,~] = detectdipole(I0,lambda1,znozzle,tp,R0,gas,q,t,zmax,nres,figure(1));close(1);
	
	%2) Calculate the phase matching for all z and r:
    [~,Dp,~,~,ethaft] = phase_matching(I0,V,P,lp,profile,znozzle,gas,q,f,R0,lambda1,tp,Te,alpha,zmax,rmax,nres,t,xHe);
	
	%3) Calculate the harmonic growth along the optical axis using the dipole response and the phase matching calculated before:
    [~,Iq] = amplitude(dipqz,ethaft,Dp,P,lp,profile,q,lambda1,znozzle,zmax,nres,gas);
		
	%4) If you are interrested only in the output harmonic, just take the end of the array:
    Iqz = Iq(end);
	
	%display the result
	display(Iqz)


	%To have access to other datas like electron wavefunction, wavevector mismatch, ionization fraction... Take a look into the function .m files.
	%The pdesolver fail for too high ionization fraction -> too high pressure, intensity or too low gas velocity, (I think its when ionization reach higher than 1).

	%If you are interrested in negligible absorption or perfect phasematching calculations, you can call amplitude_PMfree.m, amplitude_ABSfree.m and amplitude_ABSPM_free.m instead of amplitude.m
	%Calculation with non zero helium fraction require the call of amplitude_He.m instead of amplitude.m