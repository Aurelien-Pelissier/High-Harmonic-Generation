# High-Harmonic-Generation
<img align="left" src="https://raw.githubusercontent.com/Aurelien-Pelissier/High-Harmonic-Generation/master/report/HHG.png" width=200>
High harmonic generation (HHG) refers to the process of creating vacuum (VUV) or extreme (XUV) ultraviolet light through a nonlinear interaction of an intense laser field with a gas target. This repository contains the source code to calculate the high harmonics amplitude accounting for all involved physical phenomenons (Quantum mechanic, Supersonic flow physics, Plasma physics and Nonlinear optics). More specifically, we consider supersonic gas flow at the nozzle outlet, ions dynamic in the plasma, quantum atomic response of freed electrons, phase matching and absorption.




## Involved parameters

#### Laser parameters

- Peak intensity
- Pulse length (FWHM)
- Wavelength
- Beam Radius
- Frequency
- Harmonic order

#### Gas parameters
- Gas Velocity
- Peak preassure
- Nozzle diameter
- Nozzle position
- Gas (Krypton, Argon, Xenon)





## Running the code

In the folder, all the MATLAB files starting with "main" are different simulations of the HHG process to study the influence of various parameters such as intensity, pressure or gas velocity. However, the basic code to run the simulation is the following:


```matlab
%General parameters:
	q =21;                  	%Harmonic order (1:51)
	t = 0;                  	%Time [s] (gaussian pulse centered in t = 0)
	alpha = 2e-14;			%Phase coefficient (taken at 2e-14 cm2/W for the whole study)
	Te = 3;	            		%Freed electron temperature [eV] (taken as 3eV for the whole study)


%Laser parameters:
	I0 = 6e13;              	%Peak intensity [W/cm2]
	tp = 130e-15;			%Pulse length (FWHM)
	lambda1 = 1050e-9;		%Fundamental wavelength [m]
	R0 = 19.6e-6;			%Beam radius [m]
	f = 60e6; 	           	%Laser frequency [Hz]


%Gas parameters:
	V = 250;          		%Gas velocity [m/s]
	P = 500;                	%Peak pressure [mbar]
	lp = 150e-6; 			%Interaction length (FWHM if gaussian)/ nozzle diameter [m]
	profile = 'gauss'; 	  	%Density profile ('squar' or 'gauss')
	znozzle = 0;	       		%Nozzle position [m]
	gas = 'Kr';     		%Gas ('Ar', 'Kr', 'Xe')
	xHe = 0;                	%Helium fraction (0:1)


%Graph parameters:
	zmax = 1e-3; 		   	%Boundaries calculation on optical axis
	rmax = 50e-6;		   	%Boundaries calculations on nozzle axis
	nres = 200; 	  	   	%Resolution
	
		%The calculation is made with:
		%z	= [-zmax : (2*zmax)/nres : zmax];
		%r	= [-rmax : (2*rmax)/nres : rmax];
```

The above parameters are realistic values that we would use when generating HHG experimentally.
Now that we have set all the important input parameters for the simulation, we can compute the harmonic output with:

```matlab

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
```
Note that this program return one number, which is the harmonic amplitude output in arbitrary unit. While this number itself has no significant meaning, it is interresting to modify the input parameters to see how the output is modified.

To have access to other datas like electron wavefunction, wavevector mismatch, ionization fraction... You can take a look into the function .m files. If you are interrested in negligible absorption or perfect phasematching calculations, you can call amplitude_PMfree.m, amplitude_ABSfree.m and amplitude_ABSPM_free.m instead of amplitude.m. Calculation with non zero helium fraction require the call of amplitude_He.m instead of amplitude.m


#### Code structure

	main
	------detectdipole					%check if dipole file already exist
	------------dipole					%calculate the dipole response
	------------------potde					%potential depth of Coulomb potential

	------phasematching					%calculate phasematching
	------------detectetha  				%check if the ionization (single pulse) file already exist
	------------------ethap 				%calculate ioniation fraction after one pulse
	------------------------Ipot				%ionization potentiel
	------------------------Nt				%pre-exponential factor of the Yudin rate

	------------detectethan  				%check if the ionization (multiple pulse) file already exist
	------------------etharz 				%calculate ioniation fraction after n pulse
	------------------------ethabp 				%compute the decay between 2 pulses
	------------------------------solvepde 			%solve the pde for a fixed z
	------------------------------------Ndens		%Number density at 1atm of gas
	------------------------------------mobility		%mobility of ions
	------------------------------------Press		%Pressure distribution

	------phase						%calculate phase matching
	------------refractive					%compute refractive index
	------------Igauss					%intensity gaussian distribution
	------------Press					%pressure distribution
	------------ioniz					%ionization fraction

	------amplitude						%calculate final harmonic amplitude
	------------absorb					%absorption by the gas













## Results




## Reports
