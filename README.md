# High-Harmonic-Generation
<img align="left" src="https://raw.githubusercontent.com/Aurelien-Pelissier/High-Harmonic-Generation/master/report/HHG.png" width=200>
High harmonic generation (HHG) refers to the process of creating vacuum (VUV) or extreme (XUV) ultraviolet light through a nonlinear interaction of an intense laser field with a gas target. This repository contains the source code to calculate the high harmonics amplitude accounting for all involved physical phenomenons (Quantum mechanic, Supersonic flow physics, Plasma physics and Nonlinear optics). More specifically, we consider supersonic gas flow at the nozzle outlet, ions dynamic in the plasma, quantum atomic response of freed electrons, phase matching and absorption.


## Running the code

#### Graphical User Interface (GUI)
The GUI can be run with `src/HHG_GUI/main.m`, It is a user-friendly graphical interface which allow users to try different values for the input parameters and check for their impact on the HHG process. Particularily, the GUI contains information about ionization, phase matching and dipole response.

Inserti image here


Note that, while the GUI might be slow at the beginning, all calculation are saved in txt files so the calculation will be much faster the next time you run the simulation with the same parameters.


#### The code
The source code for the HHG process can be found in `src/HHG_Code/`. Inside the folder, all the MATLAB files starting with `main_[..].m` are simulations of the HHG process to study the influence of various parameters such as intensity, pressure or gas velocity. They all use the same basic code to to run the simulation, which is the following:


```matlab
%General parameters:
	q = 21;                  	%Harmonic order (integer between 1 and 51)
	t = 0;                  	%Time [s] (gaussian pulse centered in t = 0)
	alpha = 2e-14;			%Phase coefficient [cm2/W] (taken at 2e-14 cm2/W for the whole study)
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
	xHe = 0;                	%Helium fraction (between 0 and 1)

%Graph parameters:
	zmax = 1e-3; 		   	%Boundaries calculation on optical axis [m]
	rmax = 50e-6;		   	%Boundaries calculations on nozzle axis [m]
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
Note that this code will return one number, which is the harmonic amplitude output in arbitrary unit. While this number itself has no significant meaning, it is interresting to modify the input parameters to see how the output is modified. If you are interrested in negligible absorption or perfect phasematching calculations, you can call amplitude_PMfree.m, amplitude_ABSfree.m and amplitude_ABSPM_free.m instead of amplitude.m.


#### Code tree structure

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

	------------phase					%calculate phase matching
	------------------refractive				%compute refractive index
	------------------Igauss				%intensity gaussian distribution
	------------------Press					%pressure distribution
	------------------ioniz					%ionization fraction

	------amplitude						%calculate final harmonic amplitude
	------------absorb					%absorption by the gas


## Main Results

Since there is more than 10 parameters involved in the HHG process that can be modified experimentally, and because the involved physical phenomenon are highly nonlinear, it is not possible to give an overview of all dependencies at the same time. This section only showcases some of the results obtained with our simulation. For details regarding the theory, implementation choices and the obtained results, you can refer to the pdf report.

<img src="https://raw.githubusercontent.com/Aurelien-Pelissier/High-Harmonic-Generation/master/report/results.png" width=900>

- On the left side we have the pressure dependance on the harmonic output for different harmonic order. 

- On the right side we plot the intensity dependance for different pressure


These two graph highlight how complex the HHG process can be. We have studied the influence of most of the input parameters, and the table below summarize how to maximize the harmonic output for each parameters.



| Parameter | How to optimize |
| --- | --- |
| Peak intensity Io | For any configuration, the optimum peak intensity is roughly 7W/cm2, which correspond to the intensity when the increasing of the dipole response does not compensate the destructive effects of the bad phasematching anymore (due to ions). Changing the gas, the beam radius or the laser frequency would modify this value. |
| Pressure P | Increasing the pressure indefinitely does not work because it degrades the phasematching, and also increase the absorption effects. The optimum pressure varies between each configurations (500mbar−1000mbar), and it is most of the time limited by the absorption rather than phasematching. |
| Interaction length lp | For our range of optimal pressure, the best interaction length is in the range of 50 μm - 120 μm which is smaller than what we can do experimentally.|
| Nozzle position znoz | For low interaction length, znoz = 0mm is optimum, but as we increase lp, moving the nozzle away from the focus may actually improve the harmonic output.|
| Time t | For high repetition rate systems, the highest harmonic power is always produced close to t = 0 fs.|




## Reports
For more details about the methods used in the simulation, or more detailed results, you can check the report and the research poster, available in the `report/` folder.
