# High-Harmonic-Generation
<img align="left" src="https://raw.githubusercontent.com/Aurelien-Pelissier/High-Harmonic-Generation/master/report/HHG.png" width=200>
High harmonic generation (HHG) refers to the process of creating vacuum (VUV) or extreme (XUV) ultraviolet light through a nonlinear interaction of an intense laser field with a gas target. This repository contains the source code to calculate the high harmonics amplitude accounting for all involved physical phenomenons (Quantum mechanic, Supersonic flow physics, Plasma physics and Nonlinear optics). More specifically, we consider supersonic gas flow at the nozzle outlet, ions dynamic in the plasma, quantum atomic response of freed electrons, phase matching and absorption.




## Involved parameters

%General parameters:
	q =21;                  %Harmonic order (1:51)
	t = 0;                  %Time [s] (gaussian pulse centered in t = 0)
	alpha = 2e-14;			    %Phase coefficient (taken at 2e-14 cm2/W for the whole study)
	Te = 3;	            	  %Freed electron temperature [eV] (taken as 3eV for the whole study)


%Laser parameters:
	I0 = 6e13;              %Peak intensity [W/cm2]
	tp = 130e-15;			      %Pulse length (FWHM)
	lambda1 = 1050e-9;		  %Fundamental wavelength [m]
	R0 = 19.6e-6;			      %Beam radius [m]
	f = 60e6; 	           	%Laser frequency [Hz]


%Gas parameters:
	V = 250;          		  %Gas velocity [m/s]
	P = 500;                %Peak pressure [mbar]
	lp = 150e-6; 			      %Interaction length (FWHM if gaussian)/ nozzle diameter [m]
	profile = 'gauss'; 	  	%Density profile ('squar' or 'gauss')
	znozzle = 0;	       	  %Nozzle position [m]
	gas = 'Kr';     		    %Gas ('Ar', 'Kr', 'Xe')
	xHe = 0;                %Helium fraction (0:1)



## Code structure

main

	detectdipole
		  dipole
			    potde				          %potential depth

	phasematching
		  detectetha  				      %check if file already exist
			    ethap 				        %ioniation fraction after one pulse
				      Ipot			        %ionization potentiel
				      Nt

		detectethan  				        %check if file already exist
			  etharz 				          %ioniation fraction after n pulse
				    ethabp 			        %between 2 pulses
					      solvepde 	      %for a fixed z
						        Ndens
						        mobility
						        Press

		phase
			  refractive
			  Igauss
			  Press
			  ioniz

	  amplitude
		    absorb





## Running the code









## Results




## Reports
