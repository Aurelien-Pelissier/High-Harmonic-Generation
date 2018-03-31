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


<p><span class="gauche">Texte aligné à gauche</span><span class="droite">Texte aligné à droite</span></p>



## Code structure

main

	detectdipole
		  dipole
			    potde		              %potential depth

	phasematching
		  detectetha  				      %check if file already exist
			    ethap 			      %ioniation fraction after one pulse
				      Ipot		      %ionization potentiel
				      Nt

		detectethan  				      %check if file already exist
			  etharz 			      %ioniation fraction after n pulse
				    ethabp 		      %between 2 pulses
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
