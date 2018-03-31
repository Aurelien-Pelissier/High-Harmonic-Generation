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

	------------detectethan  				%check if the ionization(multiple pulse) file already exist
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
