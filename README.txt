* example.R 
	Runs a MCMC estimation for a simulated pattern (i.e. does it all).

* src/helpers.cpp
	Functions used for calculations in the code, such as the kernel integral over the window
	and pmfs etc. 

* src/par_updates.cpp
	Functions for the MH-updates of the parameters.

* src/poit_updates.cpp
	Functions for updating the locations of the parents.

* R/estgtpr.R
	Wrapper to simulate a pattern and run the mcmc.

* R/estgtp.R
	The actual execution of the MCMC estimation.
	
* R/rgtp.R	
  Simulation of generalized Thomas process.