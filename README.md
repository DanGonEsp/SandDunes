# Prerequisites:

	1. Kitware Cmake
	2. C++ compiler
	3. Git

# Installation:

To install the Daniel's version of the simulator, run the following commands in a folder where you want to have UG4:

	1.	mkdir UG4; cd UG4																		## Creating the folder UG4; it will contain the entire installation
	2. 	git clone https://github.com/UG4/ughub													## Cloning ughub, the package management system for ug4
	3.	mkdir ug4; cd ug4																		## Creating the folder ug4; it contains the ug4 library
	4.	../ughub/ughub init																		## Initializing the ug4 directory
	5. 	git clone -b ConvergenceRate https://github.com/DanGonEsp/ugcore.git					## Clonning the Daniel's repository
	6.	../ughub/ughub install ugcore --ignore	
	7.	mkdir plugins; cd plugins										  		
	8.	git clone -b MultiphaseFlowLinkers https://github.com/DanGonEsp/plugin_NavierStokes.git NavierStokes
	9.	git clone -b master https://github.com/UG4/plugin_Limex.git Limex  
	10.	cd ../..
	11.	mkdir build; cd build																	## Creating the temporary build folder
	12.	cmake ../ug4 -DDIM=2 -DCPU="1;4" -DNavierStokes=ON -DLimex=ON							## Initializing the directory and setting the build arguments
																								## In case of cmake incompatibility run with additional argument:
																								## -DCMAKE_POLICY_VERSION_MINIMUM=3.5
	12.	make																					## Compiling and linking the ug4 code (consider "make -j")
			
	Installation of the ug4 library and the NavierStokes plugin completed
	
	13.	cd ..																					## Returning to the UG4 directory
	14.	git clone https://gitlab.com/DanGon/sanddunes.git
	

# Running code

To perform one simulation of the sand dune, do:

	1. Go to UG4/sanddunes
	2. ../ug4/bin/ugshell -ex StabilizationStat_relative_vel_multiphase.lua -file_name Test -numRefs 1 -numPreRefs 0 -numTimeSteps 100 -DTmax 10 -Stokes true -boolAveDiff true -boolRelativeVel true


Additional parameters for the *StabilizationStat_relative_vel_multiphase.lua* script:


	* -AbsDefect							Minimum Numerical absolute defect 
	* -RedDefect							Reduction defect
	* -file_name							Folder name where solution will be saved
	* -elem_type			[tri|quad] 		
	* -numRefs				[0,1,2,3,...]	number of grid refinements
	* -numPreRefs							number of pre-refinements 
	* -numTimeSteps
	* -output				[1,2,3,...]		output files every n steps
	* -DTmax				[seconds]		Max time step
	* -DTmin				[seconds]		Min time step
	* -dt					[seconds] 		time step
	* -modifyDT				[true|false]	Adaptive time step

	* -boolRelativeVel.		[true|false]	Activate and deactivate the relative velocity
	* -boolAveDiff			[true|false]	Activate and deactivate the Erosion diffusion effect
	* -inflow				[0 - 10,m/s]	Inflow velocity (logarithmic profile)
	* -Height				[0-2,meter]		Dune Height
	* -Stokes								For solving stokes equations

	* -visc_a				[1.48e-05]		"air kinematic viscosity"
	* -rho_a				[1.2]			"Air Density"
	* -rho_s				[2500]			"Sand Density"
	* -diameter								"Sand Particle Diameter"
	* -c_init				[0.62 - 0.63]	"initial volume fraction"

The code will automatically create a folder where some output files for visualization and integral.txt are saved. In the integral.txt file  the total volume of sand of every domain is printed for different time steps.
