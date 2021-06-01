# lmeca2300/project

This folder aims to gather the project of lmeca2300 (Advanced Numerical Method) course at UCLouvain for academic year 2020-21. We study along this year the spectral method in aim to derive function. As a project, we were asked to present an implementation of the Cahn Hilliard integrator. As an objective of the project, we were also asked to prove that spectral methods are more interesting than classical spatial differentiation.

The visiter will find in the pylib folder the python source code. The two files cahn_hilliard.py and cahn_hilliard_sd.py are respectively files for the spatial and the spectral derivation. Launching the latter will show you the nice animation we have been working on. Other files in the pylib folder are used for plotting graph we showed in the report.

The visiter will find on the other hand the c source code in src folder. For launching this part of the project, you can find more information in the compilation section.

## Dependencies

Depencdencies for visualize the project are in the deps/ folder. Yet you need to have an accessible fft algorithm. We use for that purpose the FFTW3 library for which documentation can be found on http://fftw.org. Please follow the installation rules if you did not install it before.

## Compile the C project

It is composed of several directories including ./pylib/ and ./src/ where you'll find the python and de c library, respectively.
There are several source codes about the Cahn Hilliard algorithm.
Next, the directory ./deps/ gather the dependances for the project, which is the bov library (for graphism) exclusively.
Finally, the ./build/ directory is the one where you must be when you compile and run the project. Please follow the commands down here :

	mkdir project/build/
	cd project/build/
	cmake ..
	make && ./project

Once it has been built, but you have commit changes in the source code, you may skip the first two commands and only launch :

	make && ./project

If you didn't commit any changes, you may not compile the project again and skip the first part :

	./project

The simulation will be displayed and computation times will be printed in the CLI. For stopping the simulation, simply press on ESC.

## Change simulation parameters

To make the simulation do what you want, you may open the project/src/main.c file where you may give the amount of discretization points N you want. Since we face a 2d-problem, we work on a N x N points grid.

Concerning the python part of project, the reader can modify whether the discretization and the temporal integrator in the respecting files.

## Versions

There are two possible versions. Although they are working the same way, openGL is much more efficient for plot rendering. Nevertheless, we came accross several problems and it is working only for N=256. No suitable reasons has been found for. The simpler version is using the siplified BOV library which itself depends on openGL. It is much less faster yet it is working.
