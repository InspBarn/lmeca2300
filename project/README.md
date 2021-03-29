# lmeca2300/project

This folder aims to gather the project of lmeca2300 (Advanced Numerical Method) course at UCLouvain for academic year 2020-21. We study along this year the spectral method in aim to derive function. As a project, we were asked to present an implementation of the Cahn Hilliard integrator. As an objective of the project, we were also asked to prove that spectral methods are more interesting than classical spatial differentiation.

## Compile the project

It is composed of several directories including ./pylib/ and ./src/ where you'll find the python and de c library, respectively.
There are several source codes about the Cahn Hilliard algorithm.
Next, the directory ./deps/ gather the dependances for the project, which is the bov library (for graphism) exclusively.
Finally, the ./build/ directory is the one where you must be when you compile and run the project. Please follow the commands down here :

	cd project/build/
	cmake ..
	make && ./project

Once it has been built, but you have commit changes in the source code, you may skip the first two commands and only launch :

	make && ./project

If you didn't commit any changes, you may not compile the project again and skip the first part :

	./project

## Change simulation parameters

To make the simulation do what you want, you may open the project/src/main.c file where we set the parameters of the simulation :