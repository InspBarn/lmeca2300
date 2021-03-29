##lmeca2300/project

This folder aims to gather the project of lmeca2300 course at UCLouvain.
It is composed of several directories including ./pylib/ and ./src/ where you'll find the python and de c library, respectively.
There are several source codes about the Cahn Hilliard algorithm.
Next, the directory ./deps/ gather the dependances for the project, which is the bov library (for graphism) exclusively.
Finally, the ./build/ directory is the one where you must be when you compile and run the project. Please follow the commands down here :

  cd project/build/
  cmake ..
  make && ./project

To make the simulation do what you want, you may open the project/src/main.c file where we set the parameters of the simulation :
N : the amount of points.
