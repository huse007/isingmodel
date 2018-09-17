# The Ising Model

The algorithm is explained in rapport/ising_model.pdf

## Ubuntu and Red Hat
Requirements:
```
libarmadillo-dev
```
Compile:
```
c++ -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x IsingModel.cpp -larmadillo (-Rpass only with clang)
```
Test script:
```
source testrun.sh
```
Run the code with the following arguments:
```
./ofile.x outputfile npins nMCcycles init_temp final_temp tempstep
```

Example:
```
./Ising.x Lattice 100 10000000 2.1 2.4 0.001
```

For big lattices or a high number of monte carlo cycles,
use MPI. Example:
```
mpic++ -std=c++11 -O3 -o mpi_metropolis.x metropolis.cpp
mpirun -n 8 mpi_metropolis.x outfile 100 100000 2.200 2.200 0.001
```