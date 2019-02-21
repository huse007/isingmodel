# The Ising Model

A mathematical model of ferromagnetism in statistical mechanics. The model consists of discrete variables that represent magnetic dipole moments of atomic spins that can be in one of two states (+1 or −1). The spins are arranged in a lattice, allowing each spin to interact with its neighbors. The model allows the identification of phase transitions, as a simplified model of reality. 

The algorithm and results are explained in rapport/ising_model.pdf

Note: Big runs will take several days, so it is recommended to use multiple cores or computers. 

![alt text](https://raw.githubusercontent.com/huse007/Isingmodel/master/rapport/ising.png)

## Ubuntu (apt) and Red Hat (rpm)
Requires the C++ linear algebra library:
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
./Ising.x outputfile npins nMCcycles init_temp final_temp tempstep
```

Example:
```
./Ising.x Lattice 100 10000000 2.1 2.4 0.001
```
## Parallell Computing
For big lattices or a high number of monte carlo cycles,
use MPI. Example:
```
mpic++ -std=c++11 -O3 -o mpi_metropolis.x metropolis.cpp
mpirun -n 8 mpi_metropolis.x outfile 100 100000 2.200 2.200 0.001
```
