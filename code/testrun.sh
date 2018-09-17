#c++ -O3 -std=c++11 -o timing_Ising.x timing_metropolis.cpp                                              
./Ising.x test2Lattice 100 10000 2.20 2.20 0.001
#./timing_Ising.x test2Lattice 60 100 2.20 2.20 0.001
#./timing_Ising.x test2Lattice 80 100 2.20 2.20 0.001
#./timing_Ising.x test2Lattice 100 100 2.20 2.20 0.001

#/usr/lib64/openmpi/bin/mpic++ -std=c++11 -O3 -o timing_mpi_ising.x timing_mpi_metropolis.cpp
#/usr/lib64/openmpi/bin/mpirun -n 4 timing_mpi_ising.x testL 100 100000 2.200 2.200 0.001
#/usr/lib64/openmpi/bin/mpirun -np 4 timing_mpi_ising.x testL 60 100 2.200 2.200 0.001
#/usr/lib64/openmpi/bin/mpirun -np 4 timing_mpi_ising.x testL 80 100 2.200 2.200 0.001
#/usr/lib64/openmpi/bin/mpirun -np 4 timing_mpi_ising.x testL 100 100 2.200 2.200 0.001

##/usr/lib64/openmpi/bin/mpirun -n 8 timing_mpi_ising.x testL 100 100000 2.200 2.200 0.001
#/usr/lib64/openmpi/bin/mpirun -np 8 timing_mpi_ising.x testL 60 100 2.200 2.200 0.001
#/usr/lib64/openmpi/bin/mpirun -np 8 timing_mpi_ising.x testL 80 100 2.200 2.200 0.001
#/usr/lib64/openmpi/bin/mpirun -np 8 timing_mpi_ising.x testL 100 100 2.200 2.200 0.001



