mpic++ -O3 -o mpi_ising.x mpi_metropolis.cpp
mpirun -np 4 mpi_ising.x mpi_Lattice 40 100000 2.0 2.4 0.05
#mpirun -np 4 mpi_ising.x mpi_Lattice 60 100000 2.0 2.4 0.05
#mpirun -np 4 mpi_ising.x mpi_Lattice 80 100000 2.0 2.4 0.05
#mpirun -np 4 mpi_ising.x mpi_Lattice 100 100000 2.0 2.4 0.05



#python3 xyplot.py eLattice2 0 1 Temp[kT] E/J

#echo "L=2, MC=1000, T=1.0, step=0.1"
#cat Lattice2

