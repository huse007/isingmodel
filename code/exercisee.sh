mpic++ -O3 mpi_metropolis.cpp -o mpi_ising.x
#mpirun -np 4 mpi_metropolis eelodie_mpi 10 10 2.0 2.06 0.05
mpirun -np 4 mpi_ising.x mpi_Lattice 40 1000000 2.0 2.4 0.05
#mpirun -np 4 mpi_ising.x mpi_Lattice 60 100000 2.0 2.4 0.05
#mpirun -np 4 mpi_ising.x mpi_Lattice 80 100000 2.0 2.4 0.05
#mpirun -np 4 mpi_ising.x mpi_Lattice 100 100000 2.0 2.4 0.05



#python3 xyplot.py eLattice2 0 1 Temp[kT] E/J

#echo "L=2, MC=1000, T=1.0, step=0.1"
#cat Lattice2

