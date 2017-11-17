#c++ -O3 -o Ising.x metropolis.cpp
./Ising.x Lattice 2 10 1.0 1.0 0.1
echo "MC cycles: 10"
cat Lattice2
./Ising.x Lattice 2 100 1.0 1.0 0.1
echo "MC cycles: 100"
cat Lattice2
./Ising.x Lattice 2 1000 1.0 1.0 0.1
echo "MC cycles: 1000"
cat Lattice2
./Ising.x Lattice 2 10000 1.0 1.0 0.1
echo "MC cycles: 10000"
cat Lattice2
./Ising.x Lattice 2 100000 1.0 1.0 0.1
echo "MC cycles: 100000"
cat Lattice2
./Ising.x Lattice 2 1000000 1.0 1.0 0.1
echo "MC cycles: 1000000"
cat Lattice2
./Ising.x Lattice 2 10000000 1.0 1.0 0.1
echo "MC cycles: 10000000"
cat Lattice2
