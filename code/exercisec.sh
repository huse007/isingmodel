#script exercise c)
c++ -O3 -o Ising.x metropolis.cpp
./Ising.x Lattice1a 20 1000000 1.0 1.0 0.1 > newexercisec_allup_t1p0.txt
./Ising.x Lattice2a 20 1000000 2.4 2.4 0.1 > newexercisec_allup_t2p4.txt
#./Ising.x Lattice1r 20 1000000 1.0 1.0 0.1 > newexercisec_random_t1p0.txt
#./Ising.x Lattice2r 20 1000000 2.4 2.4 0.1 > newexercisec_random_t2p4.txt



