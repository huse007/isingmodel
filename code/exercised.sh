#Change the metropolis.cpp -> Print out all the accepted energies
c++ -O3 -o Ising.x metropolis.cpp
./Ising.x Lattice 20 100000 1.0 1.0 0.1 > histogramdata_t1.txt
python3 histogram.py histogramdata_t1.txt


