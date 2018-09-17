# The Ising Model

The algorithm is explained in rapport/ising_model.pdf

## Red Hat and Ubuntu
Compile:
```
c++ -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x IsingModel.cpp -larmadillo (-Rpass only with clang)
```
Run the code with the following arguments:
```
./ofile.x outputfile npins nMCcycles init_temp final_temp tempstep
```

Example:
```
./Ising.x Lattice 100 10000000 2.1 2.4 0.001
```
