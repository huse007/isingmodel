/* SPØRSMÅL 
 * 1. c) Plot energy vs. #MC cycles? En og en kjøring eller printe ut energy ved hver cycle
 * 2. d) Histogrammet for T=1.0 blir rart (det er rikitg)
 * 3. d) Stemmer verdiene på x aksen -> peak på -700 ?? hvorfor ikke 0 (det er riktig)

/* Run as
  ./metropolis.x ofile npins nMCcycles init_temp final_temp tempstep
   Example: 
   ./metropolis.x Lattice 100 10000000 2.1 2.4 0.01
   Compile: 
   c++ -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x IsingModel.cpp -larmadillo (-Rpass only with clang)
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include <mpi.h>
using namespace  std;
using namespace arma;

ofstream ofile;

/* inline function (compiler optimalization) for PBC */
inline int periodicBoundary(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Init spin_matrix, energy and magnetization
void initLattice(int, mat &, double&, double&);
void metropolisSampling(int, int, double, vec &);  
void writeToFile(int, int, double, vec);

int main(int argc, char* argv[])
{
  // Init MPI
  int myRank, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  
  string filename;
  int nspins, mc_cycles;
  double t_init, t_final, t_step;
  if (myRank == 0 && argc <= 5) {
    cout << "Invalid command: " << argv[0] << 
      " outputfile spins MC_cycles t_init t_final t_step" << endl;
    exit(1);
  }
  if (myRank == 0 && argc > 1) {
    filename=argv[1];
    nspins = atoi(argv[2]);
    mc_cycles = atoi(argv[3]);    
    t_init = atof(argv[4]);
    t_final = atof(argv[5]);
    t_step = atof(argv[6]);
  }
  
  MPI_Bcast (&nspins, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&t_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&t_final, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&t_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  


  
  /* init outfile */
  if(myRank == 0) {
    string fileout = filename;
    string argument = to_string(nspins);
    fileout.append(argument);
    ofile.open(fileout);
  }
  
  /* loop over temperature t*/
  for (double t = t_init; t <= t_final; t+=t_step){
    vec ExpectationValues = zeros<mat>(5);
    /* Monte Carlo (put results in ExpectationValues)*/
    metropolisSampling(nspins, mc_cycles, t, ExpectationValues);
    /* write to file */
    for(int i = 0; i < 5; i++)
      MPI_Allreduce(&ExpectationValues[i], &ExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if(myRank == 0) writeToFile(nspins, mc_cycles*numProcs, t, ExpectationValues);
  }
  if(myRank == 0) ofile.close();
  MPI_Finalize();
  return 0;
}

// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void metropolisSampling(int nspins, int mc_cycles, double Temperature, vec &ExpectationValues)
{
  /* init the seed
     call the Mersienne algo 
     setup uniform dist [0,1] */
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  mat spin_matrix = zeros<mat>(nspins,nspins);   // init the lattice spin values
  double Energy = 0.;// init energy 
  double MagneticMoment = 0.;//init magnetization 
  initLattice(nspins, spin_matrix, Energy, MagneticMoment);  /* initialize array for expectation values */
  vec EnergyDifference = zeros<mat>(17); /* setup array for possible energy changes */
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  /* start Monte Carlo cycles */
  for (int cycles = 1; cycles <= mc_cycles; cycles++){
    for(int x =0; x < nspins; x++) {
      for (int y= 0; y < nspins; y++){
	int ix = (int) (RandomNumberGenerator(gen)*(double)nspins);
	int iy = (int) (RandomNumberGenerator(gen)*(double)nspins);
	int deltaE =  2*spin_matrix(ix,iy)*
	  (spin_matrix(ix,periodicBoundary(iy,nspins,-1))+
	   spin_matrix(periodicBoundary(ix,nspins,-1),iy) +
	   spin_matrix(ix,periodicBoundary(iy,nspins,1)) +
	   spin_matrix(periodicBoundary(ix,nspins,1),iy));
	if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
	  spin_matrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
	  MagneticMoment += (double) 2*spin_matrix(ix,iy);
	  Energy += (double) deltaE;
	  //cout<<Energy<<endl; //exercise d) Verdier til P(E) histogram 
	}
      }
    }
    /* update expectation values */
    ExpectationValues(0) += Energy;
    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;    
    ExpectationValues(3) += MagneticMoment*MagneticMoment; 
    ExpectationValues(4) += fabs(MagneticMoment);
  }
}

void initLattice(int nspins, mat &spin_matrix,  double& Energy, double& MagneticMoment)
{
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
  // setup spin matrix and initial magnetization
  for(int x =0; x < nspins; x++) {
    for (int y= 0; y < nspins; y++){
      spin_matrix(x,y) = 1.0; // spin orientation for the ground state
      //(RandomNumberGenerator(gen)>0.0) ? spin_matrix(x,y) = 1.0 : spin_matrix(x,y) = -1.0; //random start spin config
      MagneticMoment +=  (double) spin_matrix(x,y);
    }
  }
  spin_matrix.print();
  // setup initial energy
  for(int x =0; x < nspins; x++) {
    for (int y= 0; y < nspins; y++){
      Energy -=  (double) spin_matrix(x,y)*
	(spin_matrix(periodicBoundary(x,nspins,-1),y) +
	 spin_matrix(x,periodicBoundary(y,nspins,-1)));
    }
  }
}

void writeToFile(int nspins, int mc_cycles, double temperature, vec ExpectationValues)
{
  double norm = 1.0/((double) (mc_cycles));  // divided by  number of cycles 
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm;
  // all expectation values are per spin, divide by 1/nspins/nspins
  double Evariance = (E2_ExpectationValues- E_ExpectationValues*E_ExpectationValues)/nspins/nspins;
  double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues)/nspins/nspins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << E_ExpectationValues/nspins/nspins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << M_ExpectationValues/nspins/nspins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/nspins/nspins << endl;
}

