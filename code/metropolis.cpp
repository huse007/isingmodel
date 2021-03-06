/* Run as
  ./ofile.x outputfile npins nMCcycles init_temp final_temp tempstep
   Example: 
   ./Ising.x Lattice 100 10000000 2.1 2.4 0.001
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
  /* timing */
  clock_t start,stop;
  start = clock();
  string filename;
  int nspins, mc_cycles;
  double t_init, t_final, t_step;
  /* invalid command */
  if (argc <= 5) {
    cout << "Invalid command: " << argv[0] << 
      " outputfile spins MC_cycles t_init t_final t_step" << endl;
    exit(1);
  }
  /* place parameters (does not throw exep if ex. argc == 4 */
  if (argc > 1) {
    filename=argv[1];
    nspins = atoi(argv[2]);
    mc_cycles = atoi(argv[3]);    
    t_init = atof(argv[4]);
    t_final = atof(argv[5]);
    t_step = atof(argv[6]);
  }
  /* init outfile */
  string fileout = filename;
  string argument = to_string(nspins);
  fileout.append(argument);
  ofile.open(fileout);
  
  /* loop over temperature t */
  for (double t = t_init; t <= t_final; t+=t_step){
    vec ExpectationValues = zeros<mat>(5);
    /* Monte Carlo (put results in ExpectationValues)*/
    metropolisSampling(nspins, mc_cycles, t, ExpectationValues);
    /* write to file */
    writeToFile(nspins, mc_cycles, t, ExpectationValues);
  }
  /* end timing */
  stop = clock();
  cout<<"Algorithm time: "<<((float)(stop-start)/CLOCKS_PER_SEC)<<"s"<<endl;
  ofile.close();  
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
  /* construct a spin_matrix */
  mat spin_matrix = zeros<mat>(nspins,nspins);   
  double Energy = 0.; 
  double MagneticMoment = 0.;
  /* initialize array for expectation values */
  initLattice(nspins, spin_matrix, Energy, MagneticMoment);  
  /* setup array for possible energy changes */
  vec EnergyDifference = zeros<mat>(17); 
  /* calculate exp(-beta dE) before MC loop*/
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  /* start Monte Carlo cycles */
  for (int cycles = 1; cycles <= mc_cycles; cycles++){
    int count= 0;
    /* loop lattice */
    for(int x =0; x < nspins; x++) {
      for (int y= 0; y < nspins; y++){
	int ix = (int) (RandomNumberGenerator(gen)*(double)nspins);
	int iy = (int) (RandomNumberGenerator(gen)*(double)nspins);
	int deltaE =  2*spin_matrix(ix,iy)*
	  (spin_matrix(ix,periodicBoundary(iy,nspins,-1))+
	   spin_matrix(periodicBoundary(ix,nspins,-1),iy) +
	   spin_matrix(ix,periodicBoundary(iy,nspins,1)) +
	   spin_matrix(periodicBoundary(ix,nspins,1),iy));
	/* accept config if random<=w */
	if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
	  /* flip one spin */
	  spin_matrix(ix,iy) *= -1.0;
	  MagneticMoment += (double) 2*spin_matrix(ix,iy);
	  Energy += (double) deltaE;
	  count++;
	  /* Print values for P(E) histogram in exercise d) */
	  //cout<<Energy<<endl; 
	}
      }
    }
    /* update expectation values */
    ExpectationValues(0) += Energy;
    ExpectationValues(1) += Energy*Energy;
    ExpectationValues(2) += MagneticMoment;    
    ExpectationValues(3) += MagneticMoment*MagneticMoment; 
    ExpectationValues(4) += fabs(MagneticMoment);
    /* As function of MC cycles in exercise c) */
    /*if(cycles%100 == 0) {
      double E = ExpectationValues(0) / (double) cycles/nspins/nspins;
      double M = ExpectationValues(4) / (double) cycles/nspins/nspins;
      cout<<cycles<< " " << E << " "<<M <<" "<<count<<endl;
      }*/
  }
}

void initLattice(int nspins, mat &spin_matrix,  double& Energy, double& MagneticMoment)
{
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(-1.0,1.0);
  /* setup spin matrix and initial magnetization */
  for(int x =0; x < nspins; x++) {
    for (int y= 0; y < nspins; y++){
      /* Select spin configuration */
      //spin_matrix(x,y) = 1.0; // spin orientation for the ground state
      (RandomNumberGenerator(gen)>0.0) ? spin_matrix(x,y) = 1.0 : spin_matrix(x,y) = -1.0; //random start spin config
      MagneticMoment +=  (double) spin_matrix(x,y);
    }
  }
  //  spin_matrix.print();
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

