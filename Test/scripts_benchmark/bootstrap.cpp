#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <cmath>
#include <limits>
#include <unistd.h>
#include <gsl/gsl_rng.h>                // gsl_rng
#include <gsl/gsl_randist.h>    // gaussian and other distribution
#include <gsl/gsl_math.h>

using namespace std;

///////  MAIN
int main(int argc, char*argv[])
{

	if (argc != 3+1)
        {
                cout<<"Inputs:"<<endl;
                cout<<"1 --> trajectory (only 1 column)"<<endl;
                cout<<"2 --> seed"<<endl;
                cout<<"3 --> repetitions"<<endl;
                return 0;
        }
	string trajectory(argv[1]);
	int seed = atoi(argv[2]);
	int repetitions = atoi(argv[3]);

	cout << setprecision(10);
	//Initialize random number generator for white noise
	gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus2);
        gsl_rng_set(rng, static_cast<unsigned int>(seed));
	
	ifstream leggi;
	if (access(trajectory.c_str(), F_OK ) != -1 )
        {
                leggi.open(trajectory.c_str(), ifstream::in);
        }
	else
	{
		cerr<<"Trajectory file not found! Exiting..."<<endl;
		return 0;
	}

	//Import data
	vector <double> lista;
	double comodo;
	while (!leggi.eof())
	{
		leggi >> comodo;
		if (!leggi.eof())
			lista.push_back(comodo);
	}

	double sum = 0, sumq = 0, run_ave, run_aveq, sigma;
	double ave;
	int npoints = lista.size();
	for (int i=1; i<=repetitions; i++)
	{
		ave = 0;
		for (int j=0; j<npoints; j++)
			ave += lista[gsl_rng_uniform_int(rng, npoints)];
		ave *= 1./npoints;
		sum += ave;
		sumq += ave*ave;
		if (i>1)
		{
			run_ave = sum/i;
			run_aveq = sumq/i;
			sigma = sqrt(run_aveq-run_ave*run_ave);
			cout<<i<<" "<<run_ave<<" "<<sigma<<endl;
		}
	}

	return 1;
}
