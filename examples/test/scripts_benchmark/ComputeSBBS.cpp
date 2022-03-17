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

////// GLOBAL VARIABLES
int natoms, lseq;
double timestep;
vector <double> coordinates[3];
vector <int> tipologia;
ifstream leggi;

////// DECLARATION OF FUNCTIONS
//DUMP READER
int ReadDump(double, double, int);

//COMPUTATIONS
void cnt(double *, int);
void bpvec(double *, int, int);

//VECTORIAL OPERATIONS
void rotate(double *, double *, double, double*);
void cross_product(double *, double *, double *);
double dot_product(double *, double *, int);
void versor(double *, int, double*);
double angle_vec(double *, double *);
double angle_vec_sign(double *, double *, double*);
void cp_vec(double *, double *, int);
void diff_vec(double *, double *, double *, int);
void print_vec(double *vettore, int vecsize);

///////  MAIN
int main(int argc, char*argv[])
{
	if (argc != 2)
	{
		cout<<"Inputs:"<<endl;
		cout<<"1 --> trajectory"<<endl;
		cout<<"Output: SBBS dihedral (deg)"<<endl;
		return 0;
	}
	string trajectory(argv[1]);

	cout << setprecision(10);

	//Import first frame
	int flag;
	if (access(trajectory.c_str(), F_OK ) != -1 )
        {
                leggi.open(trajectory.c_str(), ifstream::in);
                flag = ReadDump(-1, 1e10, 1);
        }
	else
	{
		cout<<"Trajectory file not found! Exiting..."<<endl;
		return 0;
	}

	int iA, iB, iC, iD;
	double rAB[3], rBC[3], rCD[3], pvec1[3], pvec2[3], phi, segno;

	//output first line
	cout<<"t";
	for (int iseq = 1; iseq <= lseq; iseq++)
        {
		cout<<" "<<iseq;
	}
	cout<<endl;

	//Analyze trajectory
	while (flag > 0)
        {
		cout<<timestep;
                for (int iseq = 1; iseq < lseq; iseq++)
                {
			//compute involved beads
                        iA = 3*iseq - 2;
			iB = iA + 1;
			iC = natoms + 2 - iB; 
			iD = iC - 1;

			//compute displacement vectors
			for (int k1=0; k1<3; k1++)
			{
				rAB[k1] = coordinates[k1][iB] - coordinates[k1][iA];
				rBC[k1] = coordinates[k1][iC] - coordinates[k1][iB];
				rCD[k1] = coordinates[k1][iD] - coordinates[k1][iC];
			}

			//compute cross products and normalize
			cross_product(rAB, rBC, pvec1);
			cross_product(rBC, rCD, pvec2);
			versor(pvec1, 3, pvec1);
			versor(pvec2, 3, pvec2);

			//compute phi
			phi = 180/M_PI*angle_vec(pvec1, pvec2);
			segno = dot_product(pvec1, rCD, 3);
			if (segno < 0)
				phi *= -1;

			//output
			cout<<" "<<phi;
                }
		cout<<endl;

		//read next frame
                flag = ReadDump(-1, 1e10, 0);
        }
        leggi.close();

	return 1;
}



///////////////// FUNCTIONS 

int ReadDump(double t00, double t11, int build)
{
	string comstr;
        getline(leggi, comstr);
	double comodo;
	int indice, tipo;
        leggi>>comodo;
        if (!leggi.eof())
        {
                if (comodo >= t00 && comodo <= t11)
                {
                        timestep = comodo;
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        //number of atoms
                        leggi>>natoms;
			if (build == 1)
			{
				lseq = (natoms+2)/6;
                                tipologia.resize(natoms+1);
                                for (int k1=0; k1<3; k1++)
                                {
                                        coordinates[k1].resize(natoms+1);
                                }
			}
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        //size of box
			for (int ii = 1; ii<=6; ii++)
                                leggi>>comodo;
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        //coordinates
                        for (int i=1; i<=natoms; i++)
                        {
                                leggi>>indice;
                                leggi>>tipo;
                                tipologia[indice] = tipo;
				leggi>>comodo;
                                coordinates[0][indice] = comodo;
                                leggi>>comodo;
                                coordinates[1][indice] = comodo;
                                leggi>>comodo;
                                coordinates[2][indice] = comodo;

                        }
                        getline(leggi, comstr);
                }
                else
                {
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        //number of atoms
                        leggi>>natoms;
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        getline(leggi, comstr);
                        for (int i=1; i<=natoms; i++)
                                getline(leggi, comstr);
                }
                return 1;
        }
        else
                return 0;
}

void cnt(double *vettore, int nuc)
{
	int iA = 3*nuc - 2;
	int iB = natoms - iA;
	for (int k1=0; k1<3; k1++)
	    	vettore[k1] = 0.5*(coordinates[k1][iA] + coordinates[k1][iB]);
	return;
}

void bpvec(double *vettore, int nuc)
{
	int iA = 3*nuc - 2;
	int iB = natoms - iA;
	for (int k1=0; k1<3; k1++)
	    	vettore[k1] = coordinates[k1][iB] - coordinates[k1][iA];
	return;
}

void rotate(double *vettore, double *axis, double phi, double *vettoredest)
{
    //Use Rodrigues' formula
    double pscal = dot_product(axis, vettore, 3);
    double pvec[3];
    cross_product(axis, vettore, pvec);
    for (int k1=0; k1<3; k1++)
	    vettoredest[k1] = vettore[k1]*cos(phi) + pvec[k1]*sin(phi) + axis[k1]*pscal*(1-cos(phi));
    return;
}

void cross_product(double *vettore1, double *vettore2, double *vettoredest)
{
	vettoredest[0] = vettore1[1]*vettore2[2] - vettore1[2]*vettore2[1];
	vettoredest[1] = vettore1[2]*vettore2[0] - vettore1[0]*vettore2[2];
	vettoredest[2] = vettore1[0]*vettore2[1] - vettore1[1]*vettore2[0];
	return;
}

double dot_product(double *vettore1, double *vettore2, int sizevec)
{
	double ss = 0;
	for (int k1=0; k1<sizevec; k1++)
		ss += vettore1[k1]*vettore2[k1];
	return ss;
}

void versor(double *vettore, int sizevec, double *vettoredest)
{
	double nn = sqrt(dot_product(vettore, vettore, sizevec));
	for (int k1=0; k1<sizevec; k1++)
		vettoredest[k1] *= 1./nn;
}

double angle_vec(double *v1, double *v2)
{
    versor(v1,3,v1);
    versor(v2,3,v2);
    double sp = dot_product(v1,v2,3);
    if (sp > 1.0)
	    sp = 1.0;
    if (sp < -1.0)
	    sp = -1.0;
    return acos(sp);
}

double angle_vec_sign(double *v1, double *v2, double *asse)
{
	versor(v1,3,v1);
	versor(v2,3,v2);
	double theta = angle_vec(v1, v2);
	double cp[3];
	cross_product(v1, v2, cp);
	double checksegno = dot_product(cp, asse, 3);
	if (checksegno <0)
		theta *= -1;
	return theta;
}

void cp_vec(double *vettore_origine, double *vettore_destinazione, int sizevec)
{
	for (int k1=0; k1<sizevec; k1++)
		vettore_destinazione[k1] = vettore_origine[k1];
	return;
}

//destvec = vecA - vecB
void diff_vec(double *vecA, double *vecB, double *destvec, int sizevec)
{
	for (int k1=0; k1<sizevec; k1++)
		destvec[k1] = vecA[k1] - vecB[k1];
	return;
}

void print_vec(double *vettore, int vecsize)
{
        for (int i=0; i<vecsize; i++)
                cout<<vettore[i]<<" ";
        cout<<endl;
}
