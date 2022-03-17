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
int natoms, lseq, nph, niter_cyl = 10000;
double timestep;
double dt_precision = 0.05, groove_width_shift = 5.8, groove_depth_shift = 3.5, diameter, haxis[3], centro[3], tmax[2];
gsl_rng *rng=gsl_rng_alloc(gsl_rng_taus2);
vector <double> coordinates[3], Phosphates[2][3], amat[2][3], bmat[2][3], cmat[2][3], dmat[2][3], t1t2mat[2][2], dtlist[2], grooves[4], hrise, htwist;
vector <int> tipologia;
ifstream leggi;

////// DECLARATION OF FUNCTIONS
//DUMP READER
int ReadDump(double, double, int);

//COMPUTATIONS
void cyl_info(int, int);
void cnt(double *, int);
void bpvec(double *, int, int);
double compute_htwist(int);
double compute_hrise(int);
void spline(int strand, int iph);
void interpolazione(double tglob, double *punto, int strand);
void groove();
void find_minimum_right(double *widdep, int iph, double *punto);
void find_minimum_left(double *widdep, int iph, double *punto);

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

	if (argc != 4)
	{
		cout<<"The command needs 3 inputs:"<<endl;
		cout<<"1 --> trajectory"<<endl;
		cout<<"2 --> seed for random number generation"<<endl;
		cout<<"3 --> prefix for output files (if folder does not exist, you need first to create it)"<<endl;
		return 0;
	}
	string trajectory(argv[1]);
	int seed = atoi(argv[2]);
	string prefix(argv[3]);

	cout << setprecision(10);
	//Initialize random number generator for white noise
        gsl_rng_set(rng, static_cast<unsigned int>(seed));
	
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

	//Open streams
	ofstream scrivi_diametro, scrivi_htwist, scrivi_hrise;
        ofstream scrivi_major_groove_width, scrivi_minor_groove_width, scrivi_major_groove_depth, scrivi_minor_groove_depth;
        char nome[1000];
        sprintf(nome,"%s_diameter", prefix.c_str());
        scrivi_diametro.open(nome);
        sprintf(nome,"%s_htwist", prefix.c_str());
        scrivi_htwist.open(nome);
        sprintf(nome,"%s_hrise", prefix.c_str());
        scrivi_hrise.open(nome);
        sprintf(nome,"%s_major_groove_width", prefix.c_str());
        scrivi_major_groove_width.open(nome);
        sprintf(nome,"%s_minor_groove_width", prefix.c_str());
        scrivi_minor_groove_width.open(nome);
        sprintf(nome,"%s_major_groove_depth", prefix.c_str());
        scrivi_major_groove_depth.open(nome);
        sprintf(nome,"%s_minor_groove_depth", prefix.c_str());
        scrivi_minor_groove_depth.open(nome);

	//Write labels in files
	scrivi_diametro<<"t diameter"<<endl;
        scrivi_hrise<<"t";
        scrivi_htwist<<"t";
        for (int iseq = 1; iseq < lseq; iseq++)
        {
                scrivi_hrise<<" "<<iseq<<"-"<<iseq+1;
                scrivi_htwist<<" "<<iseq<<"-"<<iseq+1;
        }
        scrivi_hrise<<endl;
        scrivi_htwist<<endl;
        scrivi_major_groove_width<<"t";
        scrivi_minor_groove_width<<"t";
        scrivi_major_groove_depth<<"t";
        scrivi_minor_groove_depth<<"t";
        for (int iseq = 1; iseq <= lseq; iseq++)
        {
                scrivi_major_groove_width<<" "<<iseq;
                scrivi_minor_groove_width<<" "<<iseq;
                scrivi_major_groove_depth<<" "<<iseq;
                scrivi_minor_groove_depth<<" "<<iseq;
        }
        scrivi_major_groove_width<<endl;
        scrivi_minor_groove_width<<endl;
        scrivi_major_groove_depth<<endl;
        scrivi_minor_groove_depth<<endl;

	//Analyze trajectory
	while (flag > 0)
	{
                cyl_info(1, lseq);
                scrivi_diametro<<timestep<<" "<<diameter<<endl;

		scrivi_hrise<<timestep;
                scrivi_htwist<<timestep;
                for (int iseq = 1; iseq < lseq; iseq++)
                {
                        scrivi_hrise<<" "<<0.1*compute_hrise(iseq);
                        scrivi_htwist<<" "<<180/M_PI*compute_htwist(iseq);
                }
                scrivi_hrise<<endl;
                scrivi_htwist<<endl;

		groove();
                scrivi_major_groove_width<<timestep;
                scrivi_minor_groove_width<<timestep;
                scrivi_major_groove_depth<<timestep;
                scrivi_minor_groove_depth<<timestep;
                for (int iseq=1; iseq<=lseq; iseq++)
                {
                        scrivi_major_groove_width<<" "<<0.1*grooves[0][iseq];
                        scrivi_minor_groove_width<<" "<<0.1*grooves[1][iseq];
                        scrivi_major_groove_depth<<" "<<0.1*grooves[2][iseq];
                        scrivi_minor_groove_depth<<" "<<0.1*grooves[3][iseq];
                }
                scrivi_major_groove_width<<endl;
                scrivi_minor_groove_width<<endl;
                scrivi_major_groove_depth<<endl;
                scrivi_minor_groove_depth<<endl;
		flag = ReadDump(-1, 1e10, 0);
	}

	//Close files
	scrivi_diametro.close();
	scrivi_hrise.close();
	scrivi_htwist.close();
	scrivi_major_groove_width.close();
	scrivi_minor_groove_width.close();
	scrivi_major_groove_depth.close();
	scrivi_minor_groove_depth.close();
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
                                nph = lseq-1;
                                tipologia.resize(natoms+1);
                                for (int k1=0; k1<3; k1++)
                                {
                                        coordinates[k1].resize(natoms+1);
                                        for (int ss=0; ss<=1; ss++)
                                        {
                                                Phosphates[ss][k1].resize(nph+2);
                                                amat[ss][k1].resize(nph+2);
                                                bmat[ss][k1].resize(nph+2);
                                                cmat[ss][k1].resize(nph+2);
                                                dmat[ss][k1].resize(nph+2);
                                        }
                                }
                                for (int k1=0; k1<4; k1++)
                                        grooves[k1].resize(lseq);
				hrise.resize(lseq-1);
				htwist.resize(lseq-1);
                                t1t2mat[0][0].resize(nph+2);
                                t1t2mat[0][1].resize(nph+2);
                                t1t2mat[1][0].resize(nph+2);
                                t1t2mat[1][1].resize(nph+2);
                                dtlist[0].resize(nph+2);
                                dtlist[1].resize(nph+2);
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
			//Phosphates
			indice = 0;
                        for (int i=1; i<=natoms/2; i++)
                                if (tipologia[i] == 2)
                                {
                                        indice++;
                                        Phosphates[0][0][indice] = coordinates[0][i];
                                        Phosphates[0][1][indice] = coordinates[1][i];
                                        Phosphates[0][2][indice] = coordinates[2][i];
                                }
                        for (int i=natoms/2+1; i<=natoms; i++)
                                if (tipologia[i] == 2)
                                {
                                        Phosphates[1][0][indice] = coordinates[0][i];
                                        Phosphates[1][1][indice] = coordinates[1][i];
                                        Phosphates[1][2][indice] = coordinates[2][i];
                                        indice--;
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

void cyl_info(int nucinizio, int nucfine)
{
	//compute helical axis considering cylinder wrapping the phosphates
	//precomputation
	int strand;
	int iph1 = nucinizio;
	int iph2 = nucfine - 1;
	int nph_loc = 2*(iph2 - iph1 + 1);
	double Ave[3] = {0};
	double ph[nph_loc][3];
	for (strand = 0; strand <= 1; strand++)
		for (int k1=0; k1<3; k1++)
			for (int i=iph1; i<=iph2; i++)
				Ave[k1] += Phosphates[strand][k1][i];
	for (int k1=0; k1<3; k1++)
		Ave[k1] *= 1./nph_loc;
	int ii = 0;
	for (strand = 0; strand <= 1; strand++)
	{
		for (int i=iph1; i<=iph2; i++)
		{
			for (int k1=0; k1<3; k1++)
			{
				ph[ii][k1] = Phosphates[strand][k1][i] - Ave[k1];
			}
			ii++;
		}
	}
	
	double psi[nph_loc][6];
	for (int i=0; i<nph_loc; i++)
	{
		psi[i][0] = ph[i][0]*ph[i][0];
		psi[i][1] = 2*ph[i][0]*ph[i][1];
		psi[i][2] = 2*ph[i][0]*ph[i][2];
		psi[i][3] = ph[i][1]*ph[i][1];
		psi[i][4] = 2*ph[i][1]*ph[i][2];
		psi[i][5] = ph[i][2]*ph[i][2];
	}
	double mu[6];
	for (int k=0; k<6; k++)
		mu[k] = 0;
	for (int i=0; i<nph_loc; i++)
		for (int k=0; k<6; k++)
			mu[k] += psi[i][k];
	for (int k=0; k<6; k++)
		mu[k] *= 1./nph_loc;

	for (int i=0; i<nph_loc; i++)
                for (int k=0; k<6; k++)
			psi[i][k] -= mu[k];
	
	double F0[3][3] = {0};
	for (int i=0; i<nph_loc; i++)
		for (int k1=0; k1<3; k1++)
			for (int k2=0; k2<3; k2++)
				F0[k1][k2] += ph[i][k1]*ph[i][k2];
	for (int k1=0; k1<3; k1++)
		for (int k2=0; k2<3; k2++)
			F0[k1][k2] *= 1./nph_loc;

	double F1[3][6] = {0};
	for (int i=0; i<nph_loc; i++)
		for (int k1=0; k1<3; k1++)
			for (int k2=0; k2<6; k2++)
				F1[k1][k2] += ph[i][k1]*psi[i][k2];
	for (int k1=0; k1<3; k1++)
		for (int k2=0; k2<6; k2++)
			F1[k1][k2] *= 1./nph_loc;

	double F2[6][6] = {0};
	for (int i=0; i<nph_loc; i++)
		for (int k1=0; k1<6; k1++)
	                for (int k2=0; k2<6; k2++)
                        	F2[k1][k2] += psi[i][k1]*psi[i][k2];
	for (int k1=0; k1<6; k1++)
                for (int k2=0; k2<6; k2++)
                        F2[k1][k2] *= 1./nph_loc;

	//fitting
	double Gwin = 1e10, Rwin = -1e10;
	double theta, phi;
	double Wvec[3];
	double Smat[3][3], Pmat[3][3] = {0}, Amat[3][3] = {0}, Ahatmat[3][3] = {0}, AhatAmat[3][3] = {0};
	double pvec[6], alphavec[3] = {0}, betavec[3] = {0};
	double traccia = 0, Gtry = 0;
	for (int iiter=1; iiter<=niter_cyl; iiter++)
	{
		//initialize
		for (int k1=0; k1<3; k1++)
		{
			for (int k2=0; k2<3; k2++)
			{
				Smat[k1][k2] = 0;
				Pmat[k1][k2] = 0;
				Amat[k1][k2] = 0;
				Ahatmat[k1][k2] = 0;
				AhatAmat[k1][k2] = 0;
			}
			alphavec[k1] = 0;
			betavec[k1] = 0;
		}
		traccia = 0;
		Gtry = 0;

		//computations
		phi = 2*M_PI*gsl_rng_uniform(rng);
		theta = acos(1-2*gsl_rng_uniform(rng));
		Wvec[0] = sin(theta)*cos(phi);
		Wvec[1] = sin(theta)*sin(phi);
		Wvec[2] = cos(theta);
		Smat[0][0] = 0;  Smat[0][1] = -Wvec[2];  Smat[0][2] = Wvec[1];
		Smat[1][0] = Wvec[2];  Smat[1][1] = 0;  Smat[1][2] = -Wvec[0];
		Smat[2][0] = -Wvec[1];  Smat[2][1] = Wvec[0];  Smat[2][2] = 0;
		Pmat[0][0] = 1;  Pmat[1][1] = 1;  Pmat[2][2] = 1;
		for (int k1=0; k1<3; k1++)
			for (int k2=0; k2<3; k2++)
				Pmat[k1][k2] -= Wvec[k1]*Wvec[k2];
		for (int k1=0; k1<3; k1++)
			for (int k2=0; k2<3; k2++)
				for (int l1=0; l1<3; l1++) for (int l2=0; l2<3; l2++)
					Amat[k1][k2] += Pmat[k1][l1]*F0[l1][l2]*Pmat[l2][k2];
		for (int k1=0; k1<3; k1++)
			for (int k2=0; k2<3; k2++)
				for (int l1=0; l1<3; l1++) for (int l2=0; l2<3; l2++)
					Ahatmat[k1][k2] += Smat[k1][l1]*Amat[l1][l2]*Smat[k2][l2];
		for (int k1=0; k1<3; k1++)
			for (int k2=0; k2<3; k2++)
				for (int l1=0; l1<3; l1++)
					AhatAmat[k1][k2] += Ahatmat[k1][l1]*Amat[l1][k2];
		traccia = AhatAmat[0][0] + AhatAmat[1][1] + AhatAmat[2][2];
		for (int k1=0; k1<3; k1++)
			for (int k2=0; k2<3; k2++)
				Ahatmat[k1][k2] *= 1./traccia;
		pvec[0] = Pmat[0][0];  pvec[1] = Pmat[0][1];  pvec[2] = Pmat[0][2];
		pvec[3] = Pmat[1][1];  pvec[4] = Pmat[1][2];  pvec[5] = Pmat[2][2];
		for (int k1=0; k1<3; k1++)
			for (int l1=0; l1<6; l1++)
				alphavec[k1] += F1[k1][l1]*pvec[l1];
		for (int k1=0; k1<3; k1++)
			for (int l1=0; l1<3; l1++)
				betavec[k1] += Ahatmat[k1][l1]*alphavec[l1];
		for (int k1=0; k1<6; k1++)
                        for (int k2=0; k2<6; k2++)
					Gtry += pvec[k1]*F2[k1][k2]*pvec[k2];
		if (Gtry < Gwin)
		{
			Gwin = Gtry;
			for (int k1=0; k1<3; k1++)
			{
				haxis[k1] = Wvec[k1];
				centro[k1] = betavec[k1] + Ave[k1];
			}
			Rwin = 0;
			for (int k1=0; k1<6; k1++)
				Rwin += pvec[k1]*mu[k1];			
			for (int k1=0; k1<3; k1++)
				Rwin += betavec[k1]*betavec[k1];
			Rwin = sqrt(Rwin);
		}
	}

	double cnt1[3] = {0}, cnt2[3] = {0};
	cnt(cnt1, nucinizio);	
	cnt(cnt2, nucfine);
	double controllasegno = 0;
	for (int k1=0; k1<3; k1++)
		controllasegno += (cnt2[k1]-cnt1[k1])*haxis[k1];
	if (controllasegno < 0)
		for (int k1=0; k1<3; k1++)
			haxis[k1] *= -1;
	diameter = 0.2*Rwin;

	return;
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

double compute_htwist(int nucinizio)
{
	double vec1[3], vec2[3];
	bpvec(vec1, nucinizio);	
	bpvec(vec2, nucinizio+1);
	double sp1 = dot_product(vec1, haxis, 3);
	double sp2 = dot_product(vec2, haxis, 3);
	for (int k1=0; k1<3; k1++)
	{
		vec1[k1] -= sp1*haxis[k1];
		vec2[k1] -= sp2*haxis[k1];
	}
	return angle_vec_sign(vec1, vec2, haxis);
}

double compute_hrise(int nucinizio)
{
    //h-rise is projection of dvec on helical axis
    double dvec[3], vec1[3], vec2[3];
    cnt(vec1, nucinizio);
    cnt(vec2, nucinizio+1);
    for (int k1=0; k1<3; k1++)
            dvec[k1] = vec2[k1] - vec1[k1];
    return dot_product(dvec, haxis, 3);
}

//centripetal catmull-rom spline
//the point corresponding to parameter t in [t1,t2] is avec*t**3+bvec*t**2+cvec*t+dvec; t1 and t2 are stored in t1t2
//i1 corresponds to the first point on the curve, i.e. corresponding to t1
void spline(int strand, int iph)
{
        double tvec[4] = {0};

        // determine parameters t
        double vec[3];
        double alpha = 0.5;
        for (int i=1; i<4; i++)
        {
                for (int k1=0; k1<3; k1++)
                {
                        vec[k1] = Phosphates[strand][k1][iph-1+i] - Phosphates[strand][k1][iph-1+i-1];
                }
                tvec[i] = tvec[i-1] + pow(dot_product(vec, vec, 3), 0.5*alpha);
        }

        // determine spline coefficients
        t1t2mat[strand][0][iph] = tvec[1];
        t1t2mat[strand][1][iph] = tvec[2];

        double P0[3], P1[3], P2[3], P3[3];
        double t0 = tvec[0];
        double t1 = tvec[1];
        double t2 = tvec[2];
        double t3 = tvec[3];
        for (int k1=0; k1<3; k1++)
        {
                P0[k1] = Phosphates[strand][k1][iph-1];
                P1[k1] = Phosphates[strand][k1][iph];
                P2[k1] = Phosphates[strand][k1][iph+1];
                P3[k1] = Phosphates[strand][k1][iph+2];
        }
	
	for (int k1=0; k1<3; k1++)
        {
                dmat[strand][k1][iph] = (P3[k1]*(t0 - t1)*t1*t1*(t0 - t2)*(t1 - t2)*t2 - P2[k1]*(t0 - t1)*t1*(t1 - t3)*(-(t1*t2*t3) + t0*( t2*t2 + t1*t3 - t2*t3)) +
     t2*(t2 - t3)*(-(P0[k1]*t1*(t1 - t2)*t2*(t1 - t3)) + P1[k1]*(t0 - t2)*(-( t1*t1*t3) + t0*(-(t2*t3) + t1*(t2 + t3)))))/((t0 - t1)*(t0 - t2)*pow(t1 - t2,2)*(t1 - t3)*(t2 - t3));

                cmat[strand][k1][iph] = (-(P3[k1]*(t0 - t1)*t1*(t0 - t2)*(pow(t1,2) + t1*t2 - 2*pow(t2,2))) + P2[k1]*(t0 - t1)*(t1 - t3)*(t1*t2*(-t1 + t2 - 3*t3) + t0*(pow(t1,2) + t2*(t2 - t3) + t1*(t2 + t3))) -  (t2 - t3)*(P0[k1]*t2*(-2*pow(t1,2) + t1*t2 + pow(t2,2))*(t1 - t3) + P1[k1]*(t0 - t2)*(t1*t2*(t2 - t3) - pow(t2,2)*t3 - pow(t1,2)*(t2 + t3) + t0*(-(t2*t3) + t1*(3*t2 + t3)))))/((t0 - t1)*(t0 - t2)*pow(t1 - t2,2)*(t1 - t3)*(t2 - t3));

                bmat[strand][k1][iph] = (P3[k1]*(t0 - t1)*(t0 - t2)*(2*pow(t1,2) - t1*t2 - pow(t2,2)) - P2[k1]*(t0 - t1)*(t1 - t3)*(t0*(2*t1 + t2) + t2*(t2 - 2*t3) - t1*(t2 + t3)) +  (t2 - t3)*(-(P0[k1]*(pow(t1,2) + t1*t2 - 2*pow(t2,2))*(t1 - t3)) + P1[k1]*(t0 - t2)*(-pow(t1,2) + t1*t2 + t0*(2*t1 + t2) - t1*t3 - 2*t2*t3)))/((t0 - t1)*(t0 - t2)*pow(t1 - t2,2)*(t1 - t3)*(t2 - t3));

                amat[strand][k1][iph] = (-(P3[k1]*(t0 - t1)*(t0 - t2)*(t1 - t2)) + P2[k1]*(t0 - t1)*(t0 - t3)*(t1 - t3) - (P1[k1]*(t0 - t2)*(t0 - t3) - P0[k1]*(t1 - t2)*(t1 - t3))*(t2 - t3))/ ((t0 - t1)*(t0 - t2)*pow(t1 - t2,2)*(t1 - t3)*(t2 - t3));
        }

        return;
}

void interpolazione(double tglob, double *punto, int strand)
{
        int iph = 1;
        double tt = tglob;
        while (tt > dtlist[strand][iph])
        {
                tt -= dtlist[strand][iph];
                iph++;
        }
        double tloc = tt + t1t2mat[strand][0][iph];
        for (int k1=0; k1<3; k1++)
                punto[k1] = amat[strand][k1][iph]*tloc*tloc*tloc + bmat[strand][k1][iph]*tloc*tloc + cmat[strand][k1][iph]*tloc + dmat[strand][k1][iph];
        return;
}

void groove()
{
	//spline for 0---1 is obtained by adding a -1 phosphate obtained by extrapolating from hrise and htwist of step 0-1
	double hr = compute_hrise(1);
	double ht = compute_htwist(1);
	double pscal, veccom1[3], veccom2[3];
	int strand;
	for (strand=0; strand<2; strand++)
	{
		for (int k1=0; k1<3; k1++)
			veccom1[k1] = Phosphates[strand][k1][1] - centro[k1];
		pscal = dot_product(haxis, veccom1, 3);
		for (int k1=0; k1<3; k1++)
			veccom1[k1] -= haxis[k1]*pscal;
		rotate(veccom1, haxis, -ht, veccom2);
		for (int k1=0; k1<3; k1++)
			Phosphates[strand][k1][0] = Phosphates[strand][k1][1] + veccom2[k1] - veccom1[k1] - haxis[k1]*hr;
	}

	//spline for nph-2--nph-1 is obtained by adding a nph phosphate obtained by extrapolating from hrise and htwist of corresponding step 
	hr = compute_hrise(lseq-1);
	ht = compute_htwist(lseq-1);
	for (strand=0; strand<2; strand++)
	{
		for (int k1=0; k1<3; k1++)
			veccom1[k1] = Phosphates[strand][k1][nph] - centro[k1];
		pscal = dot_product(haxis, veccom1, 3);
		for (int k1=0; k1<3; k1++)
			veccom1[k1] -= haxis[k1]*pscal;
		rotate(veccom1, haxis, -ht, veccom2);
		for (int k1=0; k1<3; k1++)
			Phosphates[strand][k1][nph+1] = Phosphates[strand][k1][nph] + veccom2[k1] - veccom1[k1] - haxis[k1]*hr;
	}
	//Interpolate phosphates
	for (int iph=0; iph<nph; iph++)
	{
                spline(0, iph);
                spline(1, iph);
	}

	//create dtlists and return tmax
	tmax[0] = 0; tmax[1] = 0;
        for (int iph=1; iph<nph; iph++)
        {
                dtlist[0][iph] = t1t2mat[0][1][iph] - t1t2mat[0][0][iph];
                dtlist[1][iph] = t1t2mat[1][1][iph] - t1t2mat[1][0][iph];
                tmax[0] += dtlist[0][iph];
                tmax[1] += dtlist[1][iph];
        }
	//compute grooves
	double punto[3];
	double tcom;
	double widdepr[2], widdepl[2];
	int iph;
	for (int iseq=1; iseq<=lseq; iseq++)
	{
		if (iseq == 1)
		{
			//extrapolate from first spline
			iph = 1;
			tcom = 1.5*t1t2mat[0][0][iph] - 0.5*t1t2mat[0][1][iph];
			for (int k1=0; k1<3; k1++)
                                punto[k1] = amat[0][k1][iph]*tcom*tcom*tcom + bmat[0][k1][iph]*tcom*tcom + cmat[0][k1][iph]*tcom + dmat[0][k1][iph];
		} 
		else if (iseq == lseq)
		{
			//extrapolate from last spline
			iph = lseq - 2;
			tcom = -0.5*t1t2mat[0][0][iph] + 1.5*t1t2mat[0][1][iph];
			for (int k1=0; k1<3; k1++)
                                punto[k1] = amat[0][k1][iph]*tcom*tcom*tcom + bmat[0][k1][iph]*tcom*tcom + cmat[0][k1][iph]*tcom + dmat[0][k1][iph];
		}
		else
		{
			iph = iseq-1;
			//determine reference point on first strand as mid-point between phosphates sandwiching the base pair
			tcom = 0.5*(t1t2mat[0][0][iph] + t1t2mat[0][1][iph]);
			for (int k1=0; k1<3; k1++)
				punto[k1] = amat[0][k1][iph]*tcom*tcom*tcom + bmat[0][k1][iph]*tcom*tcom + cmat[0][k1][iph]*tcom + dmat[0][k1][iph];
		}

		find_minimum_right(widdepr, iph, punto);
		find_minimum_left(widdepl, iph, punto);

		//determine which is which
		if (widdepr[0] < widdepl[0])
		{
			grooves[0][iseq] = widdepl[0];
			grooves[1][iseq] = widdepr[0];
			grooves[2][iseq] = widdepl[1];
			grooves[3][iseq] = widdepr[1];
		}
		else
		{
			grooves[0][iseq] = widdepr[0];
			grooves[1][iseq] = widdepl[0];
			grooves[2][iseq] = widdepr[1];
			grooves[3][iseq] = widdepl[1];
		}

	}
	
	//close to the borders there will be problems, and only one groove is computed
	//to understand which groove it is, we compare it to a reference obtained as average of the other values
	double refmaj = 0, refmin = 0;
	int refnorm = 0;
	for (int iseq=1; iseq<=lseq; iseq++)
	{
		if (grooves[0][iseq]>0 && grooves[1][iseq]>0)
		{
			refmaj += grooves[0][iseq];
			refmin += grooves[1][iseq];
			refnorm ++;
		}
	}
	refmaj *= 1./refnorm;
	refmin *= 1./refnorm;
	for (int iseq=1; iseq<=lseq; iseq++)
        {
		//select only cases in which exactly one groove was computed correctly
		//the algorithm has automatically assigned the value to the major groove... let's check and correct if necessary
                if (grooves[0][iseq]>0 && grooves[1][iseq]<0)
                {
			//we check whether our groove is closer to refmaj or refmin
			if (fabs(grooves[0][iseq]-refmin) < fabs(grooves[0][iseq]-refmaj))
			{
				grooves[1][iseq] = grooves[0][iseq];
				grooves[3][iseq] = grooves[2][iseq];
				grooves[0][iseq] = -1e10;
				grooves[2][iseq] = -1e10;
			}
                }
        }

	return;
}

void find_minimum_right(double *widdep, int iph, double *punto)
{
	//determine starting point on second strand to look for the minima corresponding to grooves widt_precisionhs
	double tcom = 0; 
	for (int jph=1; jph<iph; jph++)
                tcom += dtlist[1][jph];
        tcom += 0.5*dtlist[1][iph];
	double distvec[3], distvecold[3], distold, distnew;
	interpolazione(tcom, distvec, 1);
	diff_vec(distvec, punto, distvec, 3);
	distold = dot_product(distvec, distvec, 3);

	//look for minimum in forward direction (increasing iph)
	tcom += dt_precision;
	interpolazione(tcom, distvec, 1);
	diff_vec(distvec, punto, distvec, 3);
	distnew = dot_product(distvec, distvec, 3);
	if (distnew < distold) //in this case, we are already travelling towards a minimum
	{
		while (1)
		{
			distold = distnew;
			tcom += dt_precision;
			if (tcom > tmax[1]) //we reached the end of the road
			{
				widdep[0] = -1e10;
				widdep[1] = -1e10;
				return;
			}
			cp_vec(distvec, distvecold, 3);
			interpolazione(tcom, distvec, 1);
			diff_vec(distvec, punto, distvec, 3);
			distnew = dot_product(distvec, distvec, 3);
			if (distnew > distold)
			{
				widdep[0] = sqrt(distold) - groove_width_shift;
				//determine corresponding depth
				double puntodepth[3];
				cnt(puntodepth, iph+1); //take center of base pair
				for (int k1=0; k1<3; k1++)
					puntodepth[k1] -= punto[k1] + 0.5*distvecold[k1]; //take separation from midpoint of groove widt_precisionh vector
				widdep[1] = sqrt(dot_product(puntodepth, puntodepth, 3)) - groove_depth_shift;
				return;
			}
		}
	}
	else  //in this case, first we need to pass by the maximum, then we go to the minimum
	{
		//find maximum
		while (1)
		{
			distold = distnew;
			tcom += dt_precision;
			if (tcom > tmax[1]) // we reached the end of the road
			{
				widdep[0] = -1e10;
				widdep[1] = -1e10;
				return;
			}
			cp_vec(distvec, distvecold, 3);
                        interpolazione(tcom, distvec, 1);
			diff_vec(distvec, punto, distvec, 3);
			distnew = dot_product(distvec, distvec, 3);
			if (distnew < distold)
			{
				break;
			}
		}
		//find minimum
		while (1)
		{
			distold = distnew;
			tcom += dt_precision;
			if (tcom > tmax[1]) // we reached the end of the road
			{
				widdep[0] = -1e10;
				widdep[1] = -1e10;
				return;
			}
			cp_vec(distvec, distvecold, 3);
                        interpolazione(tcom, distvec, 1);
			diff_vec(distvec, punto, distvec, 3);
			distnew = dot_product(distvec, distvec, 3);
			if (distnew > distold)
			{
				widdep[0] = sqrt(distold) - groove_width_shift;
				//determine corresponding depth
				double puntodepth[3];
				cnt(puntodepth, iph+1); //take center of base pair
				for (int k1=0; k1<3; k1++)
					puntodepth[k1] -= punto[k1] + 0.5*distvecold[k1]; //take separation from midpoing of groove widt_precisionh vector
				widdep[1] = sqrt(dot_product(puntodepth, puntodepth, 3)) - groove_depth_shift;

				return;
			}
		}
	}
	return;
}

void find_minimum_left(double *widdep, int iph, double *punto)
{
	//determine starting point on second strand to look for the minima corresponding to grooves widt_precisionhs
	double tcom = 0;
        for (int jph=1; jph<iph; jph++)
                tcom += dtlist[1][jph];
        tcom += 0.5*dtlist[1][iph];
	double distvec[3], distvecold[3], distold, distnew;
        interpolazione(tcom, distvec, 1);
	diff_vec(distvec, punto, distvec, 3);
	distold = dot_product(distvec, distvec, 3);

	//look for minimum in backward direction (decreasing iph)
	tcom -= dt_precision;
        interpolazione(tcom, distvec, 1);
	diff_vec(distvec, punto, distvec, 3);
	distnew = dot_product(distvec, distvec, 3);
	if (distnew < distold) //in this case, we are already travelling towards a minimum
	{
		while (1)
		{
			distold = distnew;
			tcom -= dt_precision;
			if (tcom < 0) // we reached the end of the road
			{
				widdep[0] = -1e10;
				widdep[1] = -1e10;
				return;
			}
			cp_vec(distvec, distvecold, 3);
                        interpolazione(tcom, distvec, 1);
			diff_vec(distvec, punto, distvec, 3);
			distnew = dot_product(distvec, distvec, 3);
			if (distnew > distold)
			{
				widdep[0] = sqrt(distold) - groove_width_shift;
				//determine corresponding depth
				double puntodepth[3];
				cnt(puntodepth, iph+1); //take center of base pair
				for (int k1=0; k1<3; k1++)
					puntodepth[k1] -= punto[k1] + 0.5*distvecold[k1]; //take separation from midpoing of groove widt_precisionh vector
				widdep[1] = sqrt(dot_product(puntodepth, puntodepth, 3)) - groove_depth_shift;

				return;
			}
		}
	}
	else  //in this case, first we need to pass by the maximum, then we go to the minimum
	{
		//find maximum
		while (1)
		{
			distold = distnew;
			tcom -= dt_precision;
			if (tcom < 0) // we reached the end of the road
			{
				widdep[0] = -1e10;
				widdep[1] = -1e10;
				return;
			}
			cp_vec(distvec, distvecold, 3);
                        interpolazione(tcom, distvec, 1);
			diff_vec(distvec, punto, distvec, 3);
			distnew = dot_product(distvec, distvec, 3);
			if (distnew < distold)
			{
				break;
			}
		}
		//find minimum
		while (1)
		{
			distold = distnew;
			tcom -= dt_precision;
			if (tcom < 0) // we reached the end of the road
			{
				widdep[0] = -1e10;
				widdep[1] = -1e10;
				return;
			}
			cp_vec(distvec, distvecold, 3);
                        interpolazione(tcom, distvec, 1);
			diff_vec(distvec, punto, distvec, 3);
			distnew = dot_product(distvec, distvec, 3);
			if (distnew > distold)
			{
				widdep[0] = sqrt(distold) - groove_width_shift;
				//determine corresponding depth
				double puntodepth[3];
				cnt(puntodepth, iph+1); //take center of base pair
				for (int k1=0; k1<3; k1++)
					puntodepth[k1] -= punto[k1] + 0.5*distvecold[k1]; //take separation from midpoing of groove widt_precisionh vector
				widdep[1] = sqrt(dot_product(puntodepth, puntodepth, 3)) - groove_depth_shift;

				return;
			}
		}
	}
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
