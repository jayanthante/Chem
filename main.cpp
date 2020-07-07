#include<iostream>
#include<fstream>
#include<string.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<complex>
#include<iomanip>
using namespace std ;
#define pi 3.14159265358979323846

//ofstream overlap("overlap.txt", ios::out);
// ofstream meanr("rmean.txt", ios::out);
// ofstream derivative("deriv1.txt", ios::out);
// ofstream couple("coupling.txt", ios::out);
ofstream eigens("eigen.txt", ios::out);
const int netgrid = 50000 ;
int xgrid = 1024 ;
#include "numerov.cpp"
#include "time.cpp"

double H_X(double x);
double H_analytic(double x, int surf);
double dipol(double x, int surf);
double coupling(double x, int surf1, int surf2);
double hamilpolate(double xx, double* x, double* v, int grid );
void reada(double* v, double* x, int gridi, string fil);

int main(){
	double x, xmin = 1.6, xmax = 7.0;
	double h = (xmax-xmin)/(netgrid-1.0);
	double dx = (xmax-xmin)/(xgrid-1.0);
	class numerov n2 ;
	int nstates = 7;
	
	string surface[] = {"X", "b", "c", "o", "bd", "cd", "ed"} ;//names of states
	int statecount[] = {35, 25, 55, 60, 35, 55, 50}; //number of vibrational elvels; 
	int Lang[] = {0, 1, 1, 1, 1, 1};//whether excitation allowed 1 -> excitation possible
	n2.N = netgrid ;
	static double hmatrix[7][7][1024];
	for ( int i = 0 ; i < xgrid ; i++ ){
		for ( int is = 0 ; is < nstates ; is++ ){
			for ( int j = 0 ; j < nstates ; j++ ){
				hmatrix[is][j][i] = 0.0; 
			}
		}
	}

	for ( int i = 0 ; i < xgrid ; i++ ){
		x = xmin +(i*dx);
		for ( int is = 0 ; is < 7 ; is++ ){
			for ( int j = 0 ; j < 7 ; j++ ){
				if ( is != j ){
					hmatrix[is][j][i] = coupling(x,is,j);
				}
			}
		}
	}
	//X - 0, b - 1, c - 2, o - 3, bd - 4, cd - 5, ed - 6
	double cm_au = 1.0/219474.544729040;

	for ( int i = 0 ; i < xgrid ; i++ ){
		x = xmin +(i*dx);
		for ( int is = 0 ; is < 7 ; is++ ){
			if ( is == 0)hmatrix[is][is][i] = H_X(x ); //ground state
			if ( is >= 1 ){
				hmatrix[is][is][i] = H_analytic(x,is); //excited states
			}
		}
	}


	
	int netstates = 0 ; //total number of states being used
	for ( int i = 0 ; i < nstates ; i++ ) netstates = netstates + statecount[i] ;
	
	//bd, cd, ed, b, c, o - parameters from Spelsberg, Meyer
	//X - parameter from "Morse Potential for OO, NN, and NO" by Konowalow and Hirschfelder
	double* phi = new double[netgrid*netstates] ;
	// =(double*)malloc(sizeof(double)*netgrid*netstates); //the wavefunction grid for 35*7 states
	double* eigen =(double*)malloc(sizeof(double)*netstates); //the eigen values for 35*7 states
	
	//c values from Spelsberg Meyer, also for X
	double* curve = new double[netgrid] ;
	double* curved = new double[xgrid] ;
	double* xd = new double[xgrid] ;
	
	int states ;
	for ( int i = 0 ; i < nstates ; i++ ){
		states = 0 ;
		for ( int i1 = 0 ; i1 < i ; i1++ ) states = states + statecount[i1] ; //to count the vibrational levels in the states so far
		for ( int j = 0 ; j < xgrid ; j++ ){
			curved[j] = hmatrix[i][i][j];
			xd[j] = xmin + j*dx;
		}
		for ( int j = 0 ; j < netgrid ; j++ ){
			curve[j] = hamilpolate(xmin + j*h, &xd[0], &curved[0], xgrid ); //linear interpolation of the PECs to suit the current grid
		}
		cout << "# SURFACE : - " << surface[i]  << endl ;
		n2.calc( &curve[0], i, statecount[i], &phi[states*netgrid], &eigen[states]); //calling the numerov calculation for surface[i]
		eigens << surface[i] << "\t" << statecount[i] << endl ; //printing the number of eigen values
		for ( int j = 0 ; j < statecount[i] ; j++ ){
			eigens << eigen[states+j] << endl ; //printing the eigen values
		}
		eigens << endl ;
	}
	
	netstates = 0 ;
	for ( int i = 0 ; i < nstates ; i++ ){
		netstates = netstates+statecount[i] ;
		cout << "# SURFACE : - " << surface[i] << "\t" << statecount[i] << endl ;
	}
	
	//Diabatic coupling parameters from Spelsberg, Meyer
	//string surface[] = {"b", "c", "o", "bd", "cd", "ed"} ;

	double** Ematrix = new double*[netstates] ; //Hamiltonian matrix, including coupling
	double** dipole = new double*[netstates] ;
	
	for ( int i = 0 ; i < netstates ; i++ ){
		Ematrix[i] = new double[netstates] ;
		dipole[i] = new double[netstates] ;
	}
	
	for ( int i = 0 ; i < netstates ; i++ ){
		for ( int j = 0 ; j < netstates ; j++ ){
			Ematrix[i][j] = 0.0 ;//defining all to be zero, just in case
		}
	}
	
	int k1 = 0, k2 ;
	ofstream hprint("hamiltonian.txt", ios::out);
	for ( int i = 0 ; i < nstates ; i++ ){
		states = 0 ;
		for ( int i1 = 0 ; i1 < i ; i1++ )states = states+statecount[i1] ;
		for ( int j = 0 ; j < statecount[i] ; j++ ){
			k2 = states+j ;
			Ematrix[k1][k1] = eigen[k2] ; //assigning diagonal elements to the eigen values
			k1++ ;
		}
	}
	
	int c_colin, c_rowin ; //colin = column index, rowin = row index
	char statefile1[3] ;
	char statefile2[3] ;
	double cvalue, dumdip ;
	int states1 = 0 ;
	for ( int i = 0 ; i < nstates ; i++ ){ //el.state 1 index
		for ( int j = 0 ; j < statecount[i] ; j++ ){//vib levels within el state 1
			for ( int k = 0 ; k < nstates ; k++ ){//el. state 2 index
				for ( int l = 0 ; l < statecount[k] ; l++ ){ //vib levels within el state 2
					states = 0 ;
					states1 = 0 ;
					for ( int i1 = 0 ; i1 < i ; i1++ ) states = states + statecount[i1] ; //net vib levels before el.state 1
					for ( int i1 = 0 ; i1 < k ; i1++ ) states1 = states1 + statecount[i1] ; //net vib levels before el.state 2
					c_rowin = states + j ;//index of vibrational level j (within i) among all vib levels
					c_colin = states1 + l ; //index of vibrational level l (within k) among all vib levels
					if ( i == 0 && Lang[k] == 1 ){ //if one state is ground state and the other transition allowed
						dumdip = 0.0; //dummy for dipole
						for ( int i1 = 0 ; i1 < netgrid ; i1++ ){
							dumdip += phi[c_rowin*netgrid + i1 ]*dipol(xmin+i1*h,k)*phi[c_colin*netgrid + i1];
						}
						dipole[c_rowin][c_colin] = dumdip ;
					} else if ( Lang[i] == 1  && k == 0 ){
						dumdip = 0.0;
						for ( int i1 = 0 ; i1 < netgrid ; i1++ ){
							dumdip += phi[c_rowin*netgrid + i1 ]*dipol(xmin+i1*h,i)*phi[c_colin*netgrid + i1];
						}
						dipole[c_rowin][c_colin] = dumdip ;
					}else{
						dipole[c_rowin][c_colin] = 0.0 ;
					}
	
					if ( i != k ){
						cvalue = 0.0; // value of <chi_ij|H_ijkl|chi_kl> to be integrated over the whole grid
						for ( int i1 = 0 ; i1 < netgrid ; i1++ ){
							cvalue += phi[c_rowin*netgrid + i1 ]*hmatrix[i][k][i1]*phi[c_colin*netgrid + i1];
						}
						Ematrix[c_rowin][c_colin] = cvalue ;
					}
				}
			}
		}
	}
	
	for ( int i = 0 ; i < netstates ; i++ ){
		for ( int j = 0 ; j < netstates ; j++ ){
			if ( i == j && j < statecount[0] ) Ematrix[i][j] = Ematrix[i][j] - eigen[0] ; //printing the Hamiltonian
			hprint << setw(12) << setprecision(17) << scientific << Ematrix[i][j] << " " ;
		}
		hprint << endl ;
	}
	cout << "# Dynamics begin" << endl ;
	gear N2 ;
	N2.net = netstates ;
	N2.roll(&Ematrix[0], &phi[0], &surface[0], &Lang[0], &statecount[0], &dipole[0]);
	
	delete[] phi ;
	delete[] Ematrix ;

}

double H_X(double x ){
//	double xeq = 2.0676620680, w = 0.010752551, wx = 6.58748/100000.0;
//	double De = w*w/(4.0*wx);
//	double beta = w*sqrt(0.5*mass/De);
	double xeq = 2.04158, beta = 1.50765, De = 0.363965;
	double V =(1.0-exp(-beta*(x-xeq)));
	V = V*V*De;
	return V;
}

double H_analytic(double x, int surf ){
	double xeq, beta, alfa;
	double* c;
	double cdum[6][8] ={
		{0.458401, 0.0, 0.0206768, -0.0118135, 0.115635, -0.0842532, -0.00355449, 0.0412016},
		{0.471233, 0.0, 0.319474, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.477476, 0.0, 0.316443, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.470701, 0.0, 0.0542505, 0.021976, 0.0212176, 0.0111611, 0.00161939, 0.0},
		{0.472326, 0.0, 0.338606, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.522845, 0.0, 0.319334, 0.0, 0.0, 0.0, 0.0, 0.0},
	} ;
	switch(surf ){
	case 1: //b
		xeq	= 2.43702;
		beta	= 2.06197;
		alfa	= 0.5;
		c	= cdum[0];
		break;
	case 2: //c
		xeq	= 2.11006; 
		beta	= 1.43490;
		alfa	= 0.0;
		c	= cdum[1];
		break;
	case 3: //o
		xeq	= 2.22913; 
		beta	= 1.23318;
		alfa	= 0.0;
		c	= cdum[2];
		break;
	case 4: //bd
		xeq	= 2.73132; 
		beta	= 1.27484;
		alfa	= 0.1;
		c	= cdum[3];
		break;
	case 5: //cd
		xeq	= 2.12246;
		beta	= 1.36036;
		alfa	= 0.0;
		c	= cdum[4];
		break;
	case 6: //ed
		xeq	= 2.10945;
		beta	= 1.42744;
		alfa	= 0.0;
		c	= cdum[5];
		break;
	}
	double r,V ;
	r =(1.0-exp(-beta*(x-xeq)))/(1.0+(alfa*exp(-beta*(x-xeq))));
	V = 0.0 ;
	for ( int i = 0 ; i < 8 ; i++ )V = V+(c[i]*pow(r,i));
	if ( surf == 1 && x <= 1.819 ){
		V = 0.0 ; 	//b LJ patch
		double depth = 0.0779 ;
		double assimptot = 0.536 ;
		double n = 8 ;
		double m = 1.8 ;
		double r_patch = 1.819 ;
		V =(depth/(n-m))*(m*pow(xeq/x,n)- n*pow(xeq/x,m)) + assimptot ;
	}
	if ( surf == 4 && x <= 1.7 ){
		V = 0.0 ; 	//bd LJ patch
		double depth = 0.1097 ;
		double assimptot = 0.5804 ;
		double n = 2.71 ;
		double m = 1.8 ;
		double r_patch = 1.7 ;
		V = (depth/(n-m))*(m*pow(xeq/x,n)- n*pow(xeq/x,m)) + assimptot ;
	}
	return V ;
}

double coupling(double x, int surf1, int surf2 ){
	double v = 0.0 ;
	double xeq = 0.0, beta = 0.0, c = 0.0 ;
	if ( (surf1 == 1 && surf2 == 2)||(surf2 == 1 && surf1 == 2) ){//b-c
		xeq	= 1.71141 ;
		beta	= 2.21929 ;
		c	= -5.64067e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq));
	}else if ( (surf1 == 1 && surf2 == 3)||(surf2 == 1 && surf1 == 3) ){//b-o
		xeq	= 1.82855 ;
		beta	= 3.65381 ;
		c	= 2.18468e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq));
	}else if ( (surf1 == 2 && surf2 == 3)||(surf2 == 2 && surf1 == 3) ){//c-o
		xeq	= 1.26651 ;
		beta	= 1.90404 ;
		c	= 1.49604e-2 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq));
	}else if ( (surf1 == 4 && surf2 == 5)||(surf2 == 4 && surf1 == 5) ){//bd-cd
		xeq	= 1.85220 ;
		beta	= 5.78569 ;
		c	= -1.31280e-2 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq));
	}else if ( (surf1 == 4 && surf2 == 6)||(surf2 == 4 && surf1 == 6) ){//bd-ed
		xeq	= 1.85220 ;
		beta	= 5.78569 ;
		c	= -4.60098e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq));
	}else if ( (surf1 == 5 && surf2 == 6)||(surf2 == 5 && surf1 == 6) ){//cd-ed
		xeq	= 2.02162 ; 
		beta	= 3.39115 ; 
		c	= 2.85983e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq));
	}
	return v ;
}


void reada( double* v, double* x, int gridi, string fil ){
	char file[3];
	int gir;
	if ( fil == "X_0" ){
		gir = xgrid;
	}else{
		gir = 2000;
	}
	for ( int i = 0 ; i < gir ; i++ ){
		v[i] = 0.0;
		x[i] = 0.0;
	}
	ifstream rea(fil.c_str(), ios::in);
	for ( int i = 0 ; i < gridi ; i++ ){
		rea >> x[i] >> v[i] ;
	}
}

double hamilpolate(double xx, double* x, double* v, int grid ){
	//linear interpolation
	double energy = 0.0 ;
	if ( xx >= x[0] && xx <= x[grid-1] ){
		for ( int j = 0 ; j < grid ; j++ ){
			if ( xx >= x[j] && xx <= x[j+1] ){
				energy = (v[j+1] - v[j])*(xx - x[j])/(x[j+1] - x[j]) + v[j] ;
				return energy;
			}
		}
	}else{
		return energy;
	}
}

double dipol(double x, int surf){
	double moment = 0.0 ;
	if ( surf == 1 ){ //b
		moment = (0.4375*x)-1.5375 ;
	}
	if ( surf == 2 ){ //c
		moment = (0.125*x)-0.525 ;
	}
	if ( surf == 3 ){ //o
		moment = (0.25*x)-0.05 ;
	}
	if ( surf == 4 ){ //bd
		 moment = 2.4*x*x-11.04*x+11.4 ;
	}
	if ( surf == 5 ){ //cd
		 if(x<2.05){ moment = -0.6 ;}
		 if(x>=2.05){ moment = -0.4*x+0.1 ;}
	}
	if ( surf == 6 ){ //ed
		 if(x<2.1){ moment = -0.2 ;}
		 if(x>=2.1){ moment = -0.5*x+0.4 ;}
	}
	return moment ;
}
