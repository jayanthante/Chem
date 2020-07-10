#include<iostream>
#include<fstream>
#include<fftw3.h>
#include<math.h>
#include<complex>
#include<stdlib.h>
#include <omp.h>
using namespace std;

// OMP
#define CHUNKSIZE 43
#define NTHREADS 24

// all units in au unless mentioned otherwise
const int xgrid = 1024;
double xmin = 1.5, xmax = 8.0;
double dx = (xmax - xmin)/1023.0;
double ixgrid = 1.0/1024.0; 
double Rexit = 5.5 ;
int iexit = ((Rexit - xmin)/dx) - 1 ;
int nstates = 7;
double cm_au = 1.0/219474.544729040; //converts cm-1 to au
double nm_au = 45.56337 ; //converst nm to au
double mass = 12763.0205 ;
int T = 8500;

//functions
void kmatrixcalc(complex<double>** KE); //calculate KE operator
double H_analytic ( double x, int surf); //PEC for excited states
double H_X ( double x); //PEC for ground state
double pulse(double t);
void reada( double* v, double* x, int gridi, string fil); //to read x and f(x)
double linpolate(double xx, double* x, double* v, int grid); //linear interpolation
void corr(double dt, double net, complex<double>* psi, complex<double>* psi1, complex<double>* psi2, complex<double>* psi3, complex<double>* psi4, complex<double>* prodpsi); //correcot routeine
void pred(double dt, double net, complex<double>* psi, complex<double>* psi1, complex<double>* psi2, complex<double>* psi3, complex<double>* psi4); //predictor rountine
double coupling(double x, int surf1, int surf2); //calculates analytical V_ij
complex<double> vimag(double x); //absorption potential
double dipole(double x, int surf); //calculates dipole

int main(){
	int nstates = 7; //number of states
	complex<double>** kmatrix = new complex<double>*[xgrid]; //to calculate KE operator matrix
	for ( int i = 0 ; i < xgrid ; i++ ) kmatrix[i] = new complex<double>[xgrid];
	kmatrixcalc(&kmatrix[0]);
//	cout << q1 << "\t" << xp << "\t" <<  << endl;
	complex<double>* phi=new complex<double>[xgrid*nstates]() ;
	complex<double> dipoleout[7][7];

	double x ;
	double vmin = 0.0, vmax = 1.50 ;
	double norm = 0.0, norm2 = 0.0;

	double* x_X = new double[1024]() ;
	double* f_X = new double[1024]() ;
	reada( &f_X[0], &x_X[0], 1024, "X_0") ; //reads the ground state wf
	for ( int i = 0 ; i < xgrid ; i++ ){
		x = xmin + (i*dx) ;
		real(phi[i*nstates]) = linpolate( x, &x_X[0], &f_X[0], 1024) ;
		norm += real(phi[i*nstates])*real(phi[i*nstates]);
	}
	for ( int i = 0 ; i < xgrid ; i++ ) real(phi[i*nstates]) = real(phi[i*nstates])/sqrt(norm);
	norm = 0.0;
	delete[] x_X ;
	delete[] f_X ;
	double dt = 0.001; //time step in au
	int nshot = 1/dt;

	double field = 0.0 ;
	double field2 = 0.0 ;

	static double vmatrix[7][7][1024];
	for ( int i = 0 ; i < xgrid ; i++ ) {
		for ( int is = 0 ; is < nstates ; is++ ) {
			for ( int j = 0 ; j < nstates ; j++ ) {
				vmatrix[is][j][i] = 0.0; 
			}
		}
	}

	for ( int i = 0 ; i < xgrid ; i++ ) {
		x = xmin + (i*dx) ;
		for ( int is = 0 ; is < 7 ; is++ ) {
			for ( int j = 0 ; j < 7 ; j++ ) {
				if ( is != j){
					vmatrix[is][j][i] = coupling(x,is,j) ;
				}
			}
		}
	}


	for ( int i = 0 ; i < xgrid ; i++ ) {
		x = xmin + (i*dx) ;
		for ( int is = 0 ; is < 7 ; is++ ) {
			if ( is == 0 ) vmatrix[is][is][i] = H_X ( x );
			if ( is >= 1 ) vmatrix[is][is][i] = H_analytic(x,is);
		}
	}
	//absorption potential does NOT work well with this method
	complex<double>* vabs= new complex<double>[xgrid]; //absorption potential
	for ( int i = 0 ; i < xgrid ; i++ ){
		vabs[i] = vimag( xmin + (i*dx) );
		//cout << (xmin + (i*dx)) << "\t" << real(vabs[i]) << "\t" << imag(vabs[i]) << endl;
	}

	complex<double>* psi = new complex<double>[xgrid*nstates];
	
	//i*hbar*dPpsi/dt,i*hbar*dQ1psi/dt,i*hbar*dQ2psi/dt,
	complex<double>* dpsi = new complex<double>[xgrid*nstates]; // 1st derivative
	complex<double>* d2psi = new complex<double>[xgrid*nstates]; // second derivative
	complex<double>* d3psi = new complex<double>[xgrid*nstates]; // third  derivative
	complex<double>* d4psi = new complex<double>[xgrid*nstates]; // fourth  derivative
	complex<double>* prodpsi = new complex<double>[xgrid*nstates]; // H*psi

	double* popa = new double[nstates]; // 

	ofstream popu("pop.txt",ios::out);
	ofstream dout("douta.dat",ios::out) ;
	popu.precision(12);
	dout.precision(10);
	for ( int i = 0 ; i < xgrid ; i++ ){
		psi[i*nstates] = phi[i*nstates];
		for ( int j = 1 ; j < nstates ; j++ ){
			dpsi[i*nstates + j] = complex<double>(0.0,0.0);
			d2psi[i*nstates + j] = complex<double>(0.0,0.0);
			d3psi[i*nstates + j] = complex<double>(0.0,0.0);
			d3psi[i*nstates + j] = complex<double>(0.0,0.0);
		}
	}
	delete[] phi;
	//TIME
	double tym = 0;

	for ( int a = 0 ; a*dt < T ; a++ ){
		tym = a*dt;
		field = pulse(tym);
		if ( a == 0 ) cout << "# Dynamics begin" << endl;
		if ( a != 0 ){
			pred(dt,nstates*xgrid,&psi[0],&dpsi[0],&d2psi[0],&d3psi[0],&d4psi[0]); // PREDICTOR after 1st step
		}

		for ( int iter = 0 ; iter < 1 ; iter++ ){ // run corrector more times; no iteration was necessary
			//OMP parallelization
			#pragma omp parallel private(x) num_threads(NTHREADS)
			{
			#pragma omp for schedule(dynamic,CHUNKSIZE)
			for ( int i = 0 ; i < xgrid ; i++ ){
			        x = xmin + (i)*dx;

			        for ( int j = 0 ; j < nstates ; j++ ){
			       	 for ( int i2 = 0 ; i2 < xgrid ; i2 ++ ){
			       		 prodpsi[nstates*i + j] += kmatrix[i][i2]*psi[nstates*i2 + j]; //T*psi
			       	 }
			       	 //prodpsi[nstates*i + j] += vabs[i]*psi[nstates*i + j];

			       	 for ( int k = 0 ; k < nstates ; k++){
			       		 prodpsi[i*nstates + j] += vmatrix[j][k][i]*psi[nstates*i + k] ; //T*psi + V*psi
			       		 if ( k == 0 && j != 0 ) prodpsi[nstates*i + j] -= dipole(x,j)*(field )*psi[nstates*i + k]; //T*psi + V*psi - mu.E(t)*psi
			       		 if ( j == 0 && k != 0 ) prodpsi[nstates*i + j] -= dipole(x,k)*(field )*psi[nstates*i + k];
			       	 }

			       	 prodpsi[nstates*i + j] = prodpsi[nstates*i + j]*complex<double>(0.0,-1.0); //-i*H*psi
			        }
			}
			}

			if ( a == 0 && iter == 0 ){ // PREDICTOR for 1st step
			        for ( int i = 0 ; i < xgrid ; i++ ){
			       	 for ( int j = 0 ; j < nstates ; j++ ){ 
			       		 dpsi[nstates*i + j] = prodpsi[nstates*i + j]; 
			       	 }
			        }
			        pred(dt,nstates*xgrid,&psi[0],&dpsi[0],&d2psi[0],&d3psi[0],&d4psi[0]);
			}

			corr(dt,nstates*xgrid,&psi[0],&dpsi[0],&d2psi[0],&d3psi[0],&d4psi[0],&prodpsi[0]); // CORRECTOR

			for ( int i = 0 ; i < xgrid ; i++ ){
			        for ( int j = 0 ; j < nstates ; j++ ){
			       	 real(prodpsi[i*nstates + j]) = 0.0;
			       	 imag(prodpsi[i*nstates + j]) = 0.0;
			        }
			}
		}

		if ( a%nshot == 0 ){
			//to calculate dipole after every nshot steps
			for ( int i1 = 0 ; i1 < 7 ; i1++ ) {
				for ( int i2 = 0 ; i2 < 7 ; i2++ ) {
					dipoleout[i1][i2] = complex<double>(0.0,0.0);
				}
			}
			for ( int i1 = 0 ; i1 < 7 ; i1++ ) {
				for ( int i2 = 0 ; i2 < 7 ; i2++ ) {
					if ( i1 > 0 && i2 == 0 ) {
						for ( int j = 0 ; j < xgrid ; j++ ){
							dipoleout[i1][i2] += (conj(psi[j*nstates + i1])*dipole(xmin + j*dx, i1)*psi[j*nstates + i2]) ;
						}
					} else if ( i2 > 0 && i1 == 0 ){
						for ( int j = 0 ; j < xgrid ; j++ ){
							dipoleout[i1][i2] += (conj(psi[j*nstates + i1])*dipole(xmin + j*dx, i2)*psi[j*nstates + i2]) ;
						}
					}
				}
			}
			//dipole(t) printing
			dout << tym;
			for ( int i = 1 ; i < 7 ; i++ ){
				dout << "\t" << real(dipoleout[0][i] + dipoleout[i][0]) << "\t" << imag(dipoleout[0][i] + dipoleout[i][0]);
			}
			dout << endl;
			//population calculation
			popu << tym << "\t" << field;
			for ( int j = 0 ; j < nstates ; j++ ) popa[j] = 0.0;
			for ( int i = 0 ; i < xgrid ; i++ ){
				for ( int j = 0 ; j < nstates ; j++ ){
					popa[j] += (real(psi[nstates*i + j])*real(psi[nstates*i + j])) + (imag(psi[nstates*i + j])*imag(psi[nstates*i + j]));
				}
			}
			norm = 0.0;
			for ( int j = 0 ; j < nstates ; j++ ) norm += popa[j];
			for ( int j = 0 ; j < nstates ; j++ ){
				popu << "\t" << popa[j] ;
				popa[j] = 0.0;
			}
			popu << "\t" << norm << endl; //this method didn't require renormalization
/*			if ( a%(nshot*10) == 0 ){
				char outfile[10];
				sprintf(outfile,"%6.2f",tym);
				const char* la = outfile;
				ofstream tout(la,ios::out);
				for ( int j = 0 ; j < xgrid ; j++ ) {
					tout << xmin + (j*dx) ;
					for ( int i = 0 ; i < nstates ; i++ ){
						tout << "\t" << real(psi[j*nstates + i]) << "\t" << imag(psi[j*nstates + i]) ;
					}
					tout << endl;
				}
			}*/
		}
	}// time loop
	delete[] popa;
	delete[] psi;
	delete[] prodpsi;
	delete[] dpsi;
	delete[] d2psi;
	delete[] d3psi;
	delete[] d4psi;
}
//==================================


double pulse(double t){
	double e, omega, t0, sigma, sigmad ;
	double e0 = sqrt(2.0*4.0*M_PI*(2.84945177972484E-06)/137.0);
	omega = 0.5 ;
	sigma = 0.5/27.21138386 ;
	sigmad = 1.0/sigma ;
	t0 = 2.0*41.34 + 10.0*41.34 ;
	e = e0*exp(-(t-t0)*(t-t0)*0.5/(sigmad*sigmad))*cos(omega*t) ;
	return e ;
}


void reada( double* v, double* x, int gridi, string fil){
	char file[3];
	ifstream rea(fil.c_str(), ios::in) ;
	for ( int i = 0 ; i < gridi ; i++ ){
		rea >> x[i] >> v[i] ;
	}
}

double linpolate(double xx, double* x, double* v, int grid){ //linear interpolation 
	double energy = 0.0 ;
	if ( xx >= x[0] && xx <= x[grid-1]){
		for ( int j = 0 ; j < grid ; j++ ){
			if ( xx >= x[j] && xx <= x[j+1] ){
				energy = (v[j+1]-v[j])*(xx-x[j])/(x[j+1]-x[j]) + v[j] ;
				return energy;
			}
		}
	}else{
		return energy;
	}
}

void pred(double dt, double net, complex<double>* psi, complex<double>* psi1, complex<double>* psi2, complex<double>* psi3, complex<double>* psi4){
	//net = nstates*p OR nstates*q1 OR nstates*
	double c1, c2 ,c3, c4;
	//PREDICTOR ROUTINE
	c1 = dt;
	c2 = c1*dt*0.5;
	c3 = c2*dt/3.0;
	c4 = c3*dt*0.25;
	for ( int i = 0; i < net; i++){
		psi[i] = psi[i]+ (c1*psi1[i]) + (c2*psi2[i]) + (c3*psi3[i]) + (c4*psi4[i]);
		psi1[i] = psi1[i]+ (c1*psi2[i]) + (c2*psi3[i]) + (c3*psi4[i]);
		psi2[i] = psi2[i]+ (c1*psi3[i]) + (c2*psi4[i]);
		psi3[i] = psi3[i]+ (c1*psi4[i]);
	}
/*	for ( int i = 0; i < net; i++){
		psi[i] = psi[i] + (c1*psi1[i]) + (c2*psi2[i]) ;
		psi1[i] = psi1[i] + (c1*psi2[i]) ;
	}*/
}


void corr(double dt, double net, complex<double>* psi, complex<double>* psi1, complex<double>* psi2, complex<double>* psi3, complex<double>* psi4, complex<double>* prodpsi){
	double gear0, gear1, gear2, gear3, gear4, cd0, cd1, cd2, cd3, cd4 ;
	double c1, c2 ,c3, c4;
	c1 = dt;
	c2 = c1*dt*0.5;
	c3 = c2*dt/3.0;
	c4 = c3*dt*0.25;
	//correction will be passed as correction 
	//CORRECTOR ROUTINE
/*	gear0=5.0/12.0;
	gear1=1.0;
	gear2=0.5;*/
	gear0 = 251.0/720.0;
	gear2 = 11.0/12.0;
	gear3 = 1.0/3.0;
	gear4 = 1.0/24.0;

	cd0 = gear0*dt;
	cd2 = gear2*c1/c2;
	cd3 = gear3*c1/c3;
	cd4 = gear4*c1/c4;

	complex<double> corr = 0.0;
	for ( int i = 0 ; i < net ; i++ ){
		corr = prodpsi[i] - psi1[i];
		psi[i] = psi[i] + corr*cd0;
		psi1[i] = prodpsi[i];
		psi2[i] = psi2[i] + cd2*corr;
		psi3[i] = psi3[i] + cd3*corr;
		psi4[i] = psi4[i] + cd4*corr;
	}
}

void kmatrixcalc(complex<double>** KE){
//Page 317, Introduction to Quantum Mechanics - A TD Perspective, DJ Tannor, 2007
//	5 point finite differece KE; the same matrix can be used to define Fourier method or other finite difference methods
//	complex<double> w = exp(complex<double>(0.0, 2.0*M_PI*ixgrid));
	double K = 1.0/dx;
	for ( int i = 0 ; i < xgrid ; i++ ){
		for ( int j = 0 ; j < xgrid ; j++ ){
			if ( i == j ){
				KE[i][j] = (30.0/12.0)*(0.5/mass)*K*K;
			}else if ( fabs(i-j) == 1 ){
				KE[i][j] = -(16.0/12.0)*(0.5/mass)*K*K;
			}else if ( fabs(i-j) == 2 ){
				KE[i][j] = (1.0/12.0)*(0.5/mass)*K*K;
			}
			imag(KE[i][j]) = 0.0;
		}
	}
}

double dipole(double x, int surf){
	// J. Chem. Phys. 115, 6438 (2001); https://doi.org/10.1063/1.1400139
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
	 	if(x>=2.05){ moment = -0.4*x+0.22 ;}
	}
	if ( surf == 6 ){ //ed
		if(x<2.1){ moment = -0.2 ;}
		if(x>=2.1){ moment = -0.29*x+0.4 ;}
	}
	return moment ;
}


double coupling(double x, int surf1, int surf2){
	//coupling between diabatic states
	// J. Chem. Phys. 115, 6438 (2001); https://doi.org/10.1063/1.1400139
	double v = 0.0 ;
	double xeq = 0.0, beta = 0.0, c = 0.0 ;
	if ( (surf1 == 1 && surf2 == 2) || (surf2 == 1 && surf1 == 2) ){//b-c
		xeq	= 1.71141 ;
		beta	= 2.21929 ;
		c	= -5.64067e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq)) ;
	}else if ((surf1 == 1 && surf2 == 3) || (surf2 == 1 && surf1 == 3)){//b-o
		xeq	= 1.82855 ;
		beta	= 3.65381 ;
		c	= 2.18468e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq)) ;
	}else if ((surf1 == 2 && surf2 == 3) || (surf2 == 2 && surf1 == 3)){//c-o
		xeq	= 1.26651 ;
		beta	= 1.90404 ;
		c	= 1.49604e-2 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq)) ;
	}else if ((surf1 == 4 && surf2 == 5) || (surf2 == 4 && surf1 == 5)){//bd-cd
		xeq	= 1.85220 ;
		beta	= 5.78569 ;
		c	= -1.31280e-2 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq)) ;
	}else if ((surf1 == 4 && surf2 == 6) || (surf2 == 4 && surf1 == 6)){//bd-ed
		xeq	= 1.85220 ;
		beta	= 5.78569 ;
		c	= -4.60098e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq)) ;
	}else if ((surf1 == 5 && surf2 == 6) || (surf2 == 5 && surf1 == 6)){//cd-ed
		xeq	= 2.02162 ; 
		beta	= 3.39115 ; 
		c	= 2.85983e-3 ;
		v = c*exp(-beta*(x-xeq)*(x-xeq)) ;
	}
	return v ;
}

double H_X ( double x){
	// The Physics of Fluids 4, 637 (1961); https://doi.org/10.1063/1.1706374 
	double xeq = 2.04158, beta = 1.50765, De = 0.363965;
	double V = (1.0-exp(-beta*(x-xeq)));
	V = V*V*De;
	return V;
}

double H_analytic ( double x, int surf){
	//J. Chem. Phys. 115, 6438 (2001); https://doi.org/10.1063/1.1400139
	double xeq, beta, alfa;
	double* c;
	double cdum[6][8] ={
		{0.458401, 0.0, 0.0206768, -0.0118135, 0.115635, -0.0842532, -0.00355449,0.0412016},
		{0.471233, 0.0, 0.319474, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.477476, 0.0, 0.316443, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.470701, 0.0, 0.0542505, 0.0219076, 0.0212176, 0.0111611, 0.00161939, 0.0},
		{0.472326, 0.0, 0.338606, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.522845, 0.0, 0.319334, 0.0, 0.0, 0.0, 0.0, 0.0},
	} ;
	switch (surf){
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
	r = (1.0-exp(-beta*(x-xeq)))/(1.0+(alfa*exp(-beta*(x-xeq)))) ;
	V = 0.0 ;
	for ( int i = 0 ; i < 8 ; i++ ) V=V+(c[i]*pow(r,i)) ;
	if ( surf == 1 && x <= 1.819 ){
		V = 0.0 ; 	//b LJ patch
		double depth = 0.0779 ;
		double assimptot = 0.536 ;
		double n = 8 ;
		double m = 1.8 ;
		double r_patch = 1.819 ;
		V = (depth/(n-m))*( m*pow(xeq/x,n) - n*pow(xeq/x,m) )+assimptot ;
	}
	if ( surf == 4 && x <= 1.7 ){
		V = 0.0 ; 	//bd LJ patch
		double depth = 0.1097 ;
		double assimptot = 0.5804 ;
		double n = 2.71 ;
		double m = 1.8 ;
		double r_patch = 1.7 ;
		V = (depth/(n-m))*( m*pow(xeq/x,n) - n*pow(xeq/x,m) )+assimptot ;
	}
	return V ;
}


complex<double> vimag(double x){
	double U0 = 0.2, alpha = 1.90;
	complex<double> V = complex<double>(0.0, 0.0);
	if ( x > 5.5 ){
		V = complex<double>(0.0,1.0)*U0/(cosh(alpha*(x-8.0))*cosh(alpha*(x-8.0)))- 0.61485409947493092e-4;
	}
//	V = complex<double>(0.0, 0.0);
	return V;
}
