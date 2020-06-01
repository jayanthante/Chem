#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <limits>
#include<math.h>
#include<algorithm>
#include<vector>
#include<string.h>
#include <fftw3.h>
#include<complex>

using namespace std ;
extern "C" {
void mmbsjn_( double *,  int *, double * ,  double *);
} // from Netlib imsl library


//==parameters==
const int xgrid = 1024 ; //grid size
const double ixgrid = 1.0/1024.0 ; //inverse of grid
const double mass = 12734.123 ; //mass of nitrogen
int ncheb = 0 ; // number of chebyshev polynomials
int T = 20000 ; // the duration of the simulation in au
double xmin = 1.6 ; //LHS of xgrid
double xmax = 7.0 ; // RHS of xgrid
double dx = (xmax-xmin)/(1023.0) ;//sqrt(2*M_PI/xgrid) ;

//==functions==
double pulse(double t) ;
double H_X ( double x);
double H_analytic ( double x, int surf);
double dipole(double x, int surf) ;
double coupling(double x, int surf1, int surf2) ;
double hamilpolate(double xx, double* x, double* v, int grid) ;
void reada( double* v, double* x, int gridi, string fil) ;
void initak(complex<double>* akx2);
void fft(complex<double>* in, complex<double>* out, int direc);
complex<double> vimag(double x);

//==main==
int main(){
	int nstates=11;

	complex<double>* akx2=new complex<double>[xgrid]() ;
	initak(&akx2[0]);
        ofstream ako("ak.dat",ios::out) ;
	for ( int i = 0 ; i < xgrid ; i++ ){
		akx2[i] = -0.5*akx2[i]/mass;
		ako << real(akx2[i]) << "\t" << imag(akx2[i]) << endl;
	} //Caluclates the momentum grid

	complex<double>* phi=new complex<double>[xgrid]() ; //phi is a single state wavefunction int the space grid
	complex<double>* kphi=new complex<double>[xgrid]() ; //kphi is the single state wavefunction in the momentum grid 

	vector< vector< complex<double> > > phistate( nstates, vector< complex <double> >(xgrid) ) ; //the whole wavefunction for all the states

	double x ;
	double vmin = 0.0, vmax = 1.50 ; //max and min potential energy in au

	double* x_X = new double[1024]() ;
	double* f_X = new double[1024]() ;
	reada( &f_X[0], &x_X[0], 1024, "X_0") ; //reads the ground state wf
	for ( int i = 0 ; i < xgrid ; i++ ){
		x = xmin + (i*dx) ;
		real(phi[i]) = hamilpolate( x, &x_X[0], &f_X[0], 1024) ;
		real(phistate[0][i]) = real(phi[i]);
	}
	delete[] x_X ;
	delete[] f_X ;

	double dE = (abs(akx2[xgrid/2])*abs(akx2[xgrid/2])/mass) + vmax - vmin ;
	double dt = 0.01;
	double H = dE*dt*0.5;
	complex<double> phase = exp(complex<double>(0.0,-H));
	double KE = 0.5*M_PI*M_PI/(mass*dx*dx) ; // kmax(vector)=PI/dx ; KE=k*k/2*m

	cout << "dt " << dt << " ; dE " << dE << " ; dx " << dx << " ; exp(dE) " << exp(complex<double>(0.0,dE)) << endl ;

	ncheb = int(H) + 12 ; //Calculates the number of Chebyshev polynomials required
	cout << "Ncheb " << ncheb << " ; phase " << phase << " ; xgrid " << xgrid << endl ;
	
	complex<double>* KEp= new complex<double>[xgrid] ; //
	complex<double>* chebc=new complex<double>[ncheb]() ; //Cehbyshev coefficients

	double* chebc_dum = new double[ncheb]();
	double ier;
	//Chebyshev coefficients = Chebychev coefficients * phase
	// H. TalEzer and R. Kosloff J.Chem.Phys.81, 3967(1984)
	// cout << "MMBSJN begin" << endl;
	mmbsjn_(&H, &ncheb, chebc_dum, &ier);
	// cout << "MMBSJN end" << endl;
	//MMBSJN_(double dE*dt*0.5, int ncheb, double* chebc , double ier) 
	cout << "Cheb" << endl ;
	chebc[0]=chebc_dum[0]*phase ;
	cout << chebc[0] << endl ;
	chebc[1]=chebc_dum[1]*2.0*complex<double>(0.0,1.0)*phase ;
	cout << chebc[1] << endl ;
	for (int i=2 ; i<ncheb ; i++){
		chebc[i]=chebc_dum[i]*2.0*pow(complex<double>(0.0,1.0),i)*phase ;
		cout << chebc[i] << endl ;
	}
	delete[] chebc_dum;
	double field = 0.0 ;
	double norm = 0.0 ; //normalization constant for the total wavefunction
	cout << "KE  = " << KE << " ; akx2(max) " <<  akx2[xgrid/2] << endl;

	static double hmatrix[11][11][1024];
	for ( int i = 0 ; i < xgrid ; i++ ) {
		for ( int is = 0 ; is < nstates ; is++ ) {
			for ( int j = 0 ; j < nstates ; j++ ) {
				hmatrix[is][j][i] = 0.0; 
			}
		}
	}

	for ( int i = 0 ; i < xgrid ; i++ ) {
		x = xmin + (i*dx) ;
		for ( int is = 0 ; is < 7 ; is++ ) {
			for ( int j = 0 ; j < 7 ; j++ ) {
				if ( is != j){
					hmatrix[is][j][i] = coupling(x,is,j) ; //calculates the couplings between the states is and j
				}
			}
		}
	}
	

	for ( int i = 0 ; i < xgrid ; i++ ) {
		x = xmin + (i*dx) ;
		for ( int is = 0 ; is < 7 ; is++ ) {
				hmatrix[is][is][i] = H_analytic(x,is); // calculates the potential energy part for state 'is' at distance x
			}
		}
	}

//	cout << "Hamiltonian done." << endl;

	//absorption potential
	complex<double>* vabs= new complex<double>[xgrid]; 
	for ( int i = 0 ; i < xgrid ; i++ ){
		vabs[i] = vimag( xmin + (i*dx) );
		cout << (xmin + (i*dx)) << "\t" << real(vabs[i]) << "\t" << imag(vabs[i]) << endl;
	}

    ofstream popu("pop.dat",ios::out) ; //to print the population

	double tym;
	int a_dt,corr_tym;
	complex<double> energia ; // = V*phi
	complex<double>** phicheb= new complex<double>*[ncheb] ;
	
	vector< vector< complex<double> > > phistate_new( nstates, vector< complex <double> >(xgrid) ) ;
	for ( int j = 0 ; j < ncheb ; j++ ) phicheb[j]= new complex<double>[xgrid] ;
	double* population = new double[nstates];
	double popdum;

	//TIME loop starts
	for ( int a = 0 ; a*dt < T ; a++ ){
		if ( a == 0 ) cout << "Propagator begins" << endl;
		tym = a*dt;
		a_dt = int(tym) ;
		field = pulse(a*dt) ;
		if ( a%100 == 0 ) cout << tym << "\t" << field << endl; //printing time at every 100 time steps

		for ( int is = 0 ; is < nstates ; is++ ){
			//memcpy (destination, source, size)
			for ( int i = 0 ; i < xgrid ; i++ ) phi[i] = phistate[is][i];
//			memcpy(phi,&phistate[is][0],xgrid*sizeof(phi));
//			memcpy(&phicheb[0][0],phi,xgrid*sizeof(phi));

//			calculating the momentum a/c to Fourier method by Kosloff & Kosloff
			fft(&phi[0], &kphi[0], -1);
			for ( int i = 0 ; i < xgrid ; i++ ){
				kphi[i] = akx2[i]*kphi[i]*ixgrid; 
			}
			fft(&kphi[0],&KEp[0],1);
			
			for ( int i = 0 ; i < xgrid ; i++ ){
				x = xmin + (i*dx) ;
				energia = complex<double>(0.0,0.0); //calculates H*PSI
				energia = hmatrix[is][is][i]*phistate[is][i] ;
				energia += vabs[i]*phistate[is][i] ;
				for ( int k = 0 ; k < nstates ; k++ ){
					//V*phi_is
					if ( is != k ) energia += hmatrix[is][k][i]*phistate[k][i] ;
					// -mu*dipole*phi_is
					if ( is == 0 && k != 0 ){
						energia += -field*dipole(x,k)*phistate[k][i] ;
					} else if ( k == 0 && is != 0 ){
						energia += -field*dipole(x,is)*phistate[k][i];
					}
				}
				//0th order Chebyshev polynomial
				phicheb[0][i] = phistate[is][i];
				//1st order Chebyshev polynomial
				phicheb[1][i] = (2.0*(KEp[i] + energia)/dE)-phistate[is][i] ;
			}
			//Higher order Chebyshev poynomials
			for ( int j = 2 ; j < ncheb ; j++ ){
//				memcpy(phi, &phicheb[j-1], xgrid*sizeof(phi));
//				Refer to Chebychev propagator by Kosloff et al
				for ( int i = 0 ; i < xgrid ; i++ ) phi[i] = phicheb[j-1][i];
				fft(&phi[0], &kphi[0], -1);
				for ( int i = 0 ; i < xgrid ; i++ ){
					kphi[i] = akx2[i]*kphi[i]*ixgrid; 
				}
				fft(&kphi[0],&KEp[0],1) ;
				for ( int i = 0 ; i < xgrid ; i++ ){
				 	x = xmin + (i*dx) ;
					energia = complex<double>(0.0,0.0); //to calculate the energy at x
					energia = hmatrix[is][is][i]*phi[i];
					energia += vabs[i]*phi[i];
					for ( int k = 0 ; k < nstates ; k++ ){
						if ( is != k ) energia += hmatrix[is][k][i]*phistate[k][i] ;
						if ( is == 0  && k!= 0 ){
							energia += field*dipole(x,k)*phistate[k][i] ;
						} else if ( k == 0 && is != 0 ){
							energia += field*dipole(x,is)*phistate[0][i];
						}
					}
					phicheb[j][i] = (2.0*((2.0*(KEp[i] + energia)/dE) - phicheb[j-1][i])) - phicheb[j-2][i]; //recursion
				}
			}
			for ( int i = 0 ; i < xgrid ; i++){
				for ( int j = 0 ; j < ncheb ; j++ ){
					phistate_new[is][i] += chebc[j]*phicheb[j][i] ;
				}
			}
		}//state loop ends
		phistate = phistate_new;
		for ( int is = 0 ; is < nstates ; is++ ){
			for ( int j = 0 ; j < xgrid ; j++){
				phistate_new[is][j] = complex<double>(0.0,0.0) ;
			}
		}
		norm = 0.0 ;
		

		if ( a%100 == 0 ){
			for ( int i = 0 ; i < nstates ; i++ ){
				popdum = 0.0;
				for ( int j = 0 ; j < xgrid ; j++ ){
					popdum += real(phistate[i][j])*real(phistate[i][j]) + imag(phistate[i][j])*imag(phistate[i][j]) ;
				}
				norm += popdum;
				population[i] = popdum;
			}
			popu << tym;
			double norm2 = 1.0/sqrt(norm);
			for ( int i = 0 ; i < nstates ; i++ ) {
				population[i] = population[i]/norm;
				popu << "\t" << population[i];
				for ( int j = 0 ; j < xgrid ; j++){
					phistate[i][j] = phistate[i][j]*norm2 ;
				}
			}
			popu << "\t" << norm;
			popu << endl;
		}
		delete[] population;
	} // time loop ends
	delete[] vabs ;
	delete[] phi;
	delete[] kphi;
	delete[] KEp;
	delete[] chebc;

	for ( int j = 0 ; j < ncheb ; j++ ) delete[] phicheb[j];
	delete[] phicheb;
	for ( int i = nstates-1 ; i >= 0 ; i-- ) phistate[i].clear();
	for ( int i = nstates-1 ; i >= 0 ; i-- ) phistate_new[i].clear();
}
//==================================


double pulse(double t){
	double e, e0, omega, t0, sigma, sigmad ;
	// everything to be given in au
	e0 = 0.01 ;
	omega = 0.5 ;
	sigma = 4.0/27.21138386 ;
	sigmad = 1.0/sigma ;
	t0 = 7.0*sigmad ;
	e = e0*exp(-(t-t0)*(t-t0)*0.5/(sigmad*sigmad))*cos(omega*t) ;
	return e ;
}

//FFT - FFTW
void fft(complex<double>* in, complex<double>* out, int direc){
	fftw_plan plan;
	plan = fftw_plan_dft_1d(xgrid, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), direc, FFTW_ESTIMATE);
	fftw_execute(plan); /* repeat as needed */
	fftw_destroy_plan(plan);
	fftw_cleanup();
}



double dipole(double x, int surf){

	return moment ;
}

double coupling(double x, int surf1, int surf2){

	return v;
}


void reada( double* v, double* x, int gridi, string fil){
	char file[3];
	ifstream rea(fil.c_str(), ios::in) ;
	for ( int i = 0 ; i < gridi ; i++ ){
		rea >> x[i] >> v[i] ;
	}
}

double hamilpolate(double xx, double* x, double* v, int grid){
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

void initak(complex<double>* akx2){
	ofstream ak("ak.dat",ios::out);
	double anorm = 2.0*M_PI/(1024.0*dx);
//	ak << anorm << endl;
	double c = -1.0;
	int sign = c*c ;
	akx2[0] = complex<double>(0.0, 0.0);
	double temp;
	for ( int i = 1 ; i <= 512 ; i++ ){
		temp = i*i*anorm*anorm*c;
		akx2[i] = temp;
		akx2[1024-i] = sign*temp;
	}
	ak.precision(5);
	for ( int i = 0 ; i < 1024 ; i++ ) ak << i << "\t" << real(akx2[i]) << "\t" << imag(akx2[i]) << endl;
}

complex<double> vimag(double x){
	double U0 = 0.2, alpha = 1.90;
	complex<double> V = complex<double>(0.0, 0.0);
	if ( x > 4.5 ){
		V = complex<double>(0.0,1.0)*U0/(cosh(alpha* (x-7.0))*cosh(alpha* (x-7.0)));
	}
//	V = complex<double>(0.0, 0.0);
	return  V;
}
