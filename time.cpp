#include<complex>
#include<stdio.h>
//The regular pulse 
double normal(double t){
	double e, e0, omega, t0, sigma, sigmad;
	e0 = 0.0001;
	omega = 0.5;
	sigma = 1.0/27.211;
	t0 = 4.5/sigma;
	sigmad = 1.0/sigma;
	e = e0*exp(-(t-t0)*(t-t0)*0.5/(sigmad*sigmad))*cos(omega*t);
	return e;
}

//The intended chirp pulse
double chirp(double t){
	double e, e0, omega, t0, sigma;
	e0 = 0.01;
	omega = 0.5;
	sigma = 2.0/13.6;
	t0 = 10;
	e = e0*exp(-(t-t0)*(t-t0)/(2.0*sigma*sigma))*cos((omega + (10*t))*t);
	return e;
}

ofstream pop("pop.txt",ios::out);
// ofstream cpop("cpop.txt",ios::out);
//for Gear P-C routine
//no need to use class. I was just experimenting
class gear{
	void dynamics(double (*)( double), double** , double, double, double*, string*, int*, int*, double**);
	//dynamics calculates H(t)
	void cocalc(complex<double>* coeff, complex<double>* coeff1, complex<double>* coeff2, complex<double>* coeff3, complex<double>* coeff4, double** Hmatrix, double dt, complex<double>* prodvector, string* surface, double now);
	//P-C routine
public:
	int net;
	void roll(double** Hmatrix, double* phi, string* surface, int* Lang, int* statecount, double** dipole);
	//starts the function dynamics
};

void gear::roll(double** Hmatrix, double* phi, string* surface, int* Lang, int* statecount,double** dipole){
	double t,dt;
	t = 5000.0; //total time
	dt = 0.01; //time step
	double (*pfunc)(double) = normal; //kind of the pulse
	dynamics(pfunc,&Hmatrix[0],t,dt,&phi[0],&surface[0],&Lang[0],&statecount[0],&dipole[0]);
}

void gear::dynamics(double (*pulse)(double ), double** Hmatrix, double t, double dt, double* phi, string* surface, int* Lang, int* statecount, double** dipole){
	complex<double>* coeff = new complex<double>[net]; // c_0
	complex<double>* coeff1 = new complex<double>[net]; // c_1
	complex<double>* coeff2 = new complex<double>[net]; //c_2
	complex<double>* coeff3 = new complex<double>[net]; //c_3
	complex<double>* coeff4 = new complex<double>[net]; //c_4
	complex<double>* prodvector = new complex<double>[net]; // H*C matrix

	for (int i = 0; i<net; i ++ ){
		coeff[i] = (0.0,0.0);
		prodvector[i] = (0.0,0.0);
	}
	real(coeff[0]) = 1.0;
	for (int i = 0; i<net; i ++ ){
		for (int j = 0; j<net ; j ++ ){
			prodvector[i] = prodvector[i] + (Hmatrix[i][j]*coeff[j]*complex<double>(0.0,-1.0));
			coeff1[i] = prodvector[i];
		}
	}

	int row,col,rowh,colh;
	double H_dum;
	complex<double> cosum;
	cout << "# net is " << net << endl;
	for (int a = 0; a*dt <=t; a ++ ){
		for (int i = 0; i<35; i ++ ){
			for (int j = 35; j<net; j ++ ){
				int k = 0,l = 34;
				for (int i1 = 1; j-l>0; i1 ++ ){
					l = l + statecount[i1];
					k = i1;
				}
				Hmatrix[i][j] = 0.0;
				Hmatrix[j][i] = 0.0;
				H_dum = 0.0;
				if ((k != 0) && Lang[k]==1){
	/*				for (int m = 0; m<gridsize; m ++ ){
						H_dum = H_dum + (phi[(i*gridsize) + m]*moment(surface[k],1 + (m*0.0001))*phi[(j*gridsize) + m]);
					}*/
					H_dum = dipole[i][j];
	//				cout << "# surface is " << surface[k] << "; pulse " << pulse(a*dt) << endl;
					Hmatrix[i][j] = H_dum*pulse(a*dt);
					Hmatrix[j][i] = Hmatrix[i][j];
				}
				prodvector[j] = 0.0;
			}
			prodvector[i] = 0.0;
		}
		cocalc(&coeff[0],&coeff1[0],&coeff2[0],&coeff3[0],&coeff4[0],&Hmatrix[0],dt,&prodvector[0],&surface[0],a*dt);
		cosum = (0.0,0.0);
		for (int p = 0; p<net; p ++ ) cosum = cosum + (coeff[p]*coeff[p]);
		cout << "# now = " << a*dt << "; sum = " << cosum << "\t [" << abs(cosum) << "]" << endl;
	}

	delete[] coeff;
	delete[] coeff1;
	delete[] coeff2;
	delete[] coeff3;
	delete[] coeff4;
}

void gear::cocalc(complex<double>* coeff, complex<double>* coeff1, complex<double>* coeff2, complex<double>* coeff3, complex<double>* coeff4, double** Hmatrix, double dt, complex<double>* prodvector, string* surface, double now){
	pop.precision(12);
	double c1, c2, c3, c4, gear0, gear1, gear2, gear3, gear4, cd0, cd1, cd2, cd3, cd4;
	//Gear Predictor-coerrector routine for first order equation with four coefficients ie, upto third derivative of coeff
	//Coeffs from "Computer Simulation of Liquids" - M. P. Allen, D. J. Tildesley, Appendix E

	//PREDICTOR ROUTINE
	c1 = dt;
	c2 = c1*dt/2.0;
	c3 = c2*dt/3.0;
	c4 = c3*dt/4.0;
	for (int i = 0; i<net; i ++ ){
		coeff[i] = coeff[i] + (c1*coeff1[i]) + (c2*coeff2[i]) + (c3*coeff3[i]) + (c4*coeff4[i]);
		coeff1[i] = coeff1[i] + (c1*coeff2[i]) + (c2*coeff3[i]) + (c3*coeff4[i]);
		coeff2[i] = coeff2[i] + (c1*coeff3[i]) + (c2*coeff4[i]);
		coeff3[i] = coeff3[i] + (c1*coeff4[i]);
		prodvector[i] = (0.0,0.0);
	}

	// char cofile[3];

	//To print the Hamiltonian matrix
	/*sprintf(cofile,"dynamics/p_%2.2f",now);
	const char* ouitfile=cofile;
	ofstream of2(ouitfile,ios::out);
	for (int i = 0; i<net; i ++ ){
		for (int j = 0; j<net; j ++ ){
			of2 << Hmatrix[i][j] << "\t";
		}
		of2 << endl;
	}*/

	//CORRECTOR ROUTINE
	for (int i = 0; i<net; i ++ ){
		for (int j = 0; j<net ; j ++ ){
			prodvector[i] = prodvector[i] + (Hmatrix[i][j]*coeff[j]*complex<double>(0.0,-1.0));
		}
	}
	gear0 = 251.0/720.0;
	gear1 = 1.0;
	gear2 = 11.0/12.0;
	gear3 = 1.0/3.0;
	gear4 = 1.0/24.0;

	cd0=gear0*c1;
	cd2=gear2*c1/c2;
	cd3=gear3*c1/c3;
	cd4=gear4*c1/c4;

	// sprintf(cofile,"dynamics/%2.2f",now);
	//To calculate the population
	double X_cosum=0.0, b_cosum = 0.0, c_cosum = 0.0, o_cosum = 0.0, bd_cosum = 0.0, cd_cosum = 0.0, ed_cosum = 0.0;
	complex<double> corr;
	for (int i = 0; i<net; i ++ ){
		corr = prodvector[i]-coeff1[i];
		coeff[i] = coeff[i] + (cd0*corr);
		coeff1[i] = prodvector[i];
		coeff2[i] = coeff2[i] + (cd2*corr);
		coeff3[i] = coeff3[i] + (cd3*corr);
		coeff4[i] = coeff4[i] + (cd4*corr);
		prodvector[i]=0.0;
	//	if (2*now - 2*trunc(now) ==0){
		//	sprintf(cofile,"dynamics/%2.2f",now);
		//	const char* outfile=cofile;
		//	ofstream oft(outfile,ios::out|ios::app);
		//	oft << real(coeff[i]) << "\t" << imag(coeff[i]) << "\t" << endl;
	//	}
		if (i < 35){ //vibrational levels
			X_cosum = X_cosum + (real(coeff[i])*real(coeff[i])) + (imag(coeff[i])*imag(coeff[i]));
		} else if (i < 35 + 25 && i >= 35){
			b_cosum = b_cosum + (real(coeff[i])*real(coeff[i])) + (imag(coeff[i])*imag(coeff[i]));
		} else if ( i>= 35 + 25 && i < 35 + 25 + 35){
			c_cosum = c_cosum + (real(coeff[i])*real(coeff[i])) + (imag(coeff[i])*imag(coeff[i]));
		} else if ( i>= 35 + 25 + 35 && i < 35 + 25 + 35 + 35){
			o_cosum = o_cosum + (real(coeff[i])*real(coeff[i])) + (imag(coeff[i])*imag(coeff[i]));
		} else if ( i>= 35 + 25 + 35 + 35 && i < 35 + 25 + 35 + 35 + 25){
			bd_cosum = bd_cosum + (real(coeff[i])*real(coeff[i])) + (imag(coeff[i])*imag(coeff[i]));
		} else if ( i>= 35 + 25 + 35 + 35 + 25 && i < 35 + 25 + 35 + 35 + 25 + 35){
			cd_cosum = cd_cosum + (real(coeff[i])*real(coeff[i])) + (imag(coeff[i])*imag(coeff[i]));
		} else if ( i>= 35 + 25 + 35 + 35 + 25 + 35 && i < 35 + 25 + 35 + 35 + 25 + 35 + 35){
			ed_cosum = ed_cosum + (real(coeff[i])*real(coeff[i])) + (imag(coeff[i])*imag(coeff[i]));
		}
	}
	double sum = 0.0;
	sum = X_cosum + b_cosum + c_cosum + o_cosum + bd_cosum + cd_cosum + ed_cosum;
	if (2*now - 2*trunc(now) ==0){
		pop << now << "\t" << b_cosum << "\t" << c_cosum << "\t" << o_cosum << "\t" << bd_cosum << "\t" << cd_cosum << "\t" << ed_cosum << "\t" << sum << endl;
	//	cpop << now << "\t" << c_cosum << endl;
	}

}
