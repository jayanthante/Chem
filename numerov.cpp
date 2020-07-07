//Class under which the Numerov method is operated
// No need to use class; I was experiemtning
class numerov{
	double wave ( double x, double h, double phi, double phi_old, double Ea, double Eb, double Ec ) ;//phi ( i+1 ) is calculated using the method
	void normalize ( double* phi, double h ) ;//normalizes the wave function
	void decay ( double* l_phi, int t_begin, int t_end ) ;
public:
	int N ;//grid points
	void calc( double* hmatrix, int state, int &statecount, double* l_phi, double* eigen ) ;//the iterating process
} ;

void numerov::calc( double* hmatrix, int state, int &statecount, double* l_phi, double* eigen ) {
	double Ea, Eb, Ec ; //E (prev) , E (current) , E (next) 
	double Emin, Emax, Emax_dum, E, Emaxp ; //E (minimum) , E (maximum) , Emax (dummy) , E, Emax (previous) 
	double end, begin ; //classical end point (RHS) of the PEC , classical begining point (LHS) 
	double* phi = ( double* ) malloc ( N*sizeof ( double ) ) ; //defines phi starting from RHS
	double* V = ( double* ) malloc ( N*sizeof ( double ) ) ; //defining the gri for V
	double xmax, xmin, xeq ; //integration is performed between x ( max ) and x ( min ) 
	int nodecount = 0, l_nodecount = 0 ; //number of nodes of RHS and LHS wavefunctions
	xmax = 7.0 ;
	xmin = 1.6 ;
	double h = (xmax - xmin)/N ; //step of integration
	//double E_d = ( 2*12766.3254 ) ;
	double E_d = 2.0*12734.123 ; //units for 2*m ; hbar = 1
	int slopestat = 0, t_slope ;
	double xslope ;
	Emin = hmatrix[2]*E_d;
	for ( int i = 0 ; i < N ; i++ ) {
		V[i] = hmatrix[i]*E_d ;
		if ( V[i] <= Emin && xmin + i*h < 3.4 ){
			Emin = V[i];
			xeq = xmin + i*h;
		}
	}
	
	//for ( int i = N - 1 ; i >= 0 ; i-- ) {
	//	V[i] = hmatrix[state][state][i]*E_d ;
	//	if ( i < N - 2 ) {
	//		if ( xmin+ ((i + 5)*h ) < xeq ) {
	//			if ( slopestat == 0 && ( V[i+1]-V[i] ) /h > 1.e-2 ) {
	//				xslope = ( xmin+ ( i*h ) +xeq )*0.5- ( i*h*0.1 ) ;
	//				t_slope = ( xslope-xmin ) /h ;
	//				i = t_slope ;
	//				slopestat = 1 ;
	//				cout << "# change at " << xslope << endl ;
	//			}
	//			if ( slopestat == 1 ) {
	//				V[i] = V[i+1]- ( V[t_slope+1]-V[t_slope] )*1.15 ;
	//			}
	//		}
	//	}
	//}
	//routine for appropriate vmax
	int eqpos = (xeq-xmin)/h; //position of xeq in the grid
	double maxpos;
	Emax = Emin;
	for ( int i = eqpos - 1 ; i < N; i++ ){
		if ( V[i] > Emax ){
			Emax = V[i];
			maxpos = xmin + i*h;
		}
	}
	Emax = V[N-1];

	cout << "# Emax " << Emax/E_d << "( " << maxpos << "); Emin " << Emin/E_d << "( " << xeq << " ) " << endl ;
	if ( Emax < Emin ) {
		statecount = 0;
		return ;
	}
	Emax_dum = Emax ;
	phi[0] = 0.0 ;
	phi[1] = 1.e-50 ; //infinitesimal value
	for ( int i = 2 ; i < N ; i++ ) {
		phi[i] = wave ( xmin+i*h, h, phi[i-1], phi[i-2], Emax-V[i-1], Emax-V[i], Emax-V[i+1] ) ;
		if ( phi[i] == 0 ) phi[i] = 1.e-40 ; //in order to count the nodes, a small value is given for phi if it is zero
		if ( phi[i]*phi[i-1] < 0 ) 	nodecount = nodecount+1 ; //when the wavefunction changes sign, it's counted as a node
	}
	
	int maxnode;
	if ( nodecount > 2 ) maxnode = nodecount -1 ; //number of nodes at D_e
	if ( statecount > nodecount ) statecount = maxnode ;
	// if ( maxnode > 25 ) {
	// 	statecount = 35 ;
	// 	if ( surface == "bd" ) statecount = 25 ;
	// }else{
	// 	statecount = 25 ;
	// //	statecount = maxnode-1 ;
	// }
	cout << "# Maxnode is " << maxnode << " (" << nodecount << ") ; statecount is " << statecount << endl ;
	
	int stat = 0 ; //stat = 0 -- > unfinished ; stat = 1 -- > finished
	int k ;
	for ( int j = 0 ; j < 35 && j < statecount ; j++ ) { //procedure for 35 states including v = 0
		E = (Emax + Emin)*0.5 ;
		int iterate = 0 ;
		do{
			k = 0 ;
			do{
				k++ ;
			}while ( V[k] >= E ) ;
			begin = xmin + k*h ; //finds the LH classical limit at E
			int t_begin = k ;
			
			k = 0 ;
			do{
				k++ ;
			}while ( V[N-2-k] >= E ) ;
			end = xmax - k*h ; //finds the RH classical limit at E
			int t_end = (end - xmin)/h ;
		
			if ( t_end > N ) {
				free ( l_phi ) ;free ( phi ) ;
				exit ( 1 ) ;
				//just in case, if E is in continuum
			}
		
			nodecount = 0 ;
			l_nodecount = 0 ;
		
			l_phi[j*N] = 0. ;
			l_phi[ (j*N) +1] = 1.e-100 ;
			for ( int i = ( j*N ) +2 ; i < ( j+1 )*N ; i++ ) {
				l_phi[i] = wave ( xmin+ (i%N)*h, h, l_phi[i-1], l_phi[i-2], E - V[(i-1)%N], E-V[i%N], E-V[ ( i+1 ) %N] ) ;
				if ( l_phi[i]*l_phi[i-1] <= 0 && l_phi[i-1]*l_phi[i-2] != 0 ) 	l_nodecount = l_nodecount+1 ;
			}
			stat = 0 ;
			int nodemin = 0, nodemax = 0 ;
			phi[0] = 0. ;
			phi[1] = 1.e-100 ;
			for ( int i = 2 ; i < N ; i++ ) {
				phi[i] = wave ( xmin+ ( i*h ) , h, phi[i-1], phi[i-2], Emin-V[i-1], Emin-V[i], Emin-V[i+1] ) ;
				if ( phi[i]*phi[i-1] <= 0 && phi[i-1]*phi[i-2] != 0 ) 	nodemin = nodemin+1 ;
			}
			phi[0] = 0. ;
			phi[1] = 1e-100 ;
			for ( int i = 2 ; i < N ; i++ ) {
				phi[i] = wave ( xmin+i*h, h, phi[i-1], phi[i-2], Emax-V[i-1], Emax-V[i], Emax-V[i+1] ) ;
				if ( phi[i]*phi[i-1] <= 0 && phi[i-1]*phi[i-2] != 0 ) 	nodemax = nodemax+1 ;
			}
		//	cout << "# nodemin, nodemax " << nodemin << ", " << nodemax << endl ;
			iterate++ ;
		
		//	cout << "# L-nodes : " << l_nodecount << " ; end - " << end << " ; begin - " << begin << endl ;
		//	cout << "# Emax : " << Emax/E_d << " ; Emin " << Emin/E_d << " ; E " << E/E_d << endl ;
			if ( l_nodecount == j && fabs ( Emax-Emin ) /E_d <= 1e-14 ) {
/*				decay ( &l_phi[j*N], t_begin, t_end ) ;
				normalize ( &l_phi[j*N], h ) ;
				stat = 1 ;
		//		cout.precision ( 14 ) ;
		//		cout << "# Eigenvalue ( " << j << " ) = " << E/E_d << endl ;
				eigen[j] = E/E_d ;*/
				decay(&l_phi[j*N], t_begin,t_end);
				normalize(&l_phi[j*N],h);
				stat=1;
				eigen[j] = E/E_d;
				if ( state >= 7) {
					char file[3];
					sprintf(file,"%d_%d",state,j);
					const char* outfile = file;
					ofstream ofe(outfile,ios::out|ios::binary);
					ofe.precision(10);
					for (int i=j*N; i<(j+1)*N; i++){
						ofe << xmin+((i%N)*h) << "\t" << l_phi[i] << "\t" << phi[i%N] << endl;
					}
				}
			}else if ( fabs ( Emax-Emin ) /E_d > 1e-14 && l_nodecount == j ) {
				//trying to find the energy levels
				if ( nodemin <= j ) {
					if ( nodemax == j ) {
						Emin = E ;
						E = Emax ;
						Emax = Emaxp ;
					}else{
						Emaxp = Emax ;
						Emax = ( E+Emax )*0.5 ;
						Emin = E ;
						E = ( Emax+Emin )*0.5 ;
					}
				}else if ( nodemin > j ) {
					Emin = eigen[j-1]*E_d ;
				}
			}else if ( l_nodecount > j ) {
				if ( nodemin > j && j > 0 ) Emin = eigen[j-1]*E_d ;
				Emax = E ;
				E = (Emin + E)*0.5 ;
				stat = 0 ;
			}else if ( l_nodecount < j ) {
				Emin = E ;
				E = (Emax +E )*0.5 ;
				stat = 0 ;
			}
		//	cout << "# stat ; iterate ; j & l_nodecount " << stat << " " << iterate << " " << j << " " << l_nodecount << endl ;
		// while ( fabs ( Emax-Emin ) /E_d > 1e-12 && iterate < 1000 ) ;
		}while ( stat == 0 && iterate < 400 ) ;
		Emin = eigen[j]*E_d ;
		Emax = Emax_dum ;
	}
}

//Numerov Procedure: private
double numerov::wave ( double x, double h, double phi, double phi_old, double Ea, double Eb, double Ec ) {
	double psi, denom ;
	psi = ( 2.0 - ( 5.0*h*h*Eb/6.0 ) )*phi- ( 1.0+ (h*h*Ea/12.0))*phi_old ;
	denom = 1.0 + ( h*h*Ec/12.0 ) ;
	psi = psi/denom ;
	return psi ;
}

void numerov::normalize ( double* phi, double h ) {
	double norm = 0 ;
	for ( int i = 0 ; i < N ; i++ ) {
		norm = norm+ ( phi[i]*phi[i] ) ;
	}
	norm = sqrt ( norm ) ;
	// cout << "# norm is " << norm << endl ;
	//if ( fabs ( norm ) < 10e-5 ) 	for ( int i = 0 ; i < N ; i++ ) cout << i*h << "\t" << phi[i] << endl ;
	for ( int i = 0 ; i < N ; i++ ) {
		phi[i] = phi[i]/norm ;
	}
}

void numerov::decay ( double* phi, int t_begin, int t_end ) {
	int sign, mark, change ;
	double diff, logno ;
	//FOR phi
	sign = ( phi[t_end]-phi[t_end-1] ) /fabs ( phi[t_end]-phi[t_end-1] ) ; //find the sign of slope at t_end
	change = 0 ;
	logno = ( trunc ( log10 ( fabs ( phi[t_end] ) ) ) ) ; //a scaling factor of which's tenth power by which phi is divided so as to make the normalization easier
	// cout << "#logno is " << logno << endl ;
	for ( int i = 0 ; i < N ; i++ ) {
		if ( phi[i] == 0 ) phi[i] = 1.e-100 ;
		phi[i] = phi[i]/pow ( 10, logno ) ; //phi scaled
		diff = ( phi[i]-phi[i-1] ) ;
		/*The next if construct suppresses the explosion of the wf. It is done by considering whether the wavefunction or its slope change the sign.
		The integer change is considered so as to make the wf zero after the first occurence of any of the conditions described in the previous line.*/
		if ( i > t_end ) {
			if ( ( diff/fabs( diff ) != sign || phi[i]*phi[i-1] < 0 ) && change == 0 ) {
				mark = i ;
				change = 1 ;
			}
			if ( i > mark ) phi[i] = 0 ;
		}
	}
	
}
