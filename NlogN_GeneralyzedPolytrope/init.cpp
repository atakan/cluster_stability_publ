#include "init.h"

using namespace std;

ofstream input; // file consists of star parameters at the end of that program

int N;
float n,m; // secondary n & m values

vector<double> x; // vector holding theta values
vector<double> t; // vector holding eta values
vector<double> y; // vector holding differential values

vector<double> gx_tabulated; // vector holding tabulated g(x) values for specific n&m
vector<double> hx_tabulated; // vector holding tabulated g(x) values for specific n&m

static double f(double t, double x, double y){  
	return y;
}

static double g(double t, double x, double y){
	
	double val;
	/* this is problematic for x<0 */

	if (fabs(t)<1e-12) {
		val = 0.0;
	} else if (fabs(t)<1e-3) {
	//	val =  a[1]*(2*m+2)*pow(t, 2*m+1) +
	//	       a[2]*(4*m+4)*pow(t, 4*m+3) +
	//	       a[3]*(6*m+6)*pow(t, 6*m+5);
		val = -pow(t, 2*m)*pow(x, n+m)/3.0;
	} else {
		val = -2.0*y/t - pow(t, 2*m)*pow(x, n+m);
	}
	if(!finite(val)) {
		printf("a uuu (2)\n");
		abort();
	}
	return val;
}

double incr(double k1, double k2, double k3, double k4, double k5, double k6){
       	return (25.0/216.0*k1 + 1408.0/2565.0*k3
		+ 2197.0/4104.0*k4 - 1.0/5.0*k5);
}

long binary_search(double random,vector<double> mass)  // Binary search to find closest point
{
	
	long first = 0;
	long last = mass.size() - 1;
	while(fabs(first - last) > 1)
	{
		if(random <= mass[(first + last)/2])
			last = (first + last)/2;
		else
			first = (first + last)/2;	
	}
	
	double delta = fabs(random - mass[first]);
	if(delta > fabs(random - mass[last]))
		return last;
	return first;
	
}

double gx_function(double x)  // g(x) function
{
	
	double c = ( double (1.0) ) / (n - 0.5);
	return pow( (2 - pow(x,c)) , n - 1.5 ) * pow( (1 - pow(x,c)) , 2*m + 2 ) ;
	
} 

void simpson(double a,double b) // simpson integrator for g(x)
{
	int n = 10000; // number of intervals  ??
	double h = (double(b-a)) / n;  // width of each interval
	
	double y0 = gx_function(a); // y0
	double yn = gx_function(b); // yn
	
	double x = a + h;
	int i = 1;
	double even = 0;
	double odd = 0;
	
	while( i < n)
	{
		if( i % 2 == 0) //  i is even
			even += gx_function(x);
		else
			odd += gx_function(x);
		x += h;
		i++;
	}
	
	double result = (h/3)* (y0 + yn + 2 * even + 4 * odd);
	
	//cout << b << "   " << result << endl;
	gx_tabulated.push_back(result);	
}

double hx_function(double x)  // h(x) function
{
	double b = (double(1)) / (m+1); 
	return pow( (2 - pow(x,b)) , m ); 
	
} 

void simpson2(double a,double b) // simpson integrator for h(x)
{
	int n = 10000; // number of intervals  ??
	double h = (double(b-a)) / n;  // width of each interval
	
	double y0 = hx_function(a); // y0
	double yn = hx_function(b); // yn
	
	double x = a + h;
	int i = 1;
	double even = 0;
	double odd = 0;
	
	while( i < n)
	{
		if( i % 2 == 0) //  i is even
			even += hx_function(x);
		else
			odd += hx_function(x);
		x += h;
		i++;
	}
	
	double result = (h/3)* (y0 + yn + 2 * even + 4 * odd);
	
	//cout << b << "   " << result << endl;
	hx_tabulated.push_back(result);	
}

void calc(Star sun)  //generalyzed polytrope integrator(runge kutta included)
{
	
	double h, tol = 1e-10, hmax = 0.01, error;
	double k1, k2, k3, k4, k5, k6;
	double l1, l2, l3, l4, l5, l6;
	double t2, t3, t4, t5, t6;
	double x2, x3, x4, x5, x6;
	double y2, y3, y4, y5, y6;
	double r1s, E1s; /* starred quantities */
	vector<double> rs, Mrs, Us, rhos; /* starred quantities */	
	double eta_1, dtheta_deta_1;
	long i, array_incr = 10000, array_size=array_incr, Ni;

	/* initial conditions */
	
	/// Runge Kutta Part .... ///
	h = tol;
	i = 0;
	t.push_back(0.0);	/* eta */
	x.push_back(1.0);	/* theta */
	y.push_back(0.0);	/* dx/dt = d(theta) / d(eta) */
	// printf("%e %e %e\n", a[1], a[2], a[3]);
	do {
		//printf("%e %e %e\n", t[i], x[i], y[i]);
		k1 = h*f(t[i], x[i], y[i]);
		l1 = h*g(t[i], x[i], y[i]);

		t2 = t[i] + h/4.0;
		x2 = x[i] + k1/4.0;
		if (x2 < 0) {
			h *= 0.5;
			continue;
		}
		y2 = y[i] + l1/4.0;
		k2 = h*f(t2, x2, y2);
		l2 = h*g(t2, x2, y2);

		t3 = t[i] + 3.0/8.0*h; 
		x3 = x[i] + 3.0/32.0*k1 + 9.0/32.0*k2; 
		if (x3 < 0) {
			h *= 0.5;
			continue;
		}
		y3 = y[i] + 3.0/32.0*l1 + 9.0/32.0*l2;
		k3 = h*f(t3, x3, y3);
		l3 = h*g(t3, x3, y3);

		t4 = t[i] + 12.0/13.0*h;
		x4 = x[i] + 1932.0/2197*k1 - 7200.0/2197*k2 + 7296.0/2197*k3;
		if (x4 < 0) {
			h *= 0.5;
			continue;
		}
		y4 = y[i] + 1932.0/2197*l1 - 7200.0/2197*l2 + 7296.0/2197*l3;
		k4 = h*f(t4, x4, y4);
		l4 = h*g(t4, x4, y4);

		t5 = t[i] + h;
		x5 = x[i] + 439.0/216.0*k1 - 8.0*k2 
		          + 3680.0/513.0*k3 - 845.0/4104*k4;
		if (x5 < 0) {
			h *= 0.5;
			continue;
		}
		y5 = y[i] + 439.0/216.0*l1 - 8.0*l2 
		          + 3680.0/513.0*l3 - 845.0/4104*l4;
		k5 = h*f(t5, x5, y5);
		l5 = h*g(t5, x5, y5);

		t6 = t[i] + h/2.0;
		x6 = x[i] - 8.0/27.0*k1 + 2.0*k2 - 3544.0/2565.0*k3
		          + 1859.0/4104.0*k4 - 11.0/40.0*k5;
		if (x6 < 0) {
			h *= 0.5;
			continue;
		}
		y6 = y[i] - 8.0/27.0*l1 + 2.0*l2 - 3544.0/2565.0*l3
		          + 1859.0/4104.0*l4 - 11.0/40.0*l5;
		k6 = h*f(t6, x6, y6);
		l6 = h*g(t6, x6, y6);

		error = 1.0/360.0*k1 - 128.0/4275.0*k3 - 2197.0/75240.0*k4
			+ 1.0/50.0*k5 + 2.0/55.0*k6;
		if (error>tol){
			h *= 0.5;
			continue;
		}
		
		i++;  // increment counter
		if(i>=array_size){
			array_size += array_incr;
		}

		x.push_back(x[i-1] + incr(k1, k2, k3, k4, k5, k6));
		y.push_back(y[i-1] + incr(l1, l2, l3, l4, l5, l6));
		t.push_back(t[i-1] + h);
		if (error<tol/2.0 && h<hmax/2.0){
			h *= 2.0;
		}

	} while (x[i]>1e-100);
	if (i<3) {
		printf("m=%g n=%g\n", m, n);
		printf("a uuu\n");
		abort();
	} 

	/* if we overshot at this stage, we interpolate */
	t[i] = (x[i-1]*t[i] - x[i]*t[i-1])/(x[i-1]-x[i]);
	y[i] = (x[i-1]*y[i] - x[i]*y[i-1])/(x[i-1]-x[i]);
	x[i] = 0.0;
	//printf("%e %e %e\n", t[i], x[i], y[i]);
	Ni = i-1;

	eta_1 = t[Ni-1];
	dtheta_deta_1 = y[Ni-1];
	
	r1s = (4*m+6)/(3*m-n+5);
	E1s = -1.0/r1s;
	for(i=0; i<Ni; i++){
		rs.push_back(r1s * t[i]/eta_1) ;
		Mrs.push_back(t[i]*t[i]*y[i] / (eta_1*eta_1*dtheta_deta_1)) ;
		Us.push_back(E1s - x[i]/(-eta_1*dtheta_deta_1*r1s)) ;
		rhos.push_back(1.0/(4*PI) *eta_1 /(-dtheta_deta_1)
			/(r1s*r1s*r1s) *pow(t[i], 2*m) *pow(x[i], n+m)) ;
		//printf("%e %e %e %e\n", rs[i], Mrs[i], Us[i], rhos[i]);
	}

	// Runge Kutta ends.... //
	
	double a,b,delta;   // g(x) numeric integrator
	delta = 1e-4;
	a = 0;
	b = delta; // ??
	while(b <= 1 + 1e-10)
	{
		simpson(a,b);
		b+=delta; // ??
	}
	
	double c = ( double (1.0) ) / (n - 0.5);
	double b_value = (double(1)) / (m+1);   // m = -1/2 case ' ini eklemedik...
	
	
	delta = 1e-4;  // h(x) numeric integrator
	a = 0;
	b = delta; // ??
	while(b <= 1 + 1e-10)
	{
		simpson2(a,b);
		b+=delta; // ??
	}
	
	
	//***********************************************************************************//
	
	// find r,u,v //
	
	//***********************************************************************************//
	
	int quiet[10]; // for quite start implementation
	for(int i=0 ; i < 10 ; i++)
		quiet[i] = 0;
		
	double x1;	// random number
	int limit = N / 10;
	int ind;
	
	
	
	
	for(int i=0 ; i< N ; i++)  // for each star!
	{

		//x1 = ((double) rand() / (RAND_MAX)); // for Mr calculation
		// ---- Quiet Start -----//
		while(1)
		{

			x1 = ((double) rand() / (RAND_MAX)); // for Mr calculation
			ind = x1 * 10;
			if(ind == 10)
				ind = 9;
			if(quiet[ind] < limit)
			{
				quiet[ind] += 1;
				break;
			}
		}
		
		// ----  -----//	
	
		double r = rs[binary_search(x1,Mrs)];  // r value 
	
		/*cout << x1 << endl;
		cout << Mrs[binary_search(x1,Mrs)-1] << endl;
		cout << Mrs[binary_search(x1,Mrs)] << endl;
		cout << Mrs[binary_search(x1,Mrs)+1] << endl;*/
		
		double w1 = pow((2*E1s - 2*Us[binary_search(x1,Mrs)]) , 0.5); // maximal velocity at the given point ( w1 )
	
		double x_random = ((double) rand() / (RAND_MAX)); // for g(x) integration calculation
		double target_x = (double(binary_search( x_random * gx_tabulated[gx_tabulated.size()-1] , gx_tabulated))) * delta;  // x value !!!
		//cout << target_x << endl;
	
		double w,u,v;
		if( n == 0.5) // n=1/2 case
			w = w1; // w value !!
		else
			w = w1 * (1 - pow(target_x , c)); // w value !!

		
		x_random = ((double) rand() / (RAND_MAX)); // for h(x) integration calculation
		target_x = (double(binary_search( x_random * hx_tabulated[hx_tabulated.size()-1] , hx_tabulated))) * delta;  // x value !!!
		
		double cosmu = 1 - pow(target_x , b_value);
		u = w * cosmu;
		if((rand() % 2) == 0) // with 50% probability u will change its sign
			u *= -1;
	
		double sinmu = pow((1 - cosmu*cosmu) , 0.5);
		v = w * sinmu;
	
		/*cout << w << endl;
		cout << u << endl; 
		cout << v << endl;*/
		
		if(n == 0.5 && m == -0.5)  // Appendix 2 case added here !!
		{
			v = sqrt((double(1.5 - 1.125 * r)) / 2 );
			u = v;
			if((rand() % 2) == 0) // with 50% probability u will change its sign
				u *= -1;

		}
		//radius = r , rad_vel = u , tan_vel = v;  ////////////************************************///////////////////
		
		sun.radius = r;  // assign values
		sun.rad_vel = u;
		sun.tan_vel = v;
		sun.id = i;
		if(v == 0) // avoid from purely radial orbit
		{
			i -= 1;
			quiet[ind] -= 1; // to avoid infinite loop in quite start rand() calculation
			continue;
		}
		
		input << sun.id << "   " << sun.radius << "   " << sun.rad_vel << "   " << sun.tan_vel << endl << flush; // write parameters
	}
	
	return;
}


int main(int argc, char** argv)
{
	
	input.open("input.txt");

	N = atoi(argv[1]); // reading N (number of the stars) from the terminal
	n = atof(argv[2]); // n index
	m = atof(argv[3]); // m index
	input << N << "   " << n << "   " << m << endl << flush;

	Star sun; // construct star
	srand(time(NULL)); // initialize the current time
	calc(sun); // call initiator

	input.close(); // close file

}