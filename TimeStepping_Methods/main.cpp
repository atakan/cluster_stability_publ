#include "check.cpp"

using namespace std;

int N;
double clck;
// G=1 , M=1 , R=1
float n,m; // essential parameters of system
double dt; // time-step

vector<double> ru_radii;  // ru difference testing
vector<double> ru_radial;

double first_radiis [9]; // for radiis instability
bool first = true; // check radiis
double difference = 0; // to get maximum difference

vector<Star> prev_cluster; // vector holding stars for previous time step
vector<Star> new_cluster; // vector holding stars for last state

double crossing_time = sqrt(2); // crossing time 

bool help;

void quicksort(vector<Star>& array,int left,int right) // quicksort
{
	int i,j;
	int temp;
	double pivot;
	i = left;
	j = right;
	pivot = array[(left + right) / 2].radius; // pivot value
	
	while(i <= j)
	{
		while(array[i].radius < pivot && i < right)
			i++;
		while(array[j].radius > pivot && j > left)
			j--;
	
		if(i <= j) // swap
		{
			Star p;
			p.radius = array[i].radius;
			p.rad_vel = array[i].rad_vel;
			p.tan_vel = array[i].tan_vel;
			p.id = array[i].id;

			array[i].radius = array[j].radius;
			array[i].rad_vel = array[j].rad_vel;
			array[i].tan_vel = array[j].tan_vel;
			array[i].id = array[j].id;

			array[j].radius = p.radius;
			array[j].rad_vel = p.rad_vel;
			array[j].tan_vel = p.tan_vel;
			array[j].id = p.id;

			i++;
			j--; 
		}
	}
	if(j > left)
		quicksort(array,left,j);
	if(i < right)
		quicksort(array,i,right);
	
}

int binary_search(double r, vector<Star>& array) // binary search
{
	int result;

	int start = 0;
	int finish = N - 1;
	int mid = (start + finish)/2;
	while(!(mid == start || mid == finish))
	{
		if((array[mid]).radius >= r)
			finish = mid;
		else if((array[mid]).radius <= r)
			start = mid;
		mid = (start + finish)/2;
	}
	if(r <= array[mid].radius)
		return mid;
	else
	{
		if(array[mid+1].radius >= r)
			return mid+1;
		else
			return mid+2;
	}
	return mid;

}

double iterate_single(double duration)
{
	double period = 0;
	
	while(period < duration && prev_cluster[0].radius < prev_cluster[1].radius)
	{
		double criter = sqrt(pow((prev_cluster[0]).radius , 3) / (0.5 / double(N)) ) * 1e-3 ;
		period += criter;
		if(period > duration)
		{
			period -= criter;
			break;
		}

		double local_radius = (prev_cluster[0]).radius ; // last radius of star
		double local_radvel = (prev_cluster[0]).rad_vel ; // last tangential velocity of star
		double local_tanvel = (prev_cluster[0]).tan_vel ; // last radial velocity of star

		double mass_inside = 0.5 / double(N)  ; // inside mass
		double a = mass_inside / (pow(local_radius,2)); // radial acceleration
		double l = local_tanvel * local_radius ; // l value

		(prev_cluster[0]).radius += criter * (local_radvel) + 0.5 * pow(criter,2) * ((-1)*a + double(pow(l,2)) / double(pow(local_radius,3))); // new radius
		(prev_cluster[0]).tan_vel = (local_tanvel * local_radius) / (prev_cluster[0]).radius;  // new tangential velocity
		(prev_cluster[0]).rad_vel += (criter) * ( ((-1)*a + double(pow(l,2)) / double(pow(local_radius,3))) ); // new radial velocity

	}
	return period;
}

void update() // update function for each time step
{
	int zombie;
	for (int i = 0; i< N; i++)
	{
		
		if(i == 0 && help)
		{
			zombie = prev_cluster[0].id;
			continue;
		}
			

		bool required = false;

		bool neg = false;
		double temp;
		double local_step = dt; // local time step of star !
		double duration = 0; // duration will be global time step --> dt
		
		double new_radius,new_radial,new_tangential;
		
		double local_radius = (prev_cluster[i]).radius ; // last radius of star
		double local_radvel = (prev_cluster[i]).rad_vel ; // last tangential velocity of star
		double local_tanvel = (prev_cluster[i]).tan_vel ; // last radial velocity of star

		new_radius = local_radius;
		new_radial = local_radvel;
		new_tangential = local_tanvel;

		double remain = dt;

		while(duration < dt)
		{
			temp = 0.5; // (1/2)*(own mass) added
			temp += i;
			
			double mass_inside = temp / N  ; // inside mass

			//cout << "Dynamic  " << mass_inside << endl;
		
			// mass inside could be measured from distribution function, statically, not dynamically.
		
			//mass_inside = (pow(local_radius,3) * sqrt(pow(local_radius,2) + 1)) / pow((pow(local_radius,2) + 1),2);  // from distribution function
			//cout << "Static  " << mass_inside << endl;
		
			double a = mass_inside / (pow(local_radius,2)); // radial acceleration
			double l = local_tanvel * local_radius ; // l value

			/*
			new_radial = local_radvel + local_step * ((-1)*a + double(pow(l,2)) / double(pow(local_radius,3))); // new radial velocity
			new_radius = local_radius + local_step * ((local_radvel + new_radial)/2); // new radius
			new_tangential = (local_tanvel * local_radius) / new_radius;  // new tangential velocity
			*/

			//  LEAPFROG  //
						
			new_radius = local_radius + local_step * (local_radvel) + 0.5 * pow(local_step,2) * ((-1)*a + double(pow(l,2)) / double(pow(local_radius,3))); // new radius
			new_tangential = (local_tanvel * local_radius) / new_radius;  // new tangential velocity

			/*
			temp = 0.5;
			temp += binary_search(new_radius, prev_cluster);
			if(new_radius > local_radius)  // eliminate itself from binary search result
				temp -= 1.0;

			mass_inside = temp / N  ; // new inside mass

			double a_new = mass_inside / (pow(new_radius,2)); // radial acceleration
			double l_new = new_tangential * new_radius ; // l value
			*/

			new_radial = local_radvel + (local_step  / 2.0) * ( ((-1)*a + double(pow(l,2)) / double(pow(local_radius,3))) ); // new radial velocity
				
			
			if(new_radius < 0 )  // if radius negative, local time step will be decreased
			{
				
				cout << "negative" << endl;

				neg = true;
				error <<  "err  " << prev_cluster[i].id << " " << local_radius << "  " << local_radvel << " " << new_radial << " negative "  <<  endl;
			
				local_step /= 10; // local time step  update
				new_radius = local_radius;
				//cout << local_step << endl;
												
			}
			else
			{
				
				duration += local_step; // duration update
				
				/*
				remain = dt - duration;
				local_step = remain;
				
				if(neg)
					error << " err  " << prev_cluster[i].id << " " << local_radius << "  " << local_radvel << " " << new_radial << endl;					
		
				neg = false;	
				
				//cout << duration << endl;
				//cout << "negative protection" << endl;
				local_radius = new_radius;
				local_radvel = new_radial;
				local_tanvel = new_tangential;
				required = true;
				*/
				
			}
		
		}
		
		(new_cluster[i]).update(new_radius,new_radial,new_tangential); // update the star on new state
		
	}
	
	for (int i = 0; i< N; i++)  // assign prev to new
	{
		if(i == 0 && help)
			continue;
		(prev_cluster[i]).radius = (new_cluster[i]).radius;
		(prev_cluster[i]).rad_vel = (new_cluster[i]).rad_vel;
		(prev_cluster[i]).tan_vel = (new_cluster[i]).tan_vel;
	}

	quicksort(prev_cluster,0,N-1);
	//quicksort(new_cluster,0,N-1);

	for (int i = 0; i< N; i++)
	{

		if(help && prev_cluster[i].id == zombie)
			continue;
		double local_step = dt; // local time step of star !
		double temp = 0.5;
		temp += i;
		double mass_inside = temp / N  ;

		double a_new = mass_inside / (pow((prev_cluster[i]).radius,2));
		double l_new = (prev_cluster[i]).tan_vel * (prev_cluster[i]).radius; 
		(prev_cluster[i]).rad_vel += (local_step  / 2.0) *  ((-1)*a_new + double(pow(l_new,2)) / double(pow((prev_cluster[i]).radius,3)));
	}
	
}

void initialize_system()  // create stars initially by using Henon gen poly. initialization algorithm
{
	Star sun; // construct star
	
	ifstream input;
	input.open("input.txt"); // read polytrope parameters
	input >> N >> n >> m;

	for(int i=0;i<N;i++)
	{
		input >> sun.id >> sun.radius >> sun.rad_vel >> sun.tan_vel; // read star parameters
		prev_cluster.push_back(sun); // fill cluster vectors
		new_cluster.push_back(sun);
	}
	input.close();

	quicksort(prev_cluster,0,N-1);  // Sorting --> (NLogN Goal)
	//quicksort(new_cluster,0,N-1);

	for(int i=0;i<N;i++)
	{
		cout << prev_cluster[i].id << "   " << prev_cluster[i].radius << "   " << prev_cluster[i].rad_vel << "   " << prev_cluster[i].tan_vel << endl;
	}
	
}

/*************************************************************************/

double timestep_criterion(vector<Star> & stars)
{
	double time_step = 1;
	double helper_time_step = 10;
	int id;
	for(int i=0;i<stars.size();i++)
	{
		double mass_inside = (0.5 + double(i)) / N;
		double criter = sqrt(pow((stars[i]).radius , 3) / mass_inside) ;
		/*
		if(i == 0)
			cout << "first  :  " << criter  << "   ,   ";
		if(i == 1)
			cout << "second  :  " << criter  << "   ,   ";
		*/
		if(criter < time_step)
		{
			helper_time_step = time_step;
			time_step = criter;
			id = i;
		}
		else if(criter <= helper_time_step)
		{
			helper_time_step = criter;
		}
			
	}
	//cout << id << "  " << time_step << "  " << helper_time_step << endl;
	if(id == 0 && time_step <  helper_time_step / 4)
	{
		help = true;
		return iterate_single(helper_time_step * 1e-3);
	}
		
	else
	{
		help = false;
		return time_step * 1e-3;
	}
}


int main(int argc, char** argv)    //////////// MAIN //////////////
{
	clck = 0; // clock
	
	int period,r_u_period; // periods
	
	dt = crossing_time/5000; // time step ??
	
	initialize_system();  // Creating stars by using Henon's initial method (Plummer model)
	open_files(); // open data files

	energychecking(prev_cluster); // check virial theorem
	obtain(prev_cluster); // shells' radius
	rminmaxcontrol(prev_cluster);  // rmin-rmax conservation checking
	r_u_control(prev_cluster); // check r-u state of the system
	
	int counter = 1;
  	period = 10;  // period value
	
	int r_u_counter = 1;
	r_u_period = (crossing_time / dt)*1.5; // period for checking ru plane instability
	
	clock_t start;
  	start=clock(); // fire time !
  	double duration;

  	//cout << dt << endl;

	while(clck < 30 * crossing_time) // main loop ( 30 Tcros.)
	{
		dt = timestep_criterion(prev_cluster);
		update(); // update critical values of stars
		miss(prev_cluster);

		clck += dt; // increment t by dt
		//cout << clck << endl;
		
		if(counter == period)
		{
			obtain(prev_cluster); // shells' radius
			rminmaxcontrol(prev_cluster); // rmin-rmax conservation checking
			//energychecking(prev_cluster);
			counter = 1;
		}
		else
			counter++;
			
		if(r_u_counter == r_u_period)
		{
			r_u_control(prev_cluster); // for r_u plane
			energychecking(prev_cluster);
			r_u_counter = 1;
			cout << " t = " << (clck / crossing_time) << " Tcros." << endl;
			duration = ( clock() - start ) / (double) CLOCKS_PER_SEC; // runtime of program
			cout << "running time = " << duration << " secs." <<  endl;
		}
		else
			r_u_counter++;
		
	}
	energychecking(prev_cluster);

	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC; // runtime of program
	cout << endl;
	cout << "running time = " << duration << " secs." <<  endl;
	cout << endl;
	
	//instability_test(); // automatic determination of instability in ru plane (NOT USED !)

	close_files(); // close data files

	return 0;
}
