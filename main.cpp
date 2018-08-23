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
vector<Star *> indices;

double crossing_time = sqrt(2); // crossing time 

void quicksort(vector<Star>& array,int left,int right) // indirect quicksort
{
	int i,j;
	int temp;
	double pivot;
	i = left;
	j = right;
	pivot = (*(indices[(left + right) / 2])).radius; // pivot value
	
	while(i <= j)
	{
		while( (*(indices[i])).radius < pivot && i < right)
			i++;
		while( (*(indices[j])).radius > pivot && j > left)
			j--;
	
		if(i <= j) // swap
		{
			Star *t = indices[i];
			//cout << t << "  " << indices[i] << endl;
			indices[i] = indices[j];
			indices[j] = t;

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

void update() // update function for each time step
{
	
	for (int i = 0; i< N; i++)
	{
		int ind = (*(indices[i])).id;

		bool required = false;

		bool neg = false;
		double temp;
		double local_step = dt; // local time step of star !
		double duration = 0; // duration will be global time step --> dt
		
		double new_radius,new_radial,new_tangential;
		
		double local_radius = (prev_cluster[ind]).radius ; // last radius of star
		double local_radvel = (prev_cluster[ind]).rad_vel ; // last tangential velocity of star
		double local_tanvel = (prev_cluster[ind]).tan_vel ; // last radial velocity of star

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
		
		(new_cluster[ind]).update(new_radius,new_radial,new_tangential); // update the star on new state
		
	}
	
	for (int i = 0; i< N; i++)  // assign prev to new
	{
		(prev_cluster[i]).radius = (new_cluster[i]).radius;
		(prev_cluster[i]).rad_vel = (new_cluster[i]).rad_vel;
		(prev_cluster[i]).tan_vel = (new_cluster[i]).tan_vel;
	}

	quicksort(prev_cluster,0,N-1);
	//quicksort(new_cluster,0,N-1);

	for (int i = 0; i< N; i++)
	{
		int ind = (*(indices[i])).id;

		double local_step = dt; // local time step of star !
		double temp = 0.5;
		temp += i;
		double mass_inside = temp / N  ;

		double a_new = mass_inside / (pow((prev_cluster[ind]).radius,2));
		double l_new = (prev_cluster[ind]).tan_vel * (prev_cluster[ind]).radius; 
		(prev_cluster[ind]).rad_vel += (local_step  / 2.0) *  ((-1)*a_new + double(pow(l_new,2)) / double(pow((prev_cluster[ind]).radius,3)));
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

	for(int i=0;i<N;i++)
	{
		indices.push_back( &(prev_cluster[i]) );
		//cout << (*(indices[i])).id << "   " << (*(indices[i])).radius << "   " << (*(indices[i])).rad_vel << "   " << (*(indices[i])).tan_vel << endl;
	}

	quicksort(prev_cluster,0,N-1);  // Sorting --> (NLogN Goal)
	//quicksort(new_cluster,0,N-1);

	/*
	for(int i=0;i<N;i++)
	{
		cout << prev_cluster[i].id << "   " << prev_cluster[i].radius << "   " << prev_cluster[i].rad_vel << "   " << prev_cluster[i].tan_vel << endl;
	}
	*/
	
	for(int i=0;i<N;i++)
	{
		cout << (*(indices[i])).id << "   " << (*(indices[i])).radius << "   " << (*(indices[i])).rad_vel << "   " << (*(indices[i])).tan_vel << endl;
	}
	
}

/*************************************************************************/

double timestep_criterion(vector<Star> & stars)
{
	double time_step = 1;
	for(int i=0;i<stars.size();i++)
	{
		int ind = (*(indices[i])).id;
		double mass_inside = (0.5 + double(i)) / N;
		double criter = sqrt(pow((stars[ind]).radius , 3) / mass_inside) ;
		/*
		if(i == 0)
			cout << "first  :  " << criter  << "   ,   ";
		if(i == 1)
			cout << "second  :  " << criter  << "   ,   ";
		*/
		if(criter < time_step)
			time_step = criter;
	}
	//cout << "  result :  " << time_step << endl;
	return time_step * 1e-3;
}


int main(int argc, char** argv)    //////////// MAIN //////////////
{
	clck = 0; // clock
	
	int period,r_u_period; // periods
	
	dt = crossing_time/5000; // time step ??
	
	initialize_system();  // Creating stars by using Henon's initial method (Plummer model)
	open_files(); // open data files

	energychecking(prev_cluster, indices); // check virial theorem
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
		//miss(prev_cluster);

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
			energychecking(prev_cluster,indices);
			r_u_counter = 1;
			cout << " t = " << (clck / crossing_time) << " Tcros." << endl;
			duration = ( clock() - start ) / (double) CLOCKS_PER_SEC; // runtime of program
			cout << "running time = " << duration << " secs." <<  endl;
		}
		else
			r_u_counter++;
		
	}
	energychecking(prev_cluster, indices);

	duration = ( clock() - start ) / (double) CLOCKS_PER_SEC; // runtime of program
	cout << endl;
	cout << "running time = " << duration << " secs." <<  endl;
	cout << endl;
	
	//instability_test(); // automatic determination of instability in ru plane (NOT USED !)

	close_files(); // close data files

	return 0;
}
