#include "init.h"

using namespace std;

ofstream output; // radiis.txt
ofstream minmax; // rmin-rmax controlling
ofstream r_u;  // r_u plane
ofstream error;  // to follow errors

ofstream script_output;


void obtain(vector<Star> & cluster)  // obtain Mass inside radiis 
{
	double scaling = 3 * PI / 16;
	
	double diff = 0;
	for(int i =N/10; i<N; i += N/10) // print shells' radius
	{
		output << cluster[i].radius << endl;  // not scaled
		if(first)
		{
			first_radiis[i/(N/10) - 1] = cluster[i].radius;
		}
		else
		{
			diff += fabs(cluster[i].radius - first_radiis[i/(N/10) - 1]);
		}
		
	}
	first = false;
	if(diff > difference)
		difference = diff;

	
		
}

void rminmaxcontrol(vector<Star> & cluster)  // To obtain star curves around the center
{
	double scaling = 3 * PI / 16;
	for(int i = 1; i < 11;i++)
	{
		for(int j = 0; j < cluster.size();j++)
		{
			if( (cluster[j]).id == i )
			{
				minmax << (cluster[j]).radius << endl;
				break;
			}
		}
		
	}
}

void r_u_control(vector<Star> & cluster) // result of r_u state
{
	double scaling = 3 * PI / 16;
	double av_radius, av_radial;
	av_radius = 0;
	av_radial = 0;
	int count = 0;
	for(int i =0;i<N;i++)
	{
		/*
		if((cluster[i]).radius < 1.75 and cluster[i].radius > 0.3 and abs((cluster[i]).rad_vel) < 1.0)  // there is a boundary to exclude bad asses !
		{
			r_u << (cluster[i]).radius  << endl; // x axis
			av_radius += (cluster[i]).radius;
			r_u << (cluster[i]).rad_vel << endl; // y axis
			av_radial += (cluster[i]).rad_vel;
			count++;
		}
		*/
		r_u << (cluster[i]).radius  << endl; // x axis
		av_radius += (cluster[i]).radius;
		r_u << (cluster[i]).rad_vel << endl; // y axis
		av_radial += (cluster[i]).rad_vel;
		count++;
		
	}	
	r_u << 0 << endl;  // to seperate datas

	av_radius /= count;
	av_radial /= count;
	//cout << "average radius : " << av_radius << ", average radial vel. : " << av_radial << endl;

	ru_radii.push_back(av_radius);
	ru_radial.push_back(av_radial);

}

void energychecking(vector<Star> & cluster, vector<Star *> & inds)  // check virial energy condition 
{
	double kinetic = 0;
	double potential = 0;
	int exist = 0;
	int ind;
	for(int i=0;i<N;i++)
	{
		ind = (*(inds[i])).id;
		if( 1 )
		{
			kinetic += 0.5 *  double(1)/N * (pow((cluster[ind]).tan_vel,2) + pow((cluster[ind]).rad_vel,2));  
			for(int j = i+1;j<N;j++)
			{
				potential += (double(1)/N) / (*(inds[j])).radius ;		// if radius of others larger than the star's				
			}
			potential += (i + 0.5)*(double(1)/N) / (cluster[ind]).radius ;  // if radius of others smaller than the star's	
			exist++;
		}
		
	}
	
	potential *= -0.5 * (double(1)/N) ; 
	
	cout <<  "# of survival stars:  " << exist <<   " , KE = " << kinetic << "  ,  PE = "<< potential << " ,  Tot = " << kinetic + potential << endl;  
	
	
}

void miss(vector<Star> & cluster)   // for controlling bad guys
{
	double max = 0;
	int num = 0;
	for(int i = 0;i<N;i++)
	{
		if(sqrt(pow((cluster[i]).rad_vel,2) + pow((cluster[i]).tan_vel,2)) > max)
		{
			max = sqrt(pow((cluster[i]).rad_vel,2) + pow((cluster[i]).tan_vel,2));
			num = i;
		}
	}
	error << " miss  " <<  clck << "   " << (cluster[num]).id <<  "   " <<  (cluster[num]).radius << "   " << (cluster[num]).rad_vel << "  " << (cluster[num]).tan_vel << endl;
	
}

void instability_test() // Search the radii and radial differences in time to determine whether instability occurs or not.
{
	double dif;
	double max = 0;
	for(int i=0;i<ru_radii.size()-1;i++)
	{
		dif = sqrt( pow((ru_radii[i+1] - ru_radii[i] ) , 2) + pow((ru_radial[i+1] - ru_radial[i] ) , 2) );
		/*
		if(dif > max)
			max = dif;
		*/
		max+=dif;
	}
	//cout << max << endl;
	//cout << difference << endl; // maximum radiis difference (spatial structure analogue)
}

void open_files()
{
	output.open("radiis.txt"); // open a file
	minmax.open("rminmax.txt"); // open a file
	r_u.open("r_u_plane.txt"); // open a file
	error.open("error.txt");
	script_output.open("script_output.txt" , fstream::app);
}

void close_files()
{
	output.close();
	minmax.close();
	r_u.close();
	error.close();
	script_output.close();

}