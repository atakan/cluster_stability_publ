#ifndef _STAR_H
#define _STAR_H


using namespace std;

class Star
{
    public:
  
        
        double rad_vel,tan_vel;  // radial and tangential velocities
        double radius; // distance to center
        long id;
        //double nradius,nrad_vel,ntan_vel;  // new values
        // phase = arctan(tan_vel / rad_vel) 
        // velocity = sqrt(rad_vel^2 + tan_vel^2)
        
        Star();   // constructor
        void update(double r, double rv, double tv);  // update the star state
        
};

Star::Star() // Star construction function --> based on Henon's initial state generation
{
	// G=1 , M=1 , R=1
	
	radius = 0; // assign radius
	rad_vel = 0; // assign radial velocity
	tan_vel = 0; // assign tangential velocity

}


void Star::update(double r, double rv, double tv)  // update the star state
{
	radius = r;
	rad_vel = rv;
	tan_vel = tv;
}

#endif
