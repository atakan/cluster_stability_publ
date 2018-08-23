#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <math.h> 

#include "Star.h"  // header file of class

#define PI 3.14159265

using namespace std;


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




