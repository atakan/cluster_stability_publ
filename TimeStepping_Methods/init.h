#include <iostream>
#include <fstream>
#include <stdio.h>      
#include <stdlib.h>     /* atoi calling */
#include <math.h>    
#include <vector>  /* for vectors */
#include <ctime> // for timer
#include <time.h>

#include "Star.h"  // including header of stars //

#define PI 3.1415926535897

using namespace std;

extern int N;  // number of stars in cluster
extern double clck; 

extern vector<double> ru_radii; // for check.cpp
extern vector<double> ru_radial; // for check.cpp

extern double  first_radiis [9]; // for radiis instability
extern bool first; // check first radiis
extern double difference; // to get maximum difference

//void initialize(float n , float m , vector<Star> & cluster);
