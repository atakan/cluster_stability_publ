#include <iostream>
#include <fstream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <math.h>  
#include <vector>

#define PI 3.14159265

using namespace std;

int binary_search(double r, vector<double> array) // binary search
{
	int result;

	int start = 0;
	int finish = array.size()- 1;
	int mid = (start + finish)/2;
	while( !(mid == start || mid == finish) )
	{
		if((array[mid]) >= r)
			finish = mid;
		else if((array[mid]) <= r)
			start = mid;
		mid = (start + finish)/2;
	}
	if(r <= array[mid])
		return mid;
	else
	{
		if(array[mid+1] >= r)
			return mid+1;
		else
			return mid+2;
	}
	return mid;
}

int main(int argc, char** argv)
{
	
	double r = 8.01;
	vector<double> x;
	x.push_back(1.0);
	x.push_back(2.0);
	x.push_back(3.0);
	x.push_back(4.0);
	x.push_back(5.0);
	x.push_back(6.0);
	x.push_back(7.0);
	x.push_back(8.0);
	x.push_back(9.0);
	x.push_back(10.0);
	cout << binary_search(r,x) << endl;
	//cout << sqrt(2.0)<< endl;
	/*for(int i = 0; i< 20 ;i++)
	{
		r = ((double) rand() / (RAND_MAX));
		cout << cos ( 59.0 * PI / 180.0 ) << endl;
		cout << r << endl;
	}*/
	
	
	
	/*for (int i = 0; i< 9;i++)
	{
		
		cout << 3*(PI/16)*pow((pow(a[i],double(-2)/3) - 1),double(-1)/2)<< endl;
	}*/
	
	return 0;
}
