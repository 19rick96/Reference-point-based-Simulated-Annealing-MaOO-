#include <iostream>
#include <vector>
#include <time.h>
#include <ctime>
#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "rsa.h"

using namespace std;

int main()
{
	srand(time(NULL));
	rsa optimizer;
	
	ofstream myfile;
	int i;
	////////////////////////////////////////////////////	DTLZ3	/////////////////////////////////////////////////
	
	for(i=0; i<1; i++)
	{
		optimizer.define_exp(true,false,0.99,50,91,300,12,"SBX","diff",0.75,"DTLZ3",12,3,0.5,0.2,30,20,0.1,0.0,0.25);	
		myfile.open("dtlz/DTLZ3_3_d"+IntToStr(i)+".csv");
		optimizer.optimize();
		for(int i=0; i<optimizer.arch.members.size(); i++)
		{
			int l = optimizer.arch.members[i].objective.size();
			for(int j=0; j<l; j++)
			{
				myfile<<optimizer.arch.members[i].objective[j]<<" ";
			}
			myfile<<"\n";
		}
		myfile.close();
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	
	return 0;
}

