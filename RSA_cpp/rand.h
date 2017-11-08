#ifndef _rand_H_INCLUDED_
#define _rand_H_INCLUDED_

#include <iostream>
#include <vector>
#include <time.h>
#include <string>
#include <ctime>

using namespace std;

double random(double lb,double ub)
{
	 return lb + ((static_cast<double>(rand())/RAND_MAX)*(ub - lb));
}

vector<int> sample(int l, int u, int n)
{
	vector<int> r(n,l-1);
	int cnt = 0;
	bool flag = false;
	while (cnt < n)
	{
		int a = int(random(double(l),double(u)));	
		flag = true;
		for(int i=0; i<r.size(); i++)
		{
			if(a==r[i])
			{
				flag = false;
				break;
			}			
		}
		if(flag == true)
		{
			r[cnt] = a;
			cnt++;
		}
	}
	return r;
}

#endif
