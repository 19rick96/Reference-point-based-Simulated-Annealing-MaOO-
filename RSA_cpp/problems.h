#ifndef _problems_H_INCLUDED_
#define _problems_H_INCLUDED_

#include<vector>
#include<cstdlib>
#include <cmath>
#include <limits>
#include<string>
#include "Toolkit/ExampleProblems.h"

using namespace WFG::Toolkit::Examples;
using namespace std;

const double PI = 3.1415926;

void DTLZ1(vector<double> &out, vector<double> &in)
{
	int n_var = in.size();
	int n_obj = out.size();
	int k = n_var - n_obj + 1;

	double g = 0;
	for (size_t i = n_obj-1; i < n_var; i += 1)
	{
		g += ((in[i]-0.5)*(in[i]-0.5)) - cos(20*PI*(in[i]-0.5));
	}
	g = (k + g)*100;

	for(int i=0; i<n_obj; i++)
		out[i] = 0;

	for (size_t m = 0; m < n_obj; m += 1)
	{
		double product = 0.5*(1+g);
		size_t i = 0;
		for (; n_obj >= 2+m && i <= n_obj-2-m; i += 1)
		{
			product *= in[i];
		}
		if (m > 0)
		{
			product *= (1 - in[i]);
		}
		out[m] = product;
	}
}

void DTLZ2(vector<double> &out, vector<double> &in)
{
	int n_var = in.size();
	int n_obj = out.size();
	int k = n_var - n_obj + 1;

	double g = 0;
	for (size_t i = n_obj-1; i < n_var; i += 1)
	{
		g += ((in[i]-0.5)*(in[i]-0.5));
	}

	for(int i=0; i<n_obj; i++)
		out[i] = 0;
	for (size_t m = 0; m < n_obj; m += 1)
	{
		double product = (1+g);
		size_t i=0;
		for (; i+m<=n_obj-2; i+=1)
		{
			product *= cos(in[i]*PI/2);
		}
		if (m > 0)
			product *= sin(in[i]*PI/2);
		
		out[m] = product;
	}
}

void DTLZ3(vector<double> &out, vector<double> &in)
{
	int n_var = in.size();
	int n_obj = out.size();
	int k = n_var - n_obj + 1;

	double g = 0;
	for (size_t i = n_obj-1; i <n_var; i += 1)
	{
		g += ((in[i]-0.5)*(in[i]-0.5)) - cos(20*PI*(in[i]-0.5));
	}
	g = (k + g)*100;

	for(int i=0; i<n_obj; i++)
		out[i] = 0;
	for (size_t m = 0; m < n_obj; m += 1)
	{
		double product = (1+g);
		size_t i=0;
		for (; i+m<=n_obj-2; i+=1)
		{
			product *= cos(in[i]*PI/2);
		}
		if (m > 0)
			product *= sin(in[i]*PI/2);
		
		out[m] = product;
	}
}

void DTLZ4(vector<double> &out, vector<double> &in, int a=100)
{
	int n_var = in.size();
	int n_obj = out.size();
	int k = n_var - n_obj + 1;

	double g = 0;
	for (size_t i = n_obj-1; i < n_var; i += 1)
	{
		g += ((in[i]-0.5)*(in[i]-0.5));
	}

	for(int i=0; i<n_obj; i++)
		out[i] = 0;
	for (size_t m = 0; m < n_obj; m += 1)
	{
		double product = (1+g);
		size_t i=0;
		for (; i+m<=n_obj-2; i+=1)
		{
			product *= cos( pow(in[i], a)*PI/2 );
		}
		if (m > 0)
			product *= sin( pow(in[i], a)*PI/2 );
		
		out[m] = product;
	}
}

void DTLZ5(vector<double> &out, vector<double> &in)
{
	int n_var = in.size();
	int n_obj = out.size();
	int k = n_var - n_obj + 1;

	double g = 0;
	for (size_t i = n_obj-1; i < n_var; i += 1)
	{
		g += ((in[i]-0.5)*(in[i]-0.5));
	}

	vector<double> theta(n_var);
	theta[0] = in[0]*PI/2.0; 
	for (size_t i=1; i<theta.size(); i+=1)
	{
		theta[i] = PI/(4*(1+g))*(1+(2*g*in[i]));
	}

	for(int i=0; i<n_obj; i++)
		out[i] = 0;
	for (size_t m = 0; m < n_obj; m += 1)
	{
		double product = (1+g);
		size_t i=0;
		for (; i+m<=n_obj-2; i+=1)
		{
			product *= cos(theta[i]);
		}
		if (m > 0)
			product *= sin(theta[i]);
		
		out[m] = product;
	}
}

void DTLZ6(vector<double> &out, vector<double> &in)
{
	int n_var = in.size();
	int n_obj = out.size();
	int k = n_var - n_obj + 1;

	double g = 0;
	for (size_t i = n_obj-1; i < n_var; i += 1)
	{
		g += pow(in[i], 0.1);
	}

	// the following is the same as DTLZ5
	vector<double> theta(n_var);
	theta[0] = in[0]*PI/2.0; 
	for (size_t i=1; i<theta.size(); i+=1)
	{
		theta[i] = PI/(4*(1+g))*(1+(2*g*in[i]));
	}

	for(int i=0; i<n_obj; i++)
		out[i] = 0;
	for (size_t m = 0; m < n_obj; m += 1)
	{
		double product = (1+g);
		size_t i=0;
		for (; i+m<=n_obj-2; i+=1)
		{
			product *= cos(theta[i]);
		}
		if (m > 0)
			product *= sin(theta[i]);
		
		out[m] = product;
	}
}

void DTLZ7(vector<double> &out, vector<double> &in)
{
	int n_var = in.size();
	int n_obj = out.size();
	int k = n_var - n_obj + 1;

	for (size_t m = 0; m < n_obj-1; m += 1)
	{
		out[m] = in[m];
	}

	double g = 0;
	for (size_t i = n_obj-1; i <n_var; i += 1)
	{
		g += in[i];
	}
	g = 1 + 9*g/k;

	double h = n_obj;
	for (size_t i = 0; i < n_obj-1; i += 1)
	{
		h -= out[i]/(1+g)*(1+sin(3*PI*out[i]));
	}

	out[n_obj-1] = (1+g)*h;
}

void ZDT1(vector<double> &out, vector<double> &in)
{
	double f1 = in[0];					//	30-var prob
	double s = 0.0;
	for(int i=1; i<in.size(); i++)
	{
		s += (in[i]);
	} 
	double g = 1.0 + ((9.0*s)/(in.size()-1.0));
	double f2 = g * (1.0 - sqrt(f1/g));
	out[0] = f1;
	out[1] = f2;
}

void ZDT2(vector<double> &out, vector<double> &in)
{								//	30-var prob
	double f1 = in[0];
	double s = 0.0;
	for(int i=1; i<in.size(); i++)
	{
		s += (in[i]);
	}
	double g = 1.0 + ((9.0*s)/(in.size()-1.0));
	double f2 = g * (1.0 - ((f1/g)*(f1/g)));
	out[0] = f1;
	out[1] = f2;
}

void ZDT3(vector<double> &out, vector<double> &in)
{								//	30-var prob
	double f1 = in[0];
	double s = 0.0;
	for(int i=1; i<in.size(); i++)
	{
		s += (in[i]);
	}
	double g = 1.0 + ((9.0*s)/(in.size()-1.0));
	double f2 = g * (1.0 - sqrt(f1/g) - ((f1/g)*sin(10.0*PI*f1)));
	out[0] = f1;
	out[1] = f2;
}

void ZDT4(vector<double> &out, vector<double> &in)
{								//	10-var prob
	double f1 = in[0];
	double s = 0.0;
	for(int i=1; i<in.size(); i++)
	{
		s += ((in[i]*in[i]) - (10.0*cos(4.0*PI*in[i])));
	} 
	double g = 1.0 + (10.0*(in.size())) + s;
	double f2 = g * (1.0 - sqrt(f1/g));
	out[0] = f1;
	out[1] = f2;
}

void ZDT6(vector<double> &out, vector<double> &in)
{								//	10-var prob
	double f1 = in[0];
	double s = 0.0;
	for(int i=1; i<in.size(); i++)
	{
		s += (in[i]);
	}
	double g = 1.0 + (9.0*pow((s/9.0),0.25));
	double f2 = g * (1.0 - ((f1/g)*(f1/g)));
	out[0] = f1;
	out[1] = f2;
}

void WFG1(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG1(in,k,n_obj);
}

void WFG2(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG2(in,k,n_obj);
}

void WFG3(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG3(in,k,n_obj);
}

void WFG4(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG4(in,k,n_obj);
}

void WFG5(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG5(in,k,n_obj);
}

void WFG6(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG6(in,k,n_obj);
}

void WFG7(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG7(in,k,n_obj);
}

void WFG8(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG8(in,k,n_obj);
}

void WFG9(vector<double> &out, vector<double> &in)
{
	int n_obj = out.size();
	int k = (n_obj-1);
	out = Problems::WFG9(in,k,n_obj);
}

void evaluate(string prob, vector<double> &in, vector<double> &out)
{
	if(prob == "DTLZ1")
		DTLZ1(out,in);
	else if(prob == "DTLZ2")
		DTLZ2(out,in);
	else if(prob == "DTLZ3")
		DTLZ3(out,in);
	else if(prob == "DTLZ4")
		DTLZ4(out,in);
	else if(prob == "DTLZ5")
		DTLZ5(out,in);
	else if(prob == "DTLZ6")
		DTLZ6(out,in);
	else if(prob == "ZDT1")
		ZDT1(out,in);
	else if(prob == "ZDT2")
		ZDT2(out,in);
	else if(prob == "ZDT3")
		ZDT3(out,in);
	else if(prob == "ZDT4")
		ZDT4(out,in);
	else if(prob == "ZDT6")
		ZDT6(out,in);
	else if(prob == "WFG1")
		WFG1(out,in);
	else if(prob == "WFG2")
		WFG2(out,in);
	else if(prob == "WFG3")
		WFG3(out,in);
	else if(prob == "WFG4")
		WFG4(out,in);
	else if(prob == "WFG5")
		WFG5(out,in);
	else if(prob == "WFG6")
		WFG6(out,in);
	else if(prob == "WFG7")
		WFG7(out,in);
	else if(prob == "WFG8")
		WFG8(out,in);
	else if(prob == "WFG9")
		WFG9(out,in);
	else
		DTLZ7(out,in);
}

#endif

