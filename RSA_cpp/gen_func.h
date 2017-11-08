#ifndef _gen_func_H_INCLUDED_
#define _gen_func_H_INCLUDED_

#include<vector>
#include<cstdlib>
#include <cmath>
#include <limits>
#include <string>
#include <fstream>
#include <sstream>

#include "rand.h"

using namespace std;

struct individual
{
	vector<double> decision,objective;
};

struct archive
{
	vector<individual> members;
};

string IntToStr(int i)
{
	ostringstream oss;
	oss << i;
	return oss.str();
}

const double EPS = 1.0e-14;

double dot(vector<double> &v1,vector<double> &v2)
{
	double d=0.0;
	for(int i = 0; i<v1.size(); i++)
	{
		d += (v1[i]*v2[i]);
	}
	return d;
}

void add_v(vector<double> &result, vector<double> &v1, vector<double> &v2)
{
	for(int i=0; i<v1.size(); i++)
	{
		result[i] = v1[i] + v2[i];
	}
}

void scalar_mul(vector<double> &n,vector<double> &v1,double c)
{
	for(int i=0; i<v1.size(); i++)
	{
		n[i] = c * v1[i];
	}
}

void subtract(vector<double> &sub,vector<double> &v1, vector<double> &v2)
{
	for(int i=0; i<v1.size(); i++)
	{
		sub[i] = v1[i] - v2[i];
	}
}

double l2_norm(vector<double> &v1)
{
	double norm=0.0;
	for(int i=0; i<v1.size(); i++)
	{
		norm += (v1[i]*v1[i]);
	}
	return sqrt(norm);
}

double l2_distance(vector<double> &v1, vector<double> &v2)
{
	double d = 0.0;
	for (int i = 0; i<v1.size() ; i++)
	{
		d += ((v1[i]-v2[i])*(v1[i]-v2[i]));
	}
	d = sqrt(d);
	return d;
}

int check_dominate(vector<double> &v1,vector<double> &v2)
{
	bool flag1=false,flag2=false;
	int result;
	for(int i = 0; i<v1.size() ; i++)
	{
		if (v1[i] < v2[i]) 
			flag1 = true;
		if (v1[i] > v2[i])
			flag2 = true;
	}
	if (flag1 == true && flag2 == false)
		result = 1;
	else if (flag1 == false && flag2 == true)
		result = -1;
	else
		result = 0;
	return result;
}

archive nd(archive data)
{
	vector<int> n;
	archive nd_data;
	for(int i = 0; i<data.members.size(); i++)
	{
		n.push_back(0);
		for(int j=0; j<data.members.size();j++)
		{
			if(j != i)
			{
				int d = check_dominate(data.members[i].objective,data.members[j].objective);
				if (d == -1)
					n[i] ++;
			}	
		}
		if(n[i] == 0)
			nd_data.members.push_back(data.members[i]);
	}
	return nd_data;
}

double d1(vector<double> &data,vector<double> &point)
{
	int size = data.size();
	double n1 = l2_norm(point);
	vector<double> sm(size),sub(size);
	scalar_mul(sm,point,dot(data,point)/(n1*n1));
	subtract(sub,data,sm);
	return l2_norm(sub);
}

double d2(vector<double> &data,vector<double> &point)
{
	int dim = data.size();
	vector<double> n(dim,1.0),sm(dim),sub1(dim),sub2(dim);
	scalar_mul(n,n,sqrt(dim));
	subtract(sub1,data,point);
	scalar_mul(sm,n,dot(sub1,n));
	subtract(sub2,sub1,sm);
	return l2_norm(sub2);
}

double IGD(vector< vector<double> > ref_points,archive &arch)
{
	double d = 0.0;
	for(int i=0; i<ref_points.size(); i++)
	{
		double min_d = l2_distance(ref_points[i],arch.members[0].objective);
		for(int j=1; j<arch.members.size(); j++)
		{
			double dd = l2_distance(ref_points[i],arch.members[j].objective);
			if(dd < min_d)
				min_d = dd;
		}
		d += min_d;
	}
	d = d/ref_points.size();
	return d;
}

void real_mutate(vector<double> &v_new,vector<double> &v,double b,vector<double> &min_x,vector<double> &max_x)
{
	for(int i = 0; i<v.size(); i++)
		v_new[i] = v[i];
	
	int r = int(random(0,v.size()));
	double y = v_new[r];
	double d_rnd = random(0,1);
	d_rnd = d_rnd - 0.5;
	double d_rnd_lap = 0.0;
	
	if( d_rnd < 0.0)
		d_rnd_lap = b*log(1.0 - (2.0*fabs(d_rnd)));
	else
		d_rnd_lap = -1.0*b*log(1.0 - (2.0*fabs(d_rnd)));
	y += d_rnd_lap;
	int i_count = 0;
	while(((y < min_x[r]) || (y > max_x[r])) && (i_count < 20))
	{
		y = v_new[r];
		d_rnd = random(0,1);
		d_rnd = d_rnd - 0.5;
		d_rnd_lap = 0.0;
	
		if( d_rnd < 0.0)
			d_rnd_lap = b*log(1.0 - (2.0*fabs(d_rnd)));
		else
			d_rnd_lap = -1.0*b*log(1.0 - (2.0*fabs(d_rnd)));
		y += d_rnd_lap;
		i_count++;
	}
	v_new[r] = y;
	if ( i_count == 30 )
	{
		if(v_new[r] < min_x[r])
			v_new[r] = min_x[r];
		else if(v_new[r] > max_x[r])
			v_new[r] = max_x[r];
	}
}

void polynomial_mutate(vector<double> &v_new,vector<double> &v,double mut_prob,double eta_m,vector<double> &min_x,vector<double> &max_x)
{
	for(int i=0; i<v.size(); i++)
	{
		v_new[i] = v[i];
		double r = random(0,1);
		if (r<mut_prob)
		{
			double y = v_new[i];
			double yl = min_x[i];
			double yu = max_x[i];
			
			double delta1 = (y-yl)/(yu-yl);
			double delta2 = (yu-y)/(yu-yl);

			double mut_pow = 1.0/(eta_m + 1.0);
			double deltaq = 0.0;

			double rnd = random(0,1);

			if(rnd <= 0.5)
			{
				double xy = 1.0 - delta1;
				double val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
                		deltaq =  pow(val,mut_pow) - 1.0;
			}
			else
			{
				double xy = 1.0-delta2;
                		double val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
                		deltaq = 1.0 - pow(val,mut_pow);
			}
			
			y = y + deltaq*(yu-yl);
			y = min(yu, max(yl, y));
			v_new[i] = y;
		}
	}
}

void SBX(vector<double> &child1,vector<double> &child2,vector<double> &parent1,vector<double> &parent2,double eta_c,double cross_prob,vector<double> &min_x,vector<double> &max_x)
{
	double r = random(0,1);
	if(r < cross_prob)
	{
		for(int i=0; i<parent1.size(); i++)
		{
			r = random(0,1);
			if(r <= 0.5 )
			{
				if(fabs(parent1[i]-parent2[i]) > EPS)
				{
					double y1 = min(parent1[i],parent2[i]),
					      y2 = max(parent1[i],parent2[i]);
					double yl = min_x[i], yu = max_x[i];
					r = random(0,1);
					double beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
					double alpha = 2.0 - pow(beta,(-1.0*(1.0 + eta_c)));
					double betaq = 0.0;
					if (r <= (1.0/alpha))
						betaq = pow((r*alpha),(1.0/(eta_c + 1.0)));
					else
						betaq = pow((1.0/(2.0 - (r*alpha))),(1.0/(eta_c + 1.0)));
					child1[i] = 0.5*(y1+y2-(betaq*(y2-y1)));

					beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
					alpha = 2.0 - pow(beta,(-1.0*(1.0 + eta_c)));
					if (r <= (1.0/alpha))
						betaq = pow((r*alpha),(1.0/(eta_c + 1.0)));
					else
						betaq = pow((1.0/(2.0 - (r*alpha))),(1.0/(eta_c + 1.0)));
					child2[i] = 0.5*((y1+y2)+(betaq*(y2-y1)));

					child1[i] = min(yu, max(yl, child1[i]));
					child2[i] = min(yu, max(yl, child2[i]));

					r = random(0,1);
					if (r <= 0.5)
					{
						swap(child1[i],child2[i]);
					}
				}
				else
				{
					child1[i] = parent1[i];
					child2[i] = parent2[i];
				}
			}
			else
			{
				child1[i] = parent1[i];
				child2[i] = parent2[i];
			}
		}
	}
	else
	{
		for(int i=0; i<parent1.size();i++)
		{
			child1[i] = parent1[i];
			child2[i] = parent2[i];
		}
	}
}

void SBX_mutate(vector<double> &v_new,archive &arch,vector<double> &p1,double eta_c,double cross_prob,vector<double> &min_x,vector<double> &max_x)
{
	int size = p1.size();
	int size2 = arch.members.size();
	int r = int(random(0,double(size2)));					//	Check !!
	vector<double> p2 = arch.members[r].decision;
	vector<double> c1(size),c2(size);
	SBX(c1,c2,p1,p2,eta_c,cross_prob,min_x,max_x);
	double rnd = random(0,1);
	if (rnd > 0.5)
		v_new = c1;
	else
		v_new = c2;
}

void diff_mutate(vector<double> &v_new,archive &arch,vector<double> &v,double F,double CR,vector<double> &min_x,vector<double> &max_x)
{
	int size=v.size();
	vector<int> a = sample(0,arch.members.size(),3);
	vector<double> p1 = arch.members[a[0]].decision;
	vector<double> p2 = arch.members[a[1]].decision;
	vector<double> p3 = arch.members[a[2]].decision;
	int k_rand = int(random(0,double(size)));
	vector<double> u(v);
	for(int k=0; k<size; k++)
	{
		double r = random(0,1);
		if((r<CR)||(k==k_rand))
		{
			u[k] = min(max_x[k],max(p1[k] + F*(p2[k] - p3[k]),min_x[k]));	
		}
	}
	v_new = u;
}

void decimal_mut(vector<double> &v_new,vector<double> &v,double mut_dec,vector<double> &min_x,vector<double> &max_x)
{
	int size = v.size();
	int res = 7;
	vector<int> p;
	for(int i=0; i<size; i++)
	{
		double w = (v[i] - min_x[i])/(max_x[i] - min_x[i]);
		for(int j=0; j<res; j++)
		{
			p.push_back(int(w*(pow(10,j+1)))%10);
		}		
	}
	
	int md = int(mut_dec*p.size());
	for(int i=0; i<md; i++)
	{
		int pos = int(random(0.0,double(p.size())));
		int val = int(random(0.0,10.0));
		p[pos] = val;
	}
	for(int j=0; j<size; j++)
	{
		int i = j*res;
		double w = 0.0;
		for(int k=0; k<res; k++)
		{
			w = w + (p[i+k]*pow(10,(res-k-1)));
		}
		w = w / pow(10,res);
		w = min_x[j] + w*(max_x[j] - min_x[j]);
		v_new[j] = w;
	}
}

double ASF(const vector<double> &objs, const vector<double> &weight)
{
	double max_ratio = -numeric_limits<double>::max();
	for (size_t f=0; f<objs.size(); f+=1)
	{
		double w = weight[f]?weight[f]:0.00001;
		max_ratio = std::max(max_ratio, objs[f]/w);
	}
	return max_ratio;
}

void GuassianElimination(vector<double> *px, vector< vector<double> > A, const vector<double> &b)
{
	vector<double> &x = *px;

    const size_t N = A.size();
    for (size_t i=0; i<N; i+=1)
    {
        A[i].push_back(b[i]);
    }

    for (size_t base=0; base<N-1; base+=1)
    {
        for (size_t target=base+1; target<N; target+=1)
        {
            double ratio = A[target][base]/A[base][base];
            for (size_t term=0; term<A[base].size(); term+=1)
            {
                A[target][term] -= A[base][term]*ratio;
            }
        }
    }

    x.resize(N);
    for (int i=N-1; i>=0; i-=1)
    {
        for (size_t known=i+1; known<N; known+=1)
        {
            A[i][N] -= A[i][known]*x[known];
        }
        x[i] = A[i][N]/A[i][i];
    }
}

vector<double> FindMaxObjectives(archive &arch)
{
	int dim = arch.members[0].objective.size();
	int a_size = arch.members.size();

	vector<double> max_point(arch.members[0].objective);

	for(int i=0; i<dim ; i++)
	{
		for(int j=1; j<a_size; j++)
		{
			double val = arch.members[j].objective[i];
			if(val > max_point[i])
				max_point[i] = val;
		}
	}
	return max_point;
}

void normalize(vector< vector<double> > &norm,archive &arch,vector<double> &intercepts,vector<double> &ideal)
{
	vector<double> id(arch.members[0].objective);
	ideal = id;
	int dim = arch.members[0].objective.size();
	int a_size = arch.members.size();
	
	//	Find Ideal Point	/////////////////////////////////////////////////

	for(int i=0; i<dim ; i++)
	{
		for(int j=1; j<a_size; j++)
		{
			double val = arch.members[j].objective[i];
			if(val < ideal[i])
				ideal[i] = val;
		}
		for(int j=0; j<a_size; j++)
		{
			norm[j][i] = arch.members[j].objective[i] - ideal[i];
		}
	}
	
	vector<int> extreme_points(dim);
	
	//	Minimize ASF to find the extreme points	//////////////////////////////////

	for (int f=0; f<dim; f++)
	{
		vector<double> w(dim, 0.000001);
		w[f] = 1.0;

		double min_ASF = numeric_limits<double>::max();
		int min_indv = dim;

		for (int i=0; i<a_size; i++)  
		{
			double asf = ASF(norm[i], w); 

			if ( asf < min_ASF )
			{
				min_ASF = asf;
				min_indv = i;
			}
		}
		extreme_points[f] = min_indv;
	}
	
	//	Construct Hyperplane	///////////////////////////////////////////////////
	
	bool duplicate = false;
	for (size_t i=0; !duplicate && i<extreme_points.size(); i+=1)
	{
		for (size_t j=i+1; !duplicate && j<extreme_points.size(); j+=1)
		{
			duplicate = (extreme_points[i] == extreme_points[j]);
		}
	}

	intercepts.assign(dim,0.0);

	bool negative_intercept = false;
	if (!duplicate)
	{
		// Find the equation of the hyperplane
		vector<double> b(dim, 1.0);
		vector< vector<double> > A;
		for (size_t p=0; p<dim; p+=1)
		{ 
			A.push_back(norm[extreme_points[p]]);
		}
		vector<double> x;
		GuassianElimination(&x, A, b);

		// Find intercepts
		for (size_t f=0; f<dim; f+=1)
		{
			intercepts[f] = 1.0/x[f];

			if(x[f] < 0)
			{
				negative_intercept = true;
				break;
			}
		}
	}
	
	if (duplicate || negative_intercept) 
	{
		vector<double> max_objs = FindMaxObjectives(arch);
		for (size_t f=0; f<dim; f+=1)
		{
			intercepts[f] = max_objs[f];
		}
	}	

	//	Normalize Objectives using the intercepts	////////////////////////////

	for (int i=0; i<a_size; i++)
	{
		for (int f=0; f<dim; f++)
		{
			if ( fabs(intercepts[f])>10e-10 ) // avoid the divide-by-zero error
				norm[i][f] = norm[i][f]/(intercepts[f]);
			else
				norm[i][f] = norm[i][f]/10e-10;
		}
	}
}

vector<double> associate_pt(vector<double> &pt, vector< vector<double> > &ref_pts, string dist)
{
	vector<double> data_pair(2);
	data_pair[0] = 0;
	if(dist == "d1")
		data_pair[1] = d1(pt,ref_pts[0]);
	else
		data_pair[1] = d2(pt,ref_pts[0]);
	for(int i=1; i<ref_pts.size(); i++)
	{
		double d_val;
		if(dist == "d1")
			d_val = d1(pt,ref_pts[i]);
		else
			d_val = d2(pt,ref_pts[i]);
			
		if(d_val<data_pair[1])
		{
			data_pair[1] = d_val;
			data_pair[0] = i;
		}
	}
	return data_pair;
}

void associate(vector<int> &pop_arr, vector< vector<double> > &archive_pointer, vector< vector<double> > &norm, vector< vector<double> > &ref_pts, string dist)
{
	for(int i=0;i<pop_arr.size();i++)
		pop_arr[i] = 0;
	for(int i=0; i<norm.size(); i++)
	{
		vector<double> data_pair = associate_pt(norm[i],ref_pts,dist);
		archive_pointer.push_back(data_pair);
		pop_arr[int(data_pair[0])]++ ;
	}
}
vector< vector<double> > readfile(int dim_out)
{
	vector< vector<double> > ref_pts;

	if(dim_out == 3)
	{
		string line;
		ifstream infile("DTLZ_pf/DTLZ1(3).csv");

		while (getline(infile, line))
		{
		    istringstream iss(line);
		    double a, b, c;
		    if (!(iss >> a >> b >> c)) { break; }
			vector<double> v(3);
			v[0] = 2*a;
			v[1] = 2*b;
			v[2] = 2*c;
			ref_pts.push_back(v);
		}
	}
	else if(dim_out == 5)
	{
		string line;
		ifstream infile("DTLZ_pf/DTLZ1(5).csv");

		while (getline(infile, line))
		{
		    istringstream iss(line);
		    vector<double> v(5);
		    if (!(iss >> v[0] >> v[1] >> v[2]>>v[3]>>v[4])) { break; }
		    for(int i=0; i<v.size(); i++)
			v[i] = 2*v[i];
		    ref_pts.push_back(v);
		}
	}
	else if(dim_out == 8)
	{
		string line;
		ifstream infile("DTLZ_pf/DTLZ1(8).csv");

		while (getline(infile, line))
		{
		    istringstream iss(line);
		    vector<double> v(8);
		    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7])) { break; }
			for(int i=0; i<v.size(); i++)
				v[i] = 2*v[i];
			ref_pts.push_back(v);
		}
	}
	else if(dim_out == 10)
	{
		string line;
		ifstream infile("DTLZ_pf/DTLZ1(10).csv");

		while (getline(infile, line))
		{
		    istringstream iss(line);
		    vector<double> v(10);
		    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7] >> v[8] >> v[9])) { break; }
			for(int i=0; i<v.size(); i++)
				v[i] = 2*v[i];
			ref_pts.push_back(v);
		}
	}
	else
	{
		string line;
		ifstream infile("DTLZ_pf/DTLZ1(15).csv");

		while (getline(infile, line))
		{
		    istringstream iss(line);
		    vector<double> v(15);
		    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7] >> v[8] >> v[9] >> v[10] >> v[11] >> v[12] >> v[13] >> v[14])) { break; }
			for(int i=0; i<v.size(); i++)
				v[i] = 2*v[i];
			ref_pts.push_back(v);	
		}
	}	
	return ref_pts;
}

vector< vector<double> > read_ref_pts(string func,int dim_out)
{
	vector< vector<double> > ref_pts;
	if(func == "DTLZ1")
	{
		if(dim_out == 3)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ1(3).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    double a, b, c;
			    if (!(iss >> a >> b >> c)) { break; }
				vector<double> v(3);
				v[0] = a;
				v[1] = b;
				v[2] = c;
				ref_pts.push_back(v);
			}
		}
		else if(dim_out == 5)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ1(5).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(5);
			    if (!(iss >> v[0] >> v[1] >> v[2]>>v[3]>>v[4])) { break; }
			    ref_pts.push_back(v);
			}
		}
		else if(dim_out == 8)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ1(8).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(8);
			    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7])) { break; }
				
				ref_pts.push_back(v);
			}
		}
		else if(dim_out == 10)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ1(10).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(10);
			    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7] >> v[8] >> v[9])) { break; }
				
				ref_pts.push_back(v);
			}
		}
		else
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ1(15).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(15);
			    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7] >> v[8] >> v[9] >> v[10] >> v[11] >> v[12] >> v[13] >> v[14])) { break; }
				
				ref_pts.push_back(v);	
			}
		}	
	}
	else
	{
		if(dim_out == 3)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ2(3).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    double a, b, c;
			    if (!(iss >> a >> b >> c)) { break; }
				vector<double> v(3);
				v[0] = a;
				v[1] = b;
				v[2] = c;
				ref_pts.push_back(v);
			}
		}
		else if(dim_out == 5)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ2(5).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(5);
			    if (!(iss >> v[0] >> v[1] >> v[2]>>v[3]>>v[4])) { break; }
			   
			    ref_pts.push_back(v);
			}
		}
		else if(dim_out == 8)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ2(8).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(8);
			    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7])) { break; }
				
				ref_pts.push_back(v);
			}
		}
		else if(dim_out == 10)
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ2(10).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(10);
			    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7] >> v[8] >> v[9])) { break; }
				
				ref_pts.push_back(v);
			}
		}
		else
		{
			string line;
			ifstream infile("DTLZ_pf/DTLZ2(15).csv");

			while (getline(infile, line))
			{
			    istringstream iss(line);
			    vector<double> v(15);
			    if (!(iss >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5] >> v[6] >> v[7] >> v[8] >> v[9] >> v[10] >> v[11] >> v[12] >> v[13] >> v[14])) { break; }
				
				ref_pts.push_back(v);	
			}
		}
	}
	return ref_pts;
}

class form_ref_points
{
	public:
		
	int M,div;
	vector< vector<double> > points;

	form_ref_points()
	{
		M = 2;
		div = 12;
	}

	form_ref_points(int m,int divisions)
	{
		M = m-1;
		div = divisions;
	}
	
	void recursive(vector<double> arr, int d, int l)
	{
		vector<double> arr_c(arr);
		if(d == M-1)
			points.push_back(arr_c);
		else
		{
			for(int i=0; i<l; i++)
			{
				double node_val = double(i)/double(div);
				vector<double> arr_next(arr_c);
				arr_next.push_back(node_val);
				recursive(arr_next,d+1,l-i);
			}
		}
	}

	void form()
	{
		vector<double> layer;
		for(int i=0; i<(div+1); i++)
			layer.push_back(double(i)/double(div));
		for(int i=0; i<layer.size(); i++)
		{
			vector<double> l1;
			l1.push_back(layer[i]);
			recursive(l1,0,layer.size()-i);
		}
		for(int i=0; i<points.size(); i++)
		{
			double s = 0.0;
			for(int j=0; j<points[i].size(); j++)
				s += points[i][j];
			points[i].push_back(1.0 - s);
		}
	}
	
};
	
#endif

