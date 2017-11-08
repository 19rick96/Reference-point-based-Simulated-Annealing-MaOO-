#ifndef _amosa2_H_INCLUDED_
#define _amosa2_H_INCLUDED_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>

#include "gen_func.h"
#include "problems.h"
#include "rand.h"

using namespace std;

class rsa
{
	public:

	double T,T_min,alpha,frac,F,CR,mut_prob,cross_prob,r_mut,prob_r1,prob_r2;
	string func,mut1,mut2,dist_m;
	int eta_c,eta_m,itr,hard_l,soft_l,gamma,divisions,dim_in,dim_out;
	vector<double> min_x,max_x,min_o,max_o;

	archive arch;
	vector< vector<double> > ref_pts,IGD_pts;
	vector<int> assoc_arr;
	vector<double> prob_array;
	vector<double> val_arr, intercepts, ideal;
	bool flag1,flag2;

	void define_exp(bool flg1,bool flg2,double al,int itrr,int hard,int soft,int div,string m1,string m2,double frc,string fun,int d_in,int d_out,double ff = 0.5,double ccr = 0.2,int e_c = 30,int e_m=20,double pro_r1 = 0.1,double pro_r2 = 0.1,double r_m=0.25)
	{
		flag2 = flg2;
		flag1 = flg1;

		T = 100.0;
		T_min = 0.000001;
		alpha = al;
		itr = itrr;

		hard_l = hard;
		soft_l = soft;
		divisions = div;
		gamma = 2;

		mut2 = m2;
		mut1 = m1;

		frac = frc;	
		F = ff;
		CR = ccr;
		eta_c = e_c;	
		eta_m = e_m;
		r_mut = r_m;
		prob_r1 = pro_r1;
		prob_r2 = pro_r2;

		func = fun;
		dist_m = "d1";
		dim_in = d_in;
		dim_out = d_out;

		mut_prob = 1.0/dim_in;
		cross_prob = 1.0;
		
		min_x.assign(dim_in,0);
		max_x.assign(dim_in,1);
		
		if(func == "ZDT4")
		{
			min_x.clear();
			max_x.clear();
			min_x.assign(dim_in,-5);
			max_x.assign(dim_in,5);
			min_x[0] = 0;
			max_x[0] = 1;
		}
		if(func == "WFG1" || func == "WFG2" || func == "WFG3" || func == "WFG4" || func == "WFG5" || func == "WFG6" || func == "WFG7" || func == "WFG8" || func == "WFG9")
		{
			for(int i=0;i<dim_in;i++)
				max_x[i] = max_x[i]*2*(i+1);
		}

		ref_pts = readfile(dim_out);

		IGD_pts = read_ref_pts(func,dim_out);

		arch.members.clear();
		prob_array.clear();

		for(int i=0; i<(soft_l*gamma); i++)
			arch.members.push_back(create_random_pt());

		min_o = arch.members[0].objective;
		max_o = arch.members[0].objective;

		for(int i=0; i<arch.members.size(); i++)
		{
			for(int j=0; j<arch.members[0].objective.size(); j++)
			{
				if(arch.members[i].objective[j] > max_o[j])
					max_o[j] = arch.members[i].objective[j];
				if(arch.members[i].objective[j] < min_o[j])
					min_o[j] = arch.members[i].objective[j];
			}
		}
	}

	individual create_random_pt()
	{	
		individual rpt;
		for(int i=0; i<dim_in; i++)
		{
			double r = random(min_x[i],max_x[i]);
			rpt.decision.push_back(r);
		}
		rpt.objective.resize(dim_out);
		evaluate(func,rpt.decision,rpt.objective);
		return rpt;
	}

	individual perturb(individual pt, string mut)
	{
		individual new_pt;
		new_pt.decision.resize(dim_in);
		new_pt.objective.resize(dim_out);
		if(mut == "poly_mut")
		{
			polynomial_mutate(new_pt.decision,pt.decision,mut_prob,eta_m,min_x,max_x);			
		}
		else if(mut == "SBX")
		{
			SBX_mutate(new_pt.decision,arch,pt.decision,eta_c,cross_prob,min_x,max_x);
		}
		else if(mut == "real_mut")
		{
			real_mutate(new_pt.decision,pt.decision,r_mut,min_x,max_x);
		}
		else
		{
			diff_mutate(new_pt.decision,arch,pt.decision,F,CR,min_x,max_x);
		}
		evaluate(func,new_pt.decision,new_pt.objective);
		return new_pt;
	}

	double evaluate_opt()
	{
		return IGD(IGD_pts,arch);
	}

	void dom_amt_archive(int &cnt,vector<int> &l,individual &x1)
	{
		cnt = 0;
		for(int i=0; i<arch.members.size(); i++)
		{
			double ch = check_dominate(arch.members[i].objective,x1.objective);
			if(ch == 1)
			{
				cnt++;
			}
			if(ch == -1)
				l.push_back(i);
		}
	}

	void remove_from_archive(vector<int> l)
	{
		if(arch.members.size() - l.size() > 3)
		{
			for(int i=0; i<l.size(); i++)
			{
				int ind = l[l.size()-i-1];
				arch.members.erase(arch.members.begin()+ind);
			}
		}
	}

	void update_max_min(individual &x1)
	{
		for(int i=0; i<x1.objective.size(); i++)
		{
			if(x1.objective[i] > max_o[i])
				max_o[i] = x1.objective[i];
			if(x1.objective[i] < min_o[i])
				min_o[i] = x1.objective[i];
		}
	}

	void cluster(int limit)
	{
		arch = nd(arch);
		vector< vector<double> > norm(arch.members.size(),vector<double> (arch.members[0].objective.size()));
		normalize(norm,arch,intercepts,ideal);

		vector< vector<double> > archive_ptr;
		vector<int> pop_arr(ref_pts.size());
		associate(pop_arr,archive_ptr,norm,ref_pts,dist_m);

		while(arch.members.size() > limit)
		{
			int ind = -1, max_pop = -1;
			for(int i=0; i<pop_arr.size(); i++)
			{
				if(pop_arr[i] > max_pop)
				{
					max_pop = pop_arr[i];
					ind = i;
				}
			}
			double max_d = -1;
			int ind_r = -1;
			int cnt = 0;
			for(int i=0; i<norm.size(); i++)
			{
				if((archive_ptr[i][0] == ind) && (archive_ptr[i][1]>max_d))
				{
					ind_r = i;
					max_d = archive_ptr[i][1];
					cnt++;
					if(cnt == max_pop)
						break;
				}
			}
			arch.members.erase(arch.members.begin()+ind_r);
			norm.erase(norm.begin() + ind_r);
			pop_arr[archive_ptr[ind_r][0]]--;
			archive_ptr.erase(archive_ptr.begin()+ind_r);
		}
	}

	void optimize()
	{
		individual current_pt = create_random_pt();
		int cnt=0;
		while(T > T_min)
		{
				for(int i=0; i<itr; i++)
				{
					int r = int(random(0,arch.members.size()));
					current_pt = arch.members[r];
					individual new_pt;
					if(i > (itr*frac))
					{	
						new_pt = perturb(current_pt,mut2);
						if(flag2 == true)
							new_pt = perturb(new_pt,"poly_mut");
						if(random(0,1) <= prob_r2)
							new_pt = perturb(new_pt,"real_mut");
					}
					else 
					{
						new_pt = perturb(current_pt,mut1);
						if(flag1 == true)
							new_pt = perturb(new_pt,"poly_mut");
						if(random(0,1) <= prob_r1)
							new_pt = perturb(new_pt,"real_mut");
					}
					int dom_stat = check_dominate(current_pt.objective,new_pt.objective);

					vector<int> l;
					int k;
					dom_amt_archive(k,l,new_pt);
			
					if(dom_stat != 1)
					{
						if((l.size()>=k))
						{
							arch.members.push_back(new_pt);
							if(l.size() == 0) 
							{
								if(arch.members.size() > soft_l)
									{cluster(hard_l);}
							}
							else
							{
								remove_from_archive(l);
							}
						}
						else
						{
							double prob = exp((double(l.size())-double(k))/(double(arch.members.size())*T));
							cnt++;
							prob_array.push_back(prob);
							if(random(0,1) <= prob)
							{
								arch.members.push_back(new_pt);
								if(l.size() == 0) 
								{
									if(arch.members.size() > soft_l)
										{cluster(hard_l);}
								}
							}
						}
					}
					else
					{
						double prob = exp((double(l.size())-double(k))/(double(arch.members.size())*T));
						cnt++;
						prob_array.push_back(prob);
						if(random(0,1) < prob)
						{
							arch.members.push_back(new_pt);
							if(l.size() == 0) 
							{
								if(arch.members.size() > soft_l)
									{cluster(hard_l);}
							}
						}
					}
					cout<<T<<" , "<<i<<" , "<<func<<" , "<<dim_out<<"\n";
			}
			T = alpha*T;
		}
		cluster(hard_l);
		cout<<cnt<<"\n";
	}
};

#endif



