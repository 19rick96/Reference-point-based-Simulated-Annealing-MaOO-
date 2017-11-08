import numpy as np
import copy
from math import *
import random
import time
import operator
import csv
from numpy import linalg as LA
from gen_func import *
from problems import *

class RSA(object):
	def define_exp(self,flg1,flg2,al,itrr,hard,soft,div,m1,m2,frc,fun,d_in,d_out,ff = 0.5,ccr = 0.2,e_c = 30,e_m=20,pro_r1 = 0.1,pro_r2 = 0.1,r_m=0.25):
		self.flag2 = flg2
		self.flag1 = flg1

		self.T = 100.0
		self.T_min = 0.000001
		self.alpha = al
		self.itr = itrr

		self.hard_l = hard
		self.soft_l = soft
		self.divisions = div
		self.gamma = 2

		self.mut2 = m2
		self.mut1 = m1

		self.frac = frc	
		self.F = ff
		self.CR = ccr
		self.eta_c = e_c	
		self.eta_m = e_m
		self.r_mut = r_m
		self.prob_r1 = pro_r1
		self.prob_r2 = pro_r2

		self.func = fun
		self.dim_in = d_in
		self.dim_out = d_out

		self.mut_prob = 1.0/self.dim_in
		self.cross_prob = 1.0
		
		self.min_x = np.asarray([0]*self.dim_in)
		self.max_x = np.asarray([1]*self.dim_in)
		
		self.ref_pts = readfile(self.dim_out)
		self.IGD_pts = read_ref_pts(self.func,self.dim_out)
		
		self.arch = []
		for i in xrange(0,self.soft_l*self.gamma):
			self.arch.append(self.create_random_pt())
		
		#	self.arch = nd(self.arch);

		self.min_o = self.arch[0][1]
		self.max_o = self.arch[0][1]

		for i in xrange(0,len(self.arch)):
			for j in xrange(0,len(self.arch[0][1])):
				if self.arch[i][1][j] > self.max_o[j]:	
					self.max_o[j] = self.arch[i][1][j]
				if self.arch[i][1][j] < self.min_o[j]:
					self.min_o[j] = self.arch[i][1][j]	
		self.intercepts = []
		self.ideal = []	

	def create_random_pt(self):
		decision = np.zeros(self.dim_in)
		for i in xrange(0,self.dim_in):
			r = self.min_x[i] + (random.random()*(self.max_x[i]-self.min_x[i]))
			decision[i] = r
		objective = evaluate(self.func,decision,self.dim_out)
		return [decision,objective]

	def perturb(self,pt,mut):
		decision = np.zeros(self.dim_in)
		if(mut == "poly_mut"):
			decision = polynomial_mutate(pt[0],self.mut_prob,self.eta_m,self.min_x,self.max_x)
		elif(mut == "SBX"):
			decision = SBX_mutate(self.arch,pt[0],self.eta_c,self.cross_prob,self.min_x,self.max_x)
		elif(mut == "real_mut"):
			decision = real_mutate(pt[0],self.r_mut,self.min_x,self.max_x)
		else:
			decision = diff_mutate(self.arch,pt[0],self.F,self.CR,self.min_x,self.max_x)
		objective = evaluate(self.func,decision,self.dim_out)
		return [decision,objective]

	def evaluate_opt(self):
		return IGD(self.IGD_pts,self.arch)

	def dom_amt_archive(self,indv):
		cnt = 0
		l = []
		for i in xrange(0,len(self.arch)):
			ch = check_dominate(self.arch[i][1],indv[1])
			if ch == 1 :
				cnt += 1
			if ch == -1 :
				l.append(i)
		return cnt,l
	
	def remove_from_archive(self,l):
		if (len(self.arch) - len(l)) > 3 : 
			for i in xrange(0,len(l)):
				ind = l[len(l)-i-1]
				self.arch.pop(ind)

	def update_max_min(self,indv):
		for i in xrange(0,len(indv[1])):
			if indv[1][i] > self.max_o[i] :
				self.max_o[i] = indv[1][i]
			if indv[1][i] < self.min_o[i] :
				self.min_o[i] = indv[1][i]

	def cluster(self,limit):
		self.arch = nd(self.arch)
		norm,self.intercepts,self.ideal = normalize(self.arch)
		pop_arr,archive_ptr = associate(norm,self.ref_pts)
		while(len(self.arch) > limit):
			ind = np.argmax(pop_arr)
			max_pop = pop_arr[ind]

			max_d = -1
			ind_r = -1
			cnt = 0
			for i in xrange(0,len(norm)):
				if((archive_ptr[i][0] == ind) and (archive_ptr[i][1]>max_d)):
					ind_r = i
					max_d = archive_ptr[i][1]
					cnt += 1
					if cnt == max_pop :
						break
			self.arch.pop(ind_r)
			norm.pop(ind_r)
			pop_arr[archive_ptr[ind_r][0]] -= 1
			archive_ptr.pop(ind_r)

	def iterate(self,T,itr):
		r = random.randint(0,len(self.arch)-1)
		current_pt = self.arch[r]
		new_pt = []
		if itr > (self.itr*self.frac) :
			new_pt = self.perturb(current_pt,self.mut2)
			if(self.flag2 == True):
				new_pt = self.perturb(new_pt,"poly_mut")
			if(random.random() <= self.prob_r2):
				new_pt = self.perturb(new_pt,"real_mut")
		else:
			new_pt = self.perturb(current_pt,self.mut1)
			if(self.flag1 == True):
				new_pt = self.perturb(new_pt,"poly_mut")
			if(random.random() <= self.prob_r1):
				new_pt = self.perturb(new_pt,"real_mut")
		dom_stat = check_dominate(current_pt[1],new_pt[1])
		k,l = self.dom_amt_archive(new_pt)
		if dom_stat != 1 :
			if(len(l) >= k):
				self.arch.append(new_pt)
				if(len(l) == 0):
					if len(self.arch) > self.soft_l :
						self.cluster(self.hard_l)
				else:
					self.remove_from_archive(l)
			else:
				prob = exp(((k - len(l))/len(self.arch))/T)
				if random.random() <= prob :
					self.arch.append(new_pt)
					if len(l) == 0:
						if len(self.arch) > self.soft_l:
							self.cluster(self.hard_l)
		else:
			prob = exp(((k - len(l))/len(self.arch))/T)
			if random.random() <= prob :
				self.arch.append(new_pt)
				if len(l) == 0:
					if len(self.arch) > self.soft_l:
						self.cluster(self.hard_l)

		print T," , ",itr," , ",self.func," , ",self.dim_out

	def optimize(self):
		current_pt = self.create_random_pt()
		while(self.T > self.T_min):
			for i in xrange(0,self.itr):
				self.iterate(self.T,i)
			self.T = self.alpha*self.T		
		self.cluster(self.hard_l)

