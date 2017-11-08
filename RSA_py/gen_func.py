import numpy as np
import copy
from math import *
import random
import time
import operator
import csv
from numpy import linalg as LA

EPS = 1.0e-14

def check_dominate(v1,v2):
	s = (v1<v2).sum()
	if s == len(v1):
		return 1
	elif s == 0:
		return -1
	else :
		return 0

def nd(archive):
	nd_data = []
	n = []
	for i in xrange(0,len(archive)):
		n.append(0)
		for j in xrange(0,len(archive)):
			if j != i :
				if(check_dominate(archive[i][1],archive[j][1]) == -1):
					n[i] = n[i] + 1
		if n[i] == 0 :
			nd_data.append(archive[i])
	return nd_data

def d1(data,point):
	size = len(data)
	n1 = LA.norm(point)
	point = (np.dot(data,point)/(n1**2))*point
	return LA.norm(data - point)	

def IGD(ref_points,archive):
	d = 0.0
	for i in xrange(0,len(ref_points)):
		min_d = LA.norm(ref_points[i]-archive[0][1])
		for j in xrange(1,len(archive)):
			dd = LA.norm(ref_points[i]-archive[j][1]) 
			if dd < min_d :
				min_d = dd
		d = d + min_d
	d = d/len(ref_points)
	return d

def real_mutate(v,b,min_x,max_x):
	v_new = copy.deepcopy(v)
	r = random.randint(0,len(v)-1)
	y = v_new[r]
	d_rnd = random.random() - 0.5
	d_rnd_lap = 0.0
	
	if( d_rnd < 0 ):
		d_rnd_lap = b*log(1.0 - 2.0*abs(d_rnd))
	else:
		d_rnd_lap = -1.0*b*log(1.0 - 2.0*abs(d_rnd))
	y = y + d_rnd_lap
	i_count = 0
	while(((y < min_x[r]) or (y > max_x[r])) and (i_count < 20)):
		y = v_new[r]
		d_rnd = random.random(0,1) - 0.5
		d_rnd_lap = 0.0
	
		if( d_rnd < 0.0):
			d_rnd_lap = b*log(1.0 - (2.0*abs(d_rnd)))
		else:
			d_rnd_lap = -1.0*b*log(1.0 - (2.0*abs(d_rnd)))
		y = y + d_rnd_lap
		i_count = i_count + 1
	v_new[r] = y
	if ( i_count == 20 ):
		if(v_new[r] < min_x[r]):
			v_new[r] = min_x[r]
		elif(v_new[r] > max_x[r]):
			v_new[r] = max_x[r]
	return v_new

def polynomial_mutate(v,mut_prob,eta_m,min_x,max_x):
	v_new = copy.deepcopy(v)
	for i in xrange(0,len(v)):
		r = random.random()
		if (r < mut_prob):
			y = v_new[i]
			yl = min_x[i]
			yu = max_x[i]
			
			delta1 = (y-yl)/(yu-yl)
			delta2 = (yu-y)/(yu-yl)

			mut_pow = 1.0/(eta_m + 1.0)
			deltaq = 0.0

			rnd = random.random()

			if(rnd <= 0.5):
				xy = 1.0 - delta1
				val = (2.0*rnd)+((1.0-(2.0*rnd))*(xy**(eta_m+1.0)))
                		deltaq = (val**mut_pow) - 1.0
			else:
				xy = 1.0 - delta2
                		val = (2.0*(1.0-rnd))+(2.0*(rnd-0.5)*(xy**(eta_m+1.0)))
                		deltaq = 1.0 - (val**mut_pow)
			
			y = y + deltaq*(yu-yl)
			y = min(yu, max(yl, y))
			v_new[i] = y
	return v_new

def SBX(parent1,parent2,eta_c,cross_prob,min_x,max_x):
	child1 = copy.deepcopy(parent1)
	child2 = copy.deepcopy(parent2)
	r = random.random()
	if(r < cross_prob):
		for i in xrange(0,len(parent1)):
			r = random.random()
			if(r <= 0.5 ):
				if(abs(parent1[i]-parent2[i]) > EPS):
					y1 = min(parent1[i],parent2[i])
					y2 = max(parent1[i],parent2[i])
					yl = min_x[i] 
					yu = max_x[i]
					r = random.random()
					beta = 1.0 + (2.0*(y1-yl)/(y2-y1))
					alpha = 2.0 - (beta**(-1.0*(1.0 + eta_c)))
					betaq = 0.0
					if (r <= (1.0/alpha)):
						betaq = (r*alpha)**(1.0/(eta_c + 1.0))
					else:
						betaq = (1.0/(2.0 - (r*alpha)))**(1.0/(eta_c + 1.0))
					child1[i] = 0.5*(y1+y2-(betaq*(y2-y1)))

					beta = 1.0 + (2.0*(yu-y2)/(y2-y1))
					alpha = 2.0 - (beta**(-1.0*(1.0 + eta_c)))
					if (r <= (1.0/alpha)) :
						betaq = (r*alpha)**(1.0/(eta_c + 1.0))
					else:
						betaq = (1.0/(2.0 - (r*alpha)))**(1.0/(eta_c + 1.0))
					child2[i] = 0.5*((y1+y2)+(betaq*(y2-y1)))

					child1[i] = min(yu, max(yl, child1[i]))
					child2[i] = min(yu, max(yl, child2[i]))

					r = random.random()
					if (r <= 0.5):
						sw = child1[i]
						child1[i] = child2[i]
						child2[i] = sw
	return child1,child2

def SBX_mutate(arch,p1,eta_c,cross_prob,min_x,max_x):
	size = len(p1)
	size2 = len(arch)
	r = random.randint(0,size2-1)				
	p2 = copy.deepcopy(arch[r][0])
	c1,c2 = SBX(p1,p2,eta_c,cross_prob,min_x,max_x)
	v_new = []
	rnd = random.random()
	if (rnd > 0.5):
		v_new = c1
	else:
		v_new = c2
	return v_new

def diff_mutate(archive,v,F,CR,min_x,max_x):
	size = len(v)
	a = random.sample(range(0,len(archive)),3)
	p1 = copy.deepcopy(archive[a[0]][0])
	p2 = copy.deepcopy(archive[a[1]][0])
	p3 = copy.deepcopy(archive[a[2]][0])
	
	k_rand = random.randint(0,size-1)
	u = copy.deepcopy(v)
	for k in xrange(0,size):
		r = random.random()
		if((r<CR) or (k==k_rand)):
			u[k] = min(max_x[k],max(p1[k] + F*(p2[k] - p3[k]),min_x[k]))
	return u

def ASF(objs,weight):
	w = weight[0]
	if w == 0 :
		w = 0.00001
	max_ratio = objs[0]/w
	for i in xrange(1,len(objs)):
		w = weight[i]
		if w == 0:
			w = 0.00001
		max_ratio = max(max_ratio,objs[i]/w)
	return max_ratio

def FindMaxObjectives(archive):
	dim = len(archive[0][1])
	a_size = len(archive)

	max_point = copy.deepcopy(archive[0][1])

	for i in xrange(0,dim):
		for j in xrange(1,a_size):
			val = archive[j][1][i]
			if(val > max_point[i]):
				max_point[i] = val
	return max_point

def normalize(archive):
	ideal = copy.deepcopy(archive[0][1])
	dim = len(archive[0][1])
	a_size = len(archive)

	norm = []
	for i in xrange(0,dim):
		for j in xrange(1,a_size):
			val = archive[j][1][i]
			if val < ideal[i] :
				ideal[i] = val
	for i in xrange(0,a_size):
		norm.append(archive[i][1] - ideal)
	
	extreme_points = []
	for f in xrange(0,dim):
		w = np.asarray([0.000001]*dim)
		w[f] = 1.0
		min_ASF = ASF(norm[0],w)
		min_indv = 0
		for i in xrange(1,a_size):
			asf = ASF(norm[i],w)
			if(asf < min_ASF):
				min_ASF = asf
				min_indv = i
		extreme_points.append(min_indv)

	duplicate = False
	i = 0
	while((i<len(extreme_points)) and (duplicate == False)):
		j = i + 1
		while((j<len(extreme_points)) and (duplicate == False)):
			duplicate = (extreme_points[i] == extreme_points[j])
			j += 1
		i += 1

	intercepts = np.asarray([0.0]*dim)
	negative_intercept = False
	if( duplicate == False):
		b = np.asarray([1.0]*dim)
		A = []
		for p in xrange(0,dim):
			A.append(copy.deepcopy(norm[extreme_points[p]]))
		A = np.asarray(A)
		x = np.dot(LA.inv(A),b)
		for f in xrange(0,dim):
			intercepts[f] = 1.0/x[f]
			if(x[f] < 0):
				negative_intercept = True
				break

	if( duplicate or negative_intercept):
		max_objs = FindMaxObjectives(archive)
		for f in xrange(0,dim):
			intercepts[f] = max_objs[f]

	for i in xrange(0,a_size):
		for f in xrange(0,dim):
			if(abs(intercepts[f]) > 10e-10):
				norm[i][f] = norm[i][f]/intercepts[f]
			else:
				norm[i][f] = norm[i][f]/10e-10

	return norm,intercepts,ideal

def associate_pt(pt,ref_pts):
	data_pair = []
	data_pair.append(0)
	data_pair.append(d1(pt,ref_pts[0]))
	for i in xrange(1,len(ref_pts)):
		d_val = d1(pt,ref_pts[i])
		if d_val < data_pair[1] : 
			data_pair[1] = d_val
			data_pair[0] = i
	return data_pair

def associate(norm,ref_pts):
	pop_arr = np.zeros(len(ref_pts))
	archive_pointer = []
	for i in xrange(0,len(norm)):
		data_pair = associate_pt(norm[i],ref_pts)
		archive_pointer.append(data_pair)
		pop_arr[int(data_pair[0])] += 1
	return pop_arr,archive_pointer

def readfile(dim_out):
	fname = "DTLZ/DTLZ1("+str(dim_out)+").csv"
	ref_points = []
	with open(fname) as csvfile:
		csvreader = csv.reader(csvfile,delimiter=" ")
		for line in csvreader:
		        ref_points.append(line)

	for i in range(0,len(ref_points)):
		x = []
		for j in range(0,len(ref_points[i])-1):
		        x.append(float(ref_points[i][j]))
		ref_points[i] = np.asarray(x)
	return ref_points

def read_ref_pts(func,dim_out):
	ref_points = []
	if func == "DTLZ1":
		fname = "DTLZ/DTLZ1("+str(dim_out)+").csv"
		with open(fname) as csvfile:
			csvreader = csv.reader(csvfile,delimiter=" ")
			for line in csvreader:
				ref_points.append(line)

		for i in range(0,len(ref_points)):
			x = []
			for j in range(0,len(ref_points[i])-1):
				x.append(float(ref_points[i][j]))
			ref_points[i] = np.asarray(x)	
	else:
		fname = "DTLZ/DTLZ2("+str(dim_out)+").csv"
		with open(fname) as csvfile:
			csvreader = csv.reader(csvfile,delimiter=" ")
			for line in csvreader:
				ref_points.append(line)

		for i in range(0,len(ref_points)):
			x = []
			for j in range(0,len(ref_points[i])-1):
				x.append(float(ref_points[i][j]))
			ref_points[i] = np.asarray(x)
	return ref_points

class form_ref_pts(object):
	def __init__(self,m,divisions):
		self.M = m-1
		self.div = divisions
		self.points = []

	def recursive(self,arr,d,l):
		arr_c = copy.deepcopy(arr)
		if d == self.M-1:
			self.points.append(arr_c)
		else:
			for i in xrange(0,l):
				node_val = float(i)/float(self.div)
				arr_next = copy.deepcopy(arr_c)
				arr_next.append(node_val)
				self.recursive(arr_next,d+1,l-i)

	def form(self):
		layer = []
		for i in xrange(0,div+1):
			layer.append(float(i)/float(self.div))
		for i in xrange(0,len(layer)):
			l1 = []
			l1.append(layer[i])
			self.recursive(l1,0,len(layer)-i)
		for i in xrange(0,len(self.points)):
			s = sum(self.points[i])
			self.points[i].append(1.0-s)
		self.points = np.asarray(self.points)
	
