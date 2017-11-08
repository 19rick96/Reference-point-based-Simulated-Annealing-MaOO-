import numpy as np
import copy
from math import *
import random
import time
import operator
import csv
from numpy import linalg as LA
import wfg
import PyGMO

def WFG(input_arr,n_obj,fn_name):
	position = 2*(n_obj-1)
	distance = 20
	func,bounds = wfg._wfg_problem_instance(fn_name,distance,position,n_obj)
	dec_vec = np.zeros(len(input_arr))
	for i in xrange(0,len(input_arr)):
		dec_vec[i] = bounds[0][i] + (input_arr[i]*(bounds[1][i]-bounds[0][i]))
	return np.asarray(func(dec_vec))

#################################################################################################################################

def f1(input_arr):
	f1 = input_arr[0]
	f2 = (1 + input_arr[1])/f1
	return [f1,f2]

def ZDT1(input_arr):
	f1 = input_arr[0]
	s = 0.0
	for i in range(1,len(input_arr)):
		s = s + input_arr[i]
	g = 1.0 + ((9.0*s)/(len(input_arr)-1.0))
	f2 = g * (1.0 - math.sqrt(f1/g))
	return [f1,f2]
	
def ZDT2(input_arr):
	f1 = input_arr[0]
	s = 0.0
	for i in range(1,len(input_arr)):
		s = s + input_arr[i]
	g = 1.0 + ((9.0*s)/(len(input_arr)-1.0))
	f2 = g*(1.0 - ((f1/g)**2))
	return [f1,f2]

def ZDT3(input_arr):
	f1 = input_arr[0]
	s = 0.0
	for i in range(1,len(input_arr)):
		s = s + input_arr[i]
	g = 1.0 + ((9.0*s)/(len(input_arr)-1.0))
	f2 = g*(1.0 - math.sqrt(f1/g) - ((f1/g)*math.sin(10.0*3.14*f1)))
	return [f1,f2]

def ZDT4(input_arr):
	f1 = input_arr[0]
	s = 0.0
	for i in range(1,len(input_arr)):
		s = s + ((input_arr[i]**2) - (10.0*math.cos(4*3.14*input_arr[i])))
	g = 1.0 + (10.0*(len(input_arr)-1.0)) + s
	f2 = g*(1.0 - ((f1/g)**2))
	return [f1,f2]

def ZDT6(input_arr):
	f1 = 1.0 - (math.exp(-4.0*input_arr[0])*((math.sin(6.0*3.14*input_arr[0]))**6))
	s = 0.0
	for i in range(1,len(input_arr)):
		s = s + input_arr[i]
	g = 1.0 + 9.0*((s/(len(input_arr)-1.0))**0.25)
	f2 = 1.0 - ((f1/g)**2)
	return [f1,f2]

##################################################################################################################################

def DTLZ1(input_arr,n_obj):
	n_var = len(input_arr)
	k = n_var - n_obj + 1
	f = [0.0]*n_obj
	g = 0.0
	for i in range(n_var-k,n_var):
		g = g + ((input_arr[i] - 0.5)**2) - math.cos(20*math.pi*(input_arr[i] - 0.5))
	g = 100*(k+g)
	for i in range(1,n_obj+1):
		s = 0.5 * (1 + g)
		j = n_obj - i
		while j >= 1 :
			j = j - 1
			s = s * input_arr[j]
		if i > 1:
			s = s * (1 - input_arr[n_obj - i])
		f[i-1] = s
	return np.asarray(f)
		
def DTLZ2(input_arr,n_obj):
	n_var = len(input_arr)
	k = n_var - n_obj + 1
	f = [0.0]*n_obj
	g = 0.0
	for i in range(n_var-k,n_var):
		g = g + ((input_arr[i] - 0.5)**2)
	for i in range(1,n_obj+1):
		s = (1.0 + g)
		j = n_obj - i
		while j >= 1 :
			j = j - 1
			s = s * math.cos(math.pi*input_arr[j]*0.5)
		if i > 1:
			s = s * math.sin(input_arr[n_obj - i]*math.pi/2)
		f[i-1] = s
	return np.asarray(f)
	
def DTLZ3(input_arr,n_obj):
	n_var = len(input_arr)
	k = n_var - n_obj + 1

	g = 0
	for i in xrange(n_obj-1,n_var):
		g += (((input_arr[i]-0.5)**2) - cos(20.0*pi*(input_arr[i]-0.5)))
	g = (k+g)*100
	out = np.zeros(n_obj)

	for m in xrange(0,n_obj):
		product = (1+g)
		i = 0
		while((i+m)<=n_obj-2):
			product *= cos(input_arr[i]*pi/2)	
			i += 1
		if m>0:
			product *= sin(input_arr[i]*pi/2)
		out[m] = product
	return np.asarray(out)

def DTLZ4(input_arr,n_obj,a=100):
	n_var = len(input_arr)
	k = n_var - n_obj + 1
	
	g = 0
	for i in xrange(n_obj-1,n_var):
		g += ((input_arr[i]-0.5)**2)
	out = np.zeros(n_obj)
	for m in xrange(0,n_obj):
		product = (1+g)
		i=0
		while((i+m) <= n_obj-2):
			product *= cos((input_arr[i]**a)*pi/2)
			i += 1
		if m>0 :
			product *= sin((input_arr[i]**a)*pi/2)
		out[m] = product
	return np.asarray(out)

def evaluate(func,decision,n_obj):
	if "wfg" in func:
		return WFG(decision,n_obj,func)
	elif func == "DTLZ1":
		return DTLZ1(decision,n_obj)
	elif func == "DTLZ2":
		return DTLZ2(decision,n_obj)
	elif func == "DTLZ3":
		return DTLZ3(decision,n_obj)
	elif func == "DTLZ4":
		return DTLZ4(decision,n_obj)
	elif func == "ZDT1":
		return ZDT1(decision)
	elif func == "ZDT2":
		return ZDT2(decision)
	elif func == "ZDT3":
		return ZDT3(decision)
	elif func == "ZDT4":
		return ZDT4(decision)
	else:
		return ZDT6(decision)

