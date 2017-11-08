import numpy as np
import copy
from math import *
import random
import time
import operator
import csv
from numpy import linalg as LA
from rsa import *

optimizer = RSA()

optimizer.define_exp(False,False,0.99,31,91,300,12,"SBX","diff",0.75,"DTLZ4",12,3,0.5,0.2,30,20,0.0,0.0,0.25)
optimizer.optimize()
igd = optimizer.evaluate_opt()
print igd
