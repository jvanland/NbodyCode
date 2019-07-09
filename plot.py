#!/usr/bin/python

from pylab import *
import sys

pi = 3.141592653589793

if len(sys.argv) != 2:
    print "usage: ./plot.py input_file"
    sys.exit(2)
    
filename = sys.argv[1]
f1 = open(filename, 'r')
data1 = [map(float, line.split()) for line in f1] 
c1 = array([bit[0] for bit in data1]) 
c2 = array([bit[1] for bit in data1])
c3 = array([bit[2] for bit in data1])
c4 = array([bit[3] for bit in data1])
c5 = array([bit[4] for bit in data1])
c6 = array([bit[5] for bit in data1])
c7 = array([bit[6] for bit in data1])
f1.close
time = arange(0,len(c1),1)
radius = array([sqrt(c2[i]*c2[i]+c3[i]*c3[i]+c4[i]*c4[i]) for i in range(len(c2))])
velocity = array([(c2[i]*c5[i]+c3[i]*c6[i]+c4[i]*c7[i])/radius[i] for i in range(len(c2))])

fig1 = figure()
plot(radius,velocity,'o')

show()
