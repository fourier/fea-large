#!/usr/bin/python
from sys import argv,exit
import re
import string

def create_filename(increment):
		return "stresses%d.txt" % increment

def avg_stress_yy(file):
		lines = file.readlines()
		lines = lines[1:]
		avg = 0.0
		cnt = 0
		for line in lines:
				values = [float(v) for v in line.split()]
				avg = avg + values[2]
				cnt = cnt + 1
		return avg/cnt

graphic = []
iterations = []
for increment in range(1,10000):
    name = create_filename(increment)
    try:
        f = open(name,"r")
        f.close()
    except:
        break
    name = create_filename(increment)
    f = open(name,"r")
    graphic.append(avg_stress_yy(f))
    f.close()

f = open("calculated5.txt","w+")
stress = 1
for v in graphic:
		stress = stress + 0.05/6.
		f.write("%f %f\n" % (stress,v))
f.close()

