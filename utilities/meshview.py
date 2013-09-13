#!/usr/bin/python
from base import Point
from triangles import Triangle3,Triangle6

def wrong_format():
	print "Wrong file format!"
	return []	

def import_file(fname):
	file = open(fname,"r")
	lines = file.readlines()
	file.close()
	count_args = len(lines[0].split(','))
	if count_args != 17:
		wrong_format()
	for line in lines:
		if len(line.split(',')) != 17:
			wrong_format()
	# all ok, file can be processed
	triangles = []
	for line in lines:
		vals = [float(i) for i in line.split(',')]
		triangle = Triangle3(Point(vals[0],vals[1]),Point(vals[2],vals[3]),\
			Point(vals[4],vals[5]))
		g = [[vals[14],vals[14],vals[14]],
			[vals[15],vals[15],vals[15]],
			[vals[16],vals[16],vals[16]]]
		triangle.set_values_in_nodes(g)
		del triangle.deformed
		triangle.deformed = []
		triangle.deformed.append(Point(vals[6],vals[7]))
		triangle.deformed.append(Point(vals[8],vals[9]))
		triangle.deformed.append(Point(vals[10],vals[11]))
		triangles.append(triangle)
	return triangles
