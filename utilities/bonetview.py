#!/usr/bin/python
from base import Point
from triangles import Triangle3,Triangle6
from sys import exit 
import re

def wrong_format():
	exit("Wrong file format!")
def wrong_element():
	exit("Wrong element type!")
def wrong_files():
	exit("Results and source files contain different sets of elements!")			

def getline(lines,index):
	if index >= len(lines):
		exit("Out of bounds")
	return [index + 1, lines[index]]

def check_nodes(strings,issource):
	size_check = 0
	if issource:
		size_check = 4			
	else: # results
		size_check = 6
	if len(strings) != size_check:
		print "Nodes: len(strings) = %d, but required %d" % \
				(len(strings),size_check)				
		wrong_format()

def check_triangles(strings):
	size_check = 8
	if len(strings) != size_check:
		print "Triangles: len(strings) = %d, but required %d" % \
				(len(strings),size_check)				
		wrong_format()

def check_stresses(strings):
	size_check = 4
	if len(strings) != size_check:
		print "Stresses: len(strings) = %d, but required %d" % \
				(len(strings),size_check)				
		wrong_format()

def read_section(lines,startindex,issource):
	newstartindex = 0	
	triangles = []
	result = 1
	# section size must be at least 7 lines!!!
	if startindex >= len(lines)-7:
		result = 0			
		return [result,triangles,newstartindex]	
	triangles = []
	# start parsing section
	[newstartindex,line] = getline(lines,startindex)
	# For results first section line must contain 'at increment'
	if issource == 0:
		m = re.search("at increment",line)
		if m == None:
			# does not contain necessary string
			wrong_format()
	# Second line must contain element type
	[newstartindex,line] = getline(lines,newstartindex)
	m = re.match("^tria6",line)
	if m == None:
		wrong_element()
				
	# next line contain number of nodes 
	[newstartindex,line] = getline(lines,newstartindex)
	nodes_count = int(line)
	if nodes_count > len(lines)-newstartindex:
			wrong_format()
	# start reading nodes
	points = []
	for l in range(nodes_count):
		[newstartindex,line] = getline(lines,newstartindex)
		strings = line.split()
		check_nodes(strings,issource)
		x = float(strings[2])
		y = float(strings[3])
		points.append(Point(x,y))
	# start reading triangles 
	[newstartindex,line] = getline(lines,newstartindex)
	triangles_count = int(line)
	for l in range(triangles_count):
		[newstartindex,line] = getline(lines,newstartindex)
		strings = line.split()
		check_triangles(strings)
		vals = [int(i) for i in strings]
		point1 = points[vals[0+2]-1]
		point2 = points[vals[2+2]-1]
		point3 = points[vals[4+2]-1]
		triangle = Triangle3(point1, point2, point3)
		del triangle.deformed
		triangle.deformed = []
		triangle.deformed.append(point1)
		triangle.deformed.append(point2)
		triangle.deformed.append(point3)
		# triangle formed at last, add it to
		# triangle array
		triangles.append(triangle)
	
	# Now if it is not source section,
	if issource == 0 :
		# we need t extract data from gauss points (3 points per element)
		if newstartindex > len(lines) - triangles_count*3:
			wrong_format()
		for l in range(triangles_count):
			vals = []
			for k in range(3)	:
				[newstartindex,line] = getline(lines,newstartindex)
				strings = line.split()
				check_stresses(strings)
				vals_temp = [float(i) for i in strings]
				vals.append(vals_temp)
			g = [[vals[0][0], vals[0][1], vals[0][2],vals[0][3]],
					 [vals[1][0],	vals[1][1], vals[1][2],vals[1][3]],
					 [vals[2][0], vals[2][1], vals[2][2],vals[2][3]]]
			triangles[l].set_values_in_nodes(g)
	return [result, triangles, newstartindex]

def import_file(fname,source_fname):
	file = open(fname,"r")
	lines = file.readlines()
	file.close()
	###################
	## read source file
	###################
	# set source flag in read_section function
	issource = 1 
	# set index to 0
	startindex = 0
	# result flag - needed for future parsing of
	# results file with several iterations
	result = 1
	# array of triangles from input geometry
	triangles_source = []
	if source_fname != "":
		sfile = open(source_fname,"r")
		slines = sfile.readlines()
		sfile.close()
		[result, triangles_source, startindex] = \
			read_section(slines,startindex,issource)
	########################
	## now read results file
	########################
	triangles = []
	triangles_safe = []
	# read every iteration until end of file
	# section means iteration in flagshyp results file
	result = 1
	startindex = 0
	# section will be read from output file, not from source file
	issource = 0
	while result == 1:
		triangles_safe = triangles # save triangles from prev section
		[result, triangles, startindex ] = \
			read_section(lines,startindex,issource)	# read new section		
	
	#############	#############################
	## Form correct triangles array
	###########################################
	## triangles_source - source file triangles
	## triangles_safe - results file triangles
	###########################################
	
	# if source file was in command line
	if len(triangles_source) > 0:
		if len(triangles_source) != len(triangles_safe):	
			wrong_files()
		for i in range(len(triangles_source)):
			for j in range(3):
				triangles_safe[i].p[j] = triangles_source[i].p[j]

	return triangles_safe		

