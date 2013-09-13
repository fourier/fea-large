#!/usr/bin/python

import string


# Base class for all finite elements



# Helper function which calculates factorials of numbers
def factorial(n):
	if n == 0:
		return 1
	else:
		return factorial(n-1)*n;


def isempty(line):
    return len(line) == 0 or line.isspace()

def skipwhitespaces(line):
    for i in range(len(line)):
        if line[i] == ' ' or line[i] == '\t':
            continue
        else:
            break
    return line[i:len(line)]

def iscomment(line,comment = '#'):
    line = skipwhitespaces(line)
    return line.startswith(comment)

def removeempty(lines):
    newlines = []
    for line in lines:
        if len(line) != 0 and line != ' ' \
                and line != '\t' and line != '\n':
            newlines.append(line)
    return newlines

def removecomments(lines,comment = '#'):
	newlines = []
	for line in lines:
		if not iscomment(line,comment):
			newlines.append(line)
	return newlines

def read_and_clear_lines(filename,comment = '#'):
	file = open(filename,"r")
	lines = file.readlines()
	file.close()
	lines = removeempty(lines)
	lines = removecomments(lines,comment)
	return lines

def unique(values):
	set_values = set(values)
	values = [i for i in set_values]
	return values

# parse filename and extract a name w/o extension
def parse_filename(filename):
	base_name = filename.split('.')[0]
	return base_name



# Class which represents simple 2-d plane points
# It could be either x-y orthogonal coordinates or r-z cylindrical
class Point:
	def __init__(self,x = 0, y = 0, z = 0):
		self.x = x
		self.y = y
		self.z = z
	def __str__(self):
		s = [str(self.x),str(self.y),str(self.z)]
		return string.join(s)
	def __eq__(self, other):
		if other.x == self.x and other.y == self.y and other.z == self.z:
			return True
		return False
	def __ne__(self,other):
		return not (self == other)

class Element:
	def __init__(self,p1 = Point(),p2 = Point(),p3 = Point() ):
		self.p=[]
		self.p = []
		self.deformed = []
		self.values = []
		self.p.append(p1)
		self.p.append(p2)
		self.p.append(p3)
		self.deformed.append(p1)
		self.deformed.append(p2)
		self.deformed.append(p3)
	def __str__(self):
		str = ""
		lst = []
		for i in self.p:
			lst.append(self.p.__str__())
		str = string.join(lst)
		return str
	def N( self,i,point = Point() ):
		return 0

	def set_values_in_nodes(self, val = [[]]):
#		for i in range(len(val)):
#			if len(val[i]) != len(self.p):
#				print "base.Element::set_values_in_nodes: wrong argument g =",val[i]
#				return
		if len(self.values):
			del self.values
		self.values = val
	# Calculate value of g-function by given values in nodes
	# i - degree of approximation, by default i = 3
	# and 
	def g(self,i, point = Point()):
		if i < 0 or i > len(self.values) or len(self.values) == 0:
			print "base.Element::g:wrong argument i =",i
			return
		res = 0
		for j in range(len(self.values[i])):
			res += self.N(j,point)*(self.values[i])[j]
		return res
	
