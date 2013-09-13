#!/usr/bin/python

class element1d2p:
	def __init__(self,x1 = 0, x2 = 0):
		self.x1 = x1
		self.x2 = x2
	def __str__(self):
		s = [str(self.x1),str(self.x2)]
		return string.join(s)
	def __eq__(self, other):
		if other.x1 == self.x1 and other.x2 == self.x2:
			return True
		return False
	def __ne__(self,other):
		return not (self == other)


# function for generation of one dimension mesh
# generate 1d mesh of 2-node elements by given
# first node A, last node B, number of elements N
# to matlab file 'geometry1d.m'
def one_dimension_mesh(A,B,N):
  # define step
  step = (B-A)/N
  # defile clear nodes and elements arrays
  nodes = []
  elements = []
  # fill nodes array
  for i in range(N+1):
    nodes.append(A+step*i)
  # fill elements array
  for i in range(N):
    elements.append((i,i+1))
  # open file
  f = open("geometry1d.m","w+")
  # write header
  f.write("function [nodes,elements]=geometry\n")
  # write nodes array
  line = "nodes = [\n"
  f.write(line)
  for node in nodes:
    line = "%f;\n" % node
    f.write(line)
  f.write("];\n")
  # write elements array
  f.write("elements = [\n")
  for el in elements:
    # all MATLAB/OCTAVE arrays starts from index 1, not 0
    line = "%d,%d;\n" % (el[0]+1,el[1]+1)
    f.write(line)
  f.write("];\n")
  f.close()


one_dimension_mesh(1,7,5)
