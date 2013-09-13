#!/usr/bin/python

from base import Point
from triangles import Triangle3,Triangle6
from exporter import BasePlane,export_matlab,export_geometry2d,export_flagshyp


# Program for generating 2d plane geometry triangulation and 
# simple boundary conditions
# for specified rectangles with specified number of divisions
# Example:
#
# A------------
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# |           |
# ------------B
#
# - Given rectangular geometry with base points A and B
# If we specify number of vertical bands as M and horizontal bands as N,
# where for example M = 3 and N = 4, we will get a following division to 
# rectangular areas:
#
#        
#      M=3 
#       
#   /   |   \
#
# A------------
# |   |   |   |
# |   |   |   |
# |   |   |   |  \
# |---|---|---|  
# |   |   |   |
# |   |   |   |  -  
# |   |   |   |     N=4
# |---|---|---|
# |   |   |   |  /
# |   |   |   | 
# |   |   |   |
# |---|---|---|  /
# |   |   |   |
# |   |   |   |
# |   |   |   |
# ------------B
#
#
# So the resulting triangulation will be
#
#
# A------------
# |\ /|\ /|\ /|
# | \ | \ | \ |
# |/ \|/ \|/ \|  
# |---|---|---|  
# |\ /|\ /|\ /|
# | \ | \ | \ |    
# |/ \|/ \|/ \|    
# |---|---|---|
# |\ /|\ /|\ /|  
# | \ | \ | \ | 
# |/ \|/ \|/ \|
# |---|---|---|  
# |\ /|\ /|\ /|
# | \ | \ | \ |
# |/ \|/ \|/ \|
# ------------B
#

class PlaneTriangulation(BasePlane):
	def __init__(self,
				 classname,# triangle class name
				 A = Point(1,1), # upper-left point
				 B = Point(4,7), # lower-right point
				 M = 1, # Vertical bands
				 N = 1): # Horizontal bands
		# store template class name
		self.T_ = classname
		#############################
		### Initial geometry data ###
		#############################
		self.A = A 
		self.B = B
		self.M = M
		self.N = N
		# array of unique nodes 
		self.points = []
		# index triangle array		
		self.triangles_indexed = []
		self.generate()		
	
	def triangle(self):
		generator = self.T_ + "()"
		return eval(generator)
	############################
	### Starting the program ###
	############################
	def generate(self):
		size_x = (self.B.x - self.A.x)*1.0
		size_y = (self.B.y - self.A.y)*1.0

		step_x = size_x/self.M
		step_y = size_y/self.N

		# Prepare nodes and triangles in index form
		for i in range(self.M):
		    for j in range(self.N):
				rect_nodes = []
				# Get rectangle coords together 
				# with center of rect
				# --------
				# |An  Cn|
				# | \  / |
				# |  En  |
				# | /  \ |
				# |Dn  Bn|
				# --------
				An = Point(self.A.x + i*step_x, \
					   self.A.y + j*step_y)
				Bn = Point(self.A.x + (i+1)*step_x, \
					   self.A.y + (j+1)*step_y)
				Cn = Point(self.A.x + (i+1)*step_x, \
					   self.A.y + j*step_y)
				Dn = Point(self.A.x + i*step_x, \
					   self.A.y + (j+1)*step_y)
				En = Point(self.A.x + (i+0.5)*step_x, \
					   self.A.y + (j+0.5)*step_y)

				# Create 4 triangles in rectangle
				# --------
				# |\ 1 / |
				# | \ /  |
				# |4 /\ 2|
				# | /  \ |
				# |/ 3  \|
				# --------
				triangle1 = self.triangle()
				triangle2 = self.triangle()
				triangle3 = self.triangle()
				triangle4 = self.triangle()

				triangle1.reset(An,Cn,En)
				triangle2.reset(Bn,En,Cn)
				triangle3.reset(Dn,En,Bn)				
				triangle4.reset(An,En,Dn)

				# append triangles to triangles array
				triangles = []
				triangles.append(triangle1)
				triangles.append(triangle2)
				triangles.append(triangle3)
				triangles.append(triangle4)
				
				# Now, add points to a simple array
				# loop through all 4 triangles
				# generate fake triangle( array of 6 indexes for 6-nodes 
				# triangle element
				# and add it to array of triangles with indexes
				for k in range(len(triangles)):
					self.triangles_indexed.append(\
						range(len(triangles[k].p))) 
														
					for p in range(len(triangles[k].p)): # loop through all 
														 # 6 nodes of element
						# index of node found in global unique array of points
						found = -1  
						for l in range(len(self.points)): # loop through all 
 													      							# elements of global 
																							# array of unique nodes
							if triangles[k].p[p] == self.points[l]: # found element in 
																											# global array of 
																											#	unique nodes 
								found = l # set index 'found'
								break
						if found == -1: # if node is not found in global 
														# array of unique nodes
							self.points.append(triangles[k].p[p]) # add it 			                  	   
							found = len(self.points)-1 # it has now index of 
																				 # last element in uniq array
						self.triangles_indexed[-1][p] = found # at last, set current 
																									# index of element to
																									# found value

def create_boundary_conditions(planeTriangulation,offset,fixed_x):
	upper = planeTriangulation.A.y
	lower = planeTriangulation.B.y
	boundary = []
	index = 0
	condition_type = 3
	if fixed_x == False:
		condition_type = 2
	for i in range(len(planeTriangulation.points)):
		if planeTriangulation.points[i].y == upper:
			boundary.append((i,0,0,condition_type))
		if planeTriangulation.points[i].y == lower:
			boundary.append((i,0,offset,condition_type))
	return boundary

def create_symmetry_conditions(planeTriangulation):
	left = planeTriangulation.A.x
	symmetry = []
	index = 0
	for i in range(len(planeTriangulation.points)):
		if planeTriangulation.points[i].x == left:
			symmetry.append(i)
	return symmetry



def export_to_flagshypfile(plane,offset,fnname):
	# generate array of boundary nodes with displacements conditions
	boundary = create_boundary_conditions(plane,offset,False)
	# only needed to create symmetry conditions
	symmetry = create_symmetry_conditions(plane)
	export_flagshyp(plane,boundary,symmetry,fnname)

def export_to_mfile(plane,offset,fnname):
	# generate array of boundary nodes with displacements conditions
	boundary = create_boundary_conditions(plane,offset,False)
	# only needed to create symmetry conditions
	symmetry = create_symmetry_conditions(plane)
	export_matlab(plane,boundary,symmetry,fnname)


def export_to_geometry2dfile(plane,offset,fnname):
	# generate array of boundary nodes with displacements conditions
	boundary = create_boundary_conditions(plane,offset,False)
	symmetry = create_symmetry_conditions(plane)
	export_geometry2d(plane,boundary,symmetry,fnname)


plane = PlaneTriangulation("Triangle6",Point(1,1),Point(4,7),2,4)

export_to_mfile(plane,0.05,"matlab_geometry6")
export_to_geometry2dfile(plane,0.05,"matlab_geometry")
export_to_flagshypfile(plane,0.05,"bonet")
