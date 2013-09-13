#!/usr/bin/python
from sys import argv,exit
from base import Point
from triangles import Triangle3,Triangle6

import meshview
import bonetview

class Triangle:
	def __init__(self,points,triangle):
		self.n = []
		for i in range(3):
			self.n.append(-1)
		self.element = triangle
		self.find_points(points)
	def find(self,points,point):
		i = -1
		try:
			i = points.index(point)
		except ValueError:
			i = -1
		return i
	def find_points(self,points):
		for i in range(3):
			self.n[i] = self.find(points,self.element.p[i]) 
	def edges(self):
		edges_array = []
		edges_array.append((self.n1,self.n2))
		edges_array.append((self.n2,self.n3))
		edges_array.append((self.n3,self.n1))
		return edges_array

		
class DataModel:
	def __init__(self,size = (300,400)):
		self.clear_all()
		self.size = size
	def clear_all(self):
		self.size = (300,400) 
		self.points = []
		self.deformed_points = []
		self.edges = []
		self.triangles = []
		self.min = Point()
		self.max = Point()
		self.deformed_min = Point()
		self.deformed_max = Point()
	def find_point(self,points,point):
		i = -1
		try:
			i = points.index(point)
		except ValueError:
			i = -1
		return i
	def extract_triangle_nodes(self,triangle):
		for j in range(3):
			i = self.find_point(self.points,triangle.p[j])
			if i == -1:
				self.points.append(triangle.p[j])
			i = self.find_point(self.deformed_points,triangle.deformed[j])
			if i == -1:
				self.deformed_points.append(triangle.deformed[j])
	def extract_points(self,triangles3):
		siz = len(self.points)
		if siz != 0:
			print "DataModel::extract_points: len(points) != 0, len(points) ==",siz 
			return
		for triangle in triangles3:
			self.extract_triangle_nodes(triangle)
	def find_minmax_points(self):
		if len(self.points) ==0 or len(self.deformed_points) == 0:
			return
		self.min.x = self.points[0].x
		self.min.y = self.points[0].y
		self.max.x = self.points[0].x
		self.max.y = self.points[0].y
		for p in self.points:
			if p.x > self.max.x:
				self.max.x = p.x
			if p.y > self.max.y:
				self.max.y = p.y
			if p.x < self.min.x:
				self.min.x = p.x
			if p.y < self.min.y:
				self.min.y = p.y
		self.deformed_min.x = self.deformed_points[0].x
		self.deformed_min.y = self.deformed_points[0].y
		self.deformed_max.x = self.deformed_points[0].x
		self.deformed_max.y = self.deformed_points[0].y
		for p in self.deformed_points:
			if p.x > self.deformed_max.x:
				self.deformed_max.x = p.x
			if p.y > self.deformed_max.y:
				self.deformed_max.y = p.y
			if p.x < self.deformed_min.x:
				self.deformed_min.x = p.x
			if p.y < self.deformed_min.y:
				self.deformed_min.y = p.y
	def import_triangles(self,triangles):
		self.clear_all()
		self.extract_points(triangles)
		for triangle in triangles:
			self.triangles.append(Triangle(self.points,triangle))
		self.find_minmax_points()
	def import_msh(self,filename):
		triangles = meshview.import_file(filename)
		if len(triangles) == 0:
			return
		self.import_triangles(triangles)
	def import_bonet(self,filename,sourcefilename):
		triangles = bonetview.import_file(filename,sourcefilename)
		if len(triangles) == 0:
			return			
		self.import_triangles(triangles)
		dump_msh(triangles,"results.msh")

def dump_model(data):
	print "Nodes:"
	for i in range(len(data.points)):
		print "i =",i, " point =",data.points[i].x, data.points[i].y
	print "Deformed Nodes:"
	for i in range(len(data.deformed_points)):
		print "i =",i, " point =", \
			data.deformed_points[i].x, data.deformed_points[i].y

	print "Triangles:"
	for i in range(len(data.triangles)):
		for j in range(3):
			print data.triangles[i].n[j],
		print ":",
		for j in range(3):
			print "(",data.points[data.triangles[i].n[j]],")",
		print 
	print "Min =",data.min
	print "Max =",data.max 
	print "Height =",data.max.y - data.min.y,"Width =",data.max.x - data.min.x
	print "Deformed Min =",data.deformed_min
	print "Deformed Max =",data.deformed_max 
	height = (data.max.y - data.min.y)
	deformed_height = data.deformed_max.y - data.deformed_min.y
	percent = height/100.0
	percents = abs(height - deformed_height)/percent
	print "Deformation step =", percents,"%"

def print_usage():
	print "Usage:"
	print "python datamodel.py (bonet|meshview) filename [source_data_file]"
	print "In case of bonet task one should also specify source data file"
	print "with initial data."

def dump_msh(triangles,filename):
	if len(triangles) == 0:
		return
	file = open(filename,"w+")
	for trn in triangles:	
		line = "%f,%f,%f,%f,%f,%f," % \
			(trn.p[0].x,trn.p[0].y, \
			 trn.p[1].x,trn.p[1].y, \
			 trn.p[2].x,trn.p[2].y)
		file.write(line)
		line = "%f,%f,%f,%f,%f,%f,0," % \
			(trn.deformed[0].x,trn.deformed[0].y, \
			 trn.deformed[1].x,trn.deformed[1].y, \
			 trn.deformed[2].x,trn.deformed[2].y)
		file.write(line)
		sigma_xx = 1/3. * \
			( trn.values[0][0] + trn.values[1][0] + trn.values[2][0] )
		sigma_yy = 1/3. * \
			( trn.values[0][2] + trn.values[1][2] + trn.values[2][2] )
		sigma_xy = 1/3. * \
			( trn.values[0][1] + trn.values[1][2] + trn.values[2][1] )
		h = 1/3. * ( trn.values[0][3] + trn.values[1][3] + trn.values[2][3] )
		line = "%f,%f,%f,%f\n" % (sigma_xx,sigma_yy,h,sigma_xy)
		file.write(line)
	file.close()

data = DataModel()
if len(argv) > 1:
	if argv[1] == "bonet":
		source = ""
		if len(argv) > 3:
			source = argv[3]			
		data.import_bonet(argv[2],source)
		dump_model(data)
	elif argv[1] == "meshview":
		data.import_msh(argv[2])
		dump_model(data)
	else :
		print_usage()
else:
		print_usage()		
