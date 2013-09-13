#!/usr/bin/python

import xml.dom.minidom
import os

class BasePlane:
	def __init__(self):
		# array of unique nodes 
		self.points = []
		# index triangle array		
		self.triangles_indexed = []


def export_matlab(baseplane,boundary,symmetry,fnname):
	# open a file
	filename = fnname + ".m"
	mfile = open(filename,"w+")		
	# write header
	mfile.write("function [nodes,elements,symmetry,boundary]=")
	mfile.write(fnname)
	mfile.write("\n")
	# write nodes array	
	line = "nodes = [\n"
	mfile.write(line)
	for point in baseplane.points:
		line = "%f,%f;\n" % (point.x, point.y)
		mfile.write(line)
	mfile.write("];\n")
	# write elements array - triangles from plane triangulation with indexes	
	mfile.write("elements = [\n")
	for trn in baseplane.triangles_indexed:
		line = ""
		# all indexes for MATLAB/OCTAVE should start from 1
		for i in range(len(trn)-1):
			line = line + str(trn[i]+1) + ","
		line = line + str(trn[len(trn)-1]+1) + ";\n"
		mfile.write(line)
	mfile.write("];\n")
	# write array of boundary conditions
	mfile.write("boundary = [\n")
	for bnd in boundary:
		# all indexes for MATLAB/OCTAVE should start from 1
		line = ""
		if len(bnd) == 3:
			line = "%d,%f,%f,3;\n" % (bnd[0]+1,bnd[1],bnd[2])
		elif len(bnd) == 4:
			line = "%d,%f,%f,%d;\n" % \
			    (bnd[0]+1,bnd[1],bnd[2],int(bnd[3]))
		mfile.write(line)
	mfile.write("];\n")
	# export empty symmetry condition
	mfile.write("symmetry = [");
	line = ""
	for sym in symmetry:
		if line == "":
			line = "%d" % int(sym+1)
		else:
			line = "%s,%d" %(line,sym+1)

	mfile.write(line);
	mfile.write("];\n")	
	mfile.close()

def export_geometry2d(baseplane,boundary,symmetry,fnname):
	# open a file
	filename = fnname + ".node"
	mfile = open(filename,"w+")		
	# write header
	line = "%d %d %d %d\n\n" % (len(baseplane.points),1,0,len(baseplane.triangles_indexed))
	mfile.write(line)
	# write boundary type
	line = "SymAndDisplacements %d %d\n\n" % (len(symmetry), len(boundary))
	mfile.write(line)
	# write nodes array	
	i = 0
	for point in baseplane.points:
		line = "%d %f %f 0\n" % (i, point.x, point.y)
		mfile.write(line)
		i = i+1
	mfile.write("\n")
	# write segments
	line = "0 1\n\n\n"
	mfile.write(line)
	
	# write elements array - triangles from plane triangulation with indexes	
	for trn in baseplane.triangles_indexed:
		line = "%d %d %d %d %d %d\n" % \
						(trn[0],trn[1],trn[2],trn[3],trn[4],trn[5])
		mfile.write(line)
	mfile.write("\n")

	# write symmetry array
	for i in range(len(symmetry)):
		line = "%d " % symmetry[i]
		mfile.write(line) 
	mfile.write("\n\n")

	# write array of boundary conditions
	for bnd in boundary:
		# all indexes for MATLAB/OCTAVE should start from 1
		line = ""
		if len(bnd) == 3:
			line = "%d %f %f 3\n" % (bnd[0],bnd[1],bnd[2])
		elif len(bnd) == 4:
			line = "%d %f %f %d\n" % \
			    (bnd[0],bnd[1],bnd[2],int(bnd[3]))
		mfile.write(line)
	mfile.close()

def export_flagshyp(plane,boundary,symmetry,fnname):
	# open a file
	filename = fnname + ".dat.txt"
	mfile = open(filename,"w+")		
	# write header
	mfile.write("2d plane\n")
	mfile.write("tria6\n")
	# find prescribed(non zero) displacements
	prescribed_displacements = []
	for bnd in boundary:
		if bnd[1] != 0.0:
			prescribed_displacements.append((bnd[0],1,bnd[1]))			
		if bnd[2] != 0.0:
			prescribed_displacements.append((bnd[0],2,bnd[2]))			
	# write number of nodes
	line = "%d\n" % len(plane.points)
	mfile.write(line)
	# write nodes array	
	# format: index boundary_code x y 
	for i in range(len(plane.points)):
		point = plane.points[i]	
		# find in boundary array
		boundary_code = 0
		for bnd in boundary:
			if i == bnd[0]:
				boundary_code = 3 # (x,y) prescribed
				break
		line = "%d %d %f %f\n" % (i+1, boundary_code, point.x, point.y)
		mfile.write(line)
	# write elements array - triangles from plane triangulation with indexes	
	# number of elements
	line = "%d\n" % len(plane.triangles_indexed)
	mfile.write(line)
	for i in range(len(plane.triangles_indexed)):
		trn = plane.triangles_indexed[i]
		material_number = 1
#		line = "%d %d %d %d %d %d %d %d\n" % \
#						(i+1, material_number, \
#						 trn[0]+1,trn[5]+1,trn[2]+1,trn[4]+1,trn[1]+1,trn[3]+1)
		line = "%d %d %d %d %d %d %d %d\n" % \
						(i+1, material_number, \
						 trn[0]+1,trn[1]+1,trn[2]+1,trn[3]+1,trn[4]+1,trn[5]+1)
		mfile.write(line)
	# number of different materials
	line = "%d\n" % 1
	mfile.write(line)
	# Material data: number of material, material type
	# material type 1 - plane strain compressible neo-Hookean
	# material type 6 - plane stress incompressible neo-Hookean
	line = "%d %d\n" % (1,1)
	mfile.write(line)
	# Material constants
	# constants - r,m,h - density(rho), Lame constant(mu), thikness(h)
	line = "%f %f %f\n" % (1,100,100)
	mfile.write(line)
	# Loads
	# Number of loaded nodes, number of nodes with prescribed displacements,
	# number of line or surface elements with applied pressure, 
	# gravity vector(x,y)
	line = "%d %d %d %f %f\n" % \
		(0,len(prescribed_displacements),0,0,0)
	mfile.write(line);
	# Prescribed displacements
	for disp in prescribed_displacements:
		line = "%d %d %f\n" % ( disp[0]+1, disp[1], disp[2])
		mfile.write(line)
	# Solution control parameters
	# nincr: number of load/displacement increments,
	nincr = 1
	# xlmax: max. value of load scaling parameter,
	xlmax = 1
	# dlamb: load parameter increment, 
	dlamb = 1
	# miter: maximum allowed no. of iteration per increment,
	miter = 50
	# cnorm: convergence tolerance,
	cnorm = 1e-7
	# searc: line search parameter(if 0.0 not in use),
	searc = 0.0
	# arcln: arc length parameter (if 0.0 not in use)
	arcln = 0.0
	line = "%d %f %f %d %e %f %f\n" % \
		(nincr,xlmax,dlamb,miter,cnorm,searc,arcln)
	mfile.write(line)
	# Finalize file
	mfile.close()

def export_flagshyp3d(nodes,elements,boundary,fnname):
	# open a file
	filename = fnname + ".dat.txt"
	mfile = open(filename,"w+")		
	# write header
	mfile.write("3d compressible neo-hookean task\n")
	mfile.write("tetr10\n")
	# find prescribed(non zero) displacements
	# boundary format:
	# boundary[0] - node number (1-based)
	# boundary[1]-boundary[3] - displacements in 1,2,3 directions(x,y,z)
	# boundary[4] is one of the following:
	# EFree = 0; % free;
  # EPrescribedX = 1; % : x prescribed 
  # EPrescribedY = 2; % : y prescribed 
  # EPrescribedXY = 3; % : x, y prescribed 
  # EPrescribedZ = 4; % : z prescribed
  # EPrescribedXZ = 5; % : x, z prescribed
  # EPrescribedYZ = 6; % : y, z prescribed
  # EPrescribedXYZ = 7;  % : x, y, z prescribed
	prescribed_displacements = []
	for bnd in boundary:
		if bnd[1] != 0.0:
			prescribed_displacements.append((bnd[0],1,bnd[1]))			
		if bnd[2] != 0.0:
			prescribed_displacements.append((bnd[0],2,bnd[2]))			
		if bnd[3] != 0.0:
			prescribed_displacements.append((bnd[0],3,bnd[3]))
	# write number of nodes
	line = "%d\n" % len(nodes)
	mfile.write(line)
	# write nodes array	
	# format: index boundary_code x y 
	for i in range(len(nodes)):
		# find in boundary array
		boundary_code = 0
		for bnd in boundary:
			if i == bnd[0]-1:
				boundary_code = bnd[4]
				break
		line = "%d %d %f %f %f\n" % (i+1, boundary_code, \
				nodes[i][0], nodes[i][1], nodes[i][2] )
		mfile.write(line)
	# write elements array - triangles from plane triangulation with indexes	
	# number of elements
	line = "%d\n" % len(elements)
	mfile.write(line)
	for i in range(len(elements)):
		trn = elements[i]
		material_number = 1
		line = "%d %d %d %d %d %d %d %d %d %d %d %d\n" % \
						(i+1, material_number, \
						 trn[0],trn[1],trn[2],trn[3],trn[4],trn[5],trn[6],trn[7],\
						 trn[8],trn[9])
		mfile.write(line)
	# number of different materials
	line = "%d\n" % 1
	mfile.write(line)
	# Material data: number of material, material type
	# material type 1 - plane strain compressible neo-Hookean
	# material type 6 - plane stress incompressible neo-Hookean
	line = "%d %d\n" % (1,1)
	mfile.write(line)
	# Material constants
	# constants - r,m,h - density(rho), Lame constant(mu), thikness(h)
	line = "%f %f %f\n" % (1,100,100)
	mfile.write(line)
	# Loads
	# Number of loaded nodes, number of nodes with prescribed displacements,
	# number of line or surface elements with applied pressure, 
	# gravity vector(x,y)
	line = "%d %d %d %f %f\n" % \
		(0,len(prescribed_displacements),0,0,0)
	mfile.write(line);
	# Prescribed displacements
	for disp in prescribed_displacements:
		line = "%d %d %f\n" % ( disp[0]+1, disp[1], disp[2])
		mfile.write(line)
	# Solution control parameters
	# nincr: number of load/displacement increments,
	nincr = 1
	# xlmax: max. value of load scaling parameter,
	xlmax = 1
	# dlamb: load parameter increment, 
	dlamb = 1
	# miter: maximum allowed no. of iteration per increment,
	miter = 50
	# cnorm: convergence tolerance,
	cnorm = 1e-7
	# searc: line search parameter(if 0.0 not in use),
	searc = 0.0
	# arcln: arc length parameter (if 0.0 not in use)
	arcln = 0.0
	line = "%d %f %f %d %e %f %f\n" % \
		(nincr,xlmax,dlamb,miter,cnorm,searc,arcln)
	mfile.write(line)
	# Finalize file
	mfile.close()

def export_xml(nodes,elements,boundary,fname):
	f = open(fname,"w+")
	document = xml.dom.minidom.Document()
	# root 'task' node
	task = document.createElement("task")
	document.appendChild(task)
	# add material model details
  # <model name="A5">
	model = document.createElement("model")
	model.setAttribute("name","A5")
	# <model-parameters lambda="100" mu="100"/>
	model_parameters = document.createElement("model-parameters")
	model_parameters.setAttribute("lambda","100")
	model_parameters.setAttribute("mu","100")
	# </model-parameters>
	model.appendChild(model_parameters)
	# </model>
	task.appendChild(model)
  # <solution modified-newton="yes" task-type="CARTESIAN3D" load-increments-count="100" desired-tolerance="1e-8">
	solution = document.createElement("solution")
	solution.setAttribute("modified-newton","yes")
	solution.setAttribute("task-type","CARTESIAN3D")
	solution.setAttribute("load-increments-count","100")
	solution.setAttribute("desired-tolerance","1e-8")
	# 	<element-type name="TETRAHEDRA10" nodes-count="10" gauss-nodes-count="5"/>
	element_type = document.createElement("element-type")
	element_type.setAttribute("name","TETRAHEDRA10")
	element_type.setAttribute("nodes-count","10")
	element_type.setAttribute("gauss-nodes-count","5")
	solution.appendChild(element_type)
	# <line-search max="5"/>
	line_search = document.createElement("line-search")
	line_search.setAttribute("max","5")
	solution.appendChild(line_search)
	# <arc-length max="5"/>
	arc_length = document.createElement("arc-length")
	arc_length.setAttribute("max","0")
	solution.appendChild(arc_length)
  # </solution>
	task.appendChild(solution)
  # <input-data>
	input_data = document.createElement("input-data")
	# <geometry>
	geometry = document.createElement("geometry")
	# <nodes count="200">
	geom_nodes = document.createElement("nodes")
	geom_nodes.setAttribute("count",str(len(nodes)))
	for i in range(len(nodes)):
		#  <node x="0.000000" y="1.000000" z="0.000000"/>
		node = nodes[i]
		nod = document.createElement("node")
		nod.setAttribute("id",str(i))
		nod.setAttribute("x",str(node[0]))
		nod.setAttribute("y",str(node[1]))
		nod.setAttribute("z",str(node[2]))
		geom_nodes.appendChild(nod)
	# </nodes>							
	geometry.appendChild(geom_nodes)
	# <elements count="100">
	elems = document.createElement("elements")
	elems.setAttribute("count",str(len(elements)))
	for i in range(len(elements)):
		# 	<element id="0" node1="69" node2="70" node3="22" node4="82" node5="89" node6="90" node7="91" node8="92" node9="93" node10="94"/>
		el = elements[i]
		elem = document.createElement("element")
		elem.setAttribute("id",str(i))
		for e in range(len(el)):
			nod_id = "node%d" % (e + 1)
			elem.setAttribute(nod_id,str(el[e]-1))
			
		# </element>
		elems.appendChild(elem)
	#   </elements>
	geometry.appendChild(elems)  
	# 	</geometry>
	input_data.appendChild(geometry)
	# <boundary-conditions>
	bcs = document.createElement("boundary-conditions")
	# <prescribed-displacements count="10">
	prs = document.createElement("prescribed-displacements")
	prs.setAttribute("count",str(len(boundary)))
	for i in range(len(boundary)):
		# <presc-node id="1" x="0" y="0" z="0" type="7"/>
		bnd = boundary[i]
		bnd_node = document.createElement("presc-node")
		bnd_node.setAttribute("id",str(i))
		bnd_node.setAttribute("node-id",str(bnd[0]-1))
		bnd_node.setAttribute("x",str(bnd[1]))
		bnd_node.setAttribute("y",str(bnd[2]))
		bnd_node.setAttribute("z",str(bnd[3]))
		bnd_node.setAttribute("type",str(bnd[4]))
		prs.appendChild(bnd_node)
  #	</prescribed-displacements>
	bcs.appendChild(prs)
	# </boundary-conditions>
	input_data.appendChild(bcs)
	# 	</input-data
	task.appendChild(input_data)
	
# boundary format:
	# boundary[0] - node number (1-based)
	# boundary[1]-boundary[3] - displacements in 1,2,3 directions(x,y,z)
	# boundary[4] is one of the following:
	# EFree = 0; % free;
  # EPrescribedX = 1; % : x prescribed 
  # EPrescribedY = 2; % : y prescribed 
  # EPrescribedXY = 3; % : x, y prescribed 
  # EPrescribedZ = 4; % : z prescribed
  # EPrescribedXZ = 5; % : x, z prescribed
  # EPrescribedYZ = 6; % : y, z prescribed
  # EPrescribedXYZ = 7;  % : x, y, z prescribed	
	#task.writexml(f)
	f.write(task.toprettyxml())
	f.close()


class SexpElement(object):
  @staticmethod
  def createElement(name):
    return SexpElement(name)
  def __init__(self,name):
    self._name = name
    self._children_list = []
  def setAttribute(self,name,value):
    self.__setattr__(name, value)
  def appendChild(self,child):
    self._children_list.append(child)
  def __str__(self):
    return self.tostring()
  def tostring(self):
    result = "(" + self._name
    for attr in self.__dict__:
      if attr != '_children_list' and attr != '_name':
        result = result + " :" + attr + " " + self.__dict__[attr]
    for elt in self._children_list:
      result = result + "\n"
      if isinstance(elt,list):
        s = "(" + str(elt[0])
        for e in elt[1:]:
          s = s + " " + str(e)
        s = s + ")"
        result = result + s
      else:
        result = result + elt.__str__()
    return result + ")"

def export_sexp(nodes,elements,boundary,fname):
	f = open(fname,"w+")
	f.write(";; -*- Mode: lisp; -*-\n")
	task = SexpElement.createElement("task")
	model = SexpElement.createElement("model")
	model.setAttribute("name","A5")
	model_parameters = SexpElement.createElement("model-parameters")
	model_parameters.setAttribute("lambda","100")
	model_parameters.setAttribute("mu","100")
	model.appendChild(model_parameters)
	task.appendChild(model)
	solution = SexpElement.createElement("solution")
	solution.setAttribute("modified-newton","yes")
	solution.setAttribute("task-type","CARTESIAN3D")
	solution.setAttribute("load-increments-count","100")
	solution.setAttribute("desired-tolerance","1e-8")
	element_type = SexpElement.createElement("element-type")
	element_type.setAttribute("name","TETRAHEDRA10")
	element_type.setAttribute("nodes-count","10")
	element_type.setAttribute("gauss-nodes-count","5")
	solution.appendChild(element_type)
	line_search = SexpElement.createElement("line-search")
	line_search.setAttribute("max","5")
	solution.appendChild(line_search)
	arc_length = SexpElement.createElement("arc-length")
	arc_length.setAttribute("max","0")
	solution.appendChild(arc_length)
	task.appendChild(solution)
	input_data = SexpElement.createElement("input-data")
	geometry = SexpElement.createElement("geometry")
	geom_nodes = SexpElement.createElement("nodes")
	for i in range(len(nodes)):
		node = nodes[i]
		geom_nodes.appendChild(node)
	geometry.appendChild(geom_nodes)
	elems = SexpElement.createElement("elements")
	for i in range(len(elements)):
		el = elements[i]
		elems.appendChild(map(lambda x: x-1, el))
	geometry.appendChild(elems)  
	input_data.appendChild(geometry)

	bcs = SexpElement.createElement("boundary-conditions")
	prs = SexpElement.createElement("prescribed-displacements")
	for i in range(len(boundary)):
		bnd = boundary[i]
		bnd_node = SexpElement.createElement("presc-node")
		bnd_node.setAttribute("node-id",str(bnd[0]-1))
		bnd_node.setAttribute("x",str(bnd[1]))
		bnd_node.setAttribute("y",str(bnd[2]))
		bnd_node.setAttribute("z",str(bnd[3]))
		bnd_node.setAttribute("type",str(bnd[4]))
		prs.appendChild(bnd_node)
	bcs.appendChild(prs)
	input_data.appendChild(bcs)
	task.appendChild(input_data)
	
	f.write(task.tostring())
	f.close()
	indent_str = "emacs -batch %s --eval '(indent-region (point-min) (point-max) nil)' -f save-buffer" % fname
	os.system(indent_str)
