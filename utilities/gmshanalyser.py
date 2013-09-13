#!python
from sys import argv,exit
from gmshtool import parse_file
import os
import numpy as np

def find_center(gmsh):
  nodes_y_max = max(gmsh.Nodes.nodes,key = lambda node: node[1])[1]
  nodes_y_min = min(gmsh.Nodes.nodes,key = lambda node: node[1])[1]
  return (nodes_y_max + nodes_y_min)/2.

def find_element_minmax(gmsh,element_index,coord_index):
  """
  Returns a tuple of min and max 'coord_index' coordinate values
  of the element with index 'element_index' from 'gmsh' object instance
  """
  el = gmsh.Elements.elements[element_index]
  def coord_value(idx):
    return gmsh.Nodes.nodes[idx][coord_index]
  el_min = gmsh.Nodes.nodes[min(el,key = coord_value)][coord_index]
  el_max = gmsh.Nodes.nodes[max(el,key = coord_value)][coord_index]
  return (el_min,el_max,)

def find_center_elements(gmsh,center):
  """
  Returns an array of element indicices which are on center
  """
  def element_in(el_minmax):
    return center > el_minmax[0] and center < el_minmax[1]
  elements = filter(lambda i: element_in(find_element_minmax(gmsh,i,1)),
      range(len(gmsh.Elements.elements)))
  return elements     


def mean_stress(gmsh, iteration, elements):
  """
  Calculates the mean stress tensor by given list of elements indexes,
  Gmsh with data and iteration number to take stresses from
  """
  stress_list = [gmsh.ElementData[iteration].element_data[i] for i in elements]
  stress = [[] for i in range(len(stress_list[0]))]
  for s in stress_list:
    for i in range(len(s)):
      stress[i].append(s[i])                   
  return [np.mean(x) for x in stress]

def print_strain_stress_graph(gmsh, component, elements):
  """
  Prints the stress-strain data for given Gmsh object and given
  stress component i and given elements ids list
  """
  for i in range(len(gmsh.ElementData)):
    strain = gmsh.ElementData[i].step_percents
    stress = mean_stress(gmsh,i,elements)
    print("%f %f" % (strain, stress[component]))
          
def run(filename):
  #print(os.path.splitext(filename)[0])
  gmsh = parse_file(filename)
  print("Number of iterations passed: %d" % (len(gmsh.ElementData)-1))
  #print("Number of nodes in 1st iteration: %d" % len(gmsh.NodeData[0].node_data))
  #print("Number of elements in 1st iteration: %d" % len(gmsh.ElementData[0].element_data))
  #print("2st iteration 1st node offsets on step %f: %s" % (gmsh.ElementData[1].step_percents, gmsh.ElementData[1].element_data[0]))
  center_y = find_center(gmsh)
  print "Center y coordinate is:", center_y
  elements = find_center_elements(gmsh,center_y)
  print "Number of elements in center:", len(elements)
  print "Their indiceis:", elements
  print_strain_stress_graph(gmsh, 4, elements)

  
if len(argv) > 1:
  run(argv[1])
else:
  print("Syntax: %s results.msh" % argv[0])
