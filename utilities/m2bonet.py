#!python
# Script for conversion of 3d-matlab geometry to flagshyp task file

from sys import argv,exit
import re
from base import read_and_clear_lines,parse_filename
from exporter import export_flagshyp3d,export_xml,export_sexp

def clear_whitespaces(line):
	res = ""
	for i in range(len(line)):
		if line[i] != ' ' and line[i] != '\t' \
				and  line[i] != '\n' and line[i] != '\r':
			res += line[i]
	return res


def extract_assignments(contents):
	assignments = re.findall('[a-zA-Z]+\w*=\[[\d,;\.]*\];',contents)	
	return assignments

# Class for importing generated geometry from matlab file
# and for exporting it to bonet file
class matlab_importer:
	def __init__(self):
		self.nodes = []
		self.triangles = []
		self.boundary = []
				
	def import_file(self,filename):
		lines = read_and_clear_lines(filename,'%')
		self.trim_header(lines)
		contents = self.generate_plain_contents(lines)
		assignments = extract_assignments(contents)
		for assignment in assignments:
			self.evaluate_assignment(assignment)
		if hasattr(self,'nodes'):
			print 'size of self.nodes = %d' % len(self.nodes)
		if hasattr(self,'elements'):
			print 'size of self.elements = %d' % len(self.elements)
		if hasattr(self,'boundary'):
			print 'size of self.boundary = %d' % len(self.boundary)
		if hasattr(self,'symmetry'):
			print 'size of self.symmetry = %d' % len(self.symmetry)

	def trim_header(self,lines):
		lines.remove(lines[0])

	def generate_plain_contents(self,lines):
		contents = ""
		for line in lines:
			cleared_line = clear_whitespaces(line)	
			contents += cleared_line
		return contents	

	def evaluate_assignment(self,assignment):
		name = re.findall('^[a-zA-Z]+\w*',assignment)[0]
		expr = 'self.' + name + '=[]'
		exec(expr)
		values = re.findall('([\d,\.]+)[;\]]',assignment)
		for val in values:
			expr = 'self.' + name + '.append([' + val + '])'
			exec(expr)



importer = matlab_importer()
#importer.import_file('geometry_test.m')
importer.import_file('brick_new2.m')
export_flagshyp3d(importer.nodes,importer.elements,importer.boundary,'ansys')
export_xml(importer.nodes,importer.elements,importer.boundary,'from_matlab.xml')
export_sexp(importer.nodes,importer.elements,importer.boundary,'from_matlab.sexp')

