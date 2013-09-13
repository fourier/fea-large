#!/usr/bin/python

from sys import argv,exit
from base import Point
from triangles import Triangle3,Triangle6
from exporter import BasePlane,export_matlab
from base import isempty,skipwhitespaces,iscomment,removeempty

# python allows to be referenced only class instances
# arrays and so on transferred by a value to any
# function ( not by reference )
# therefore plane will contain members
# 'symmetry' and 'boundary'
def parsenodefile(filename,plane):
    # open a file and read all lines from it
    file = open(filename,"r")
    lines = file.readlines()
    file.close()
    # clear plane nodes and triangles
    plane.points = []
    plane.triangles_indexed = []
    plane.symmetry = []
    plane.boundary = []
    # prepare for scanning
    nodes_count = 0
    segments_count = 0
    holes_count = 0
    triangles_count = 0
    # current phase of scanning
    scan_phase = 0
    # condition type:
    # "NoConditions" - then no boundary conditions at all
    # "Stresses" - then the stress boundary condition only
    # "Displacements" - then the displacements boundary condition only
    # "SymAndStresses" - then the symmetry and stress boundary 
    # condition given
    # "SymAndDisplacements" - than the symmetry and displacements 
    # boundary condition given
    condition = ''
    nodes_symmetry = 0
    nodes_boundary = 0
    # loop through lines
    for line in lines:
        # skip whitespaces and comments
        if isempty(line) or iscomment(line):
            continue
        if scan_phase == 0: # first line with number of nodes etc.
            values = removeempty(line.split(' '))
            [nodes_count,segments_count,holes_count,triangles_count] = \
                [int(i) for i in values ]
            scan_phase += 1
            continue
        elif scan_phase == 1: # 2nd line with BC information
            values = removeempty(line.split(' '))
            [condition,nodes_symmetry,nodes_boundary] = \
                [i for i in values]

            nodes_symmetry = int(nodes_symmetry)
            nodes_boundary = int(nodes_boundary)
            scan_phase += 1
            continue
        elif scan_phase == 2: # a set of lines with nodes
            values = removeempty(line.split(' '))
            x = float(values[1])
            y = float(values[2])
            plane.points.append(Point(x,y))
            # assume what at least one node in list
            nodes_count -= 1
            if nodes_count == 0:
                scan_phase += 1
            continue
        elif scan_phase == 3: # segments information
            # we don't use segments in this code
            # assume what we have at least 1 segment
            segments_count -= 1
            if segments_count == 0:
                scan_phase += 1
            continue
        elif scan_phase == 4: # triangles section 
            values = removeempty(line.split(' '))
            triangle = [int(i) for i in values]
            plane.triangles_indexed.append(triangle)
            # assume at least one triangle
            triangles_count -= 1
            if triangles_count == 0:
                scan_phase += 1
            continue
        elif scan_phase == 5: # symmetry condition
            values = removeempty(line.split(' '))
            plane.symmetry = [int(i) for i in values]
            scan_phase += 1
            continue
        else : # last phase - prescribed displacements BC
            values = removeempty(line.split(' '))
            boundary = [int(values[0]), # node index \
                        float(values[1]), # ux \
                        float(values[2])] # uy
            # if new format
            if len(values) == 4:
                boundary.append(int(values[3]))
            plane.boundary.append(boundary)
            nodes_boundary -= 1
            if nodes_boundary == 0:
                break
    # function finished


if len(argv) > 1:
    plane = BasePlane()
    values = argv[1].split('.')
    mfname = values[0]
    parsenodefile(argv[1],plane)
    export_matlab(plane,plane.boundary,plane.symmetry,mfname)




