#!/usr/bin/python

from base import Point
from triangles import Triangle3,Triangle6

A1 = Triangle3( Point(2,2), Point(1,2), Point(1,1) )
g = [[1,1,2]]
A1.set_values_in_nodes(g)
print A1.g(0,Point(1.25,1.25))
print "6 nodes element"
A2 = Triangle6( Point(2,2), Point(1,2), Point(1,1))
print len(A2.p)
print A2.p[0].x, A2.p[0].y
print A2.p[1].x, A2.p[1].y
print A2.p[2].x, A2.p[2].y
print A2.p[3].x, A2.p[3].y
print A2.p[4].x, A2.p[4].y
print A2.p[5].x, A2.p[5].y
print "Approximations of x+y:"
def g(P = Point()):
	return P.x+P.y
val = [[g(A1.p[i]) for i in range(len(A1.p))]]
print val
A1.set_values_in_nodes(val)
print A1.g(0,Point(1.2,.25))

val = [[g(A2.p[i]) for i in range(len(A2.p))]]
print val
A2.set_values_in_nodes(val)
print A1.g(0,Point(1.2,.25))
print "Dump triangles:"
print A1
