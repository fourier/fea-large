#!/usr/bin/python
from base import Point
# Program for generating 3d brick mesh
# by given geometry parameters
#    
#   |      *-----*----*
# W |   D /     /P   /|
#   |    /     /    / |   
# A |   *-----*----* ||
#   |   |          |M|*
# L | H |-----N----| /
#   |   |          |/
# L |   *----------*
#   |"""     L 
#     S 
#
#    | z 
#    |   / y
#    | /
#    /--------x 
#
#
#     M = 3 
#   /   |   \
#
# *-----------*
# |   |   |   |
# |   |   |   |  \
# |---|---|---|   |
# |   |   |   |   -
# |   |   |   |   N = 4  
# |---|---|---|   -
# |   |   |   |  /|
# |   |   |   |   |
# |---|---|---|  /
# |   |   |   |
# |   |   |   |
# *-----------*
#
# Geometry description:
# S - distance from brick to the Ox,Oy,Oz
# L - length of the brick
# D - deep of the brick
# H - height of the brick
# N,M,P - numbers of parts split by elements
# u(S)=0 - this wall of brick is fixed
# u(S+L)=Ux - prescribed displacement  

def export3dbrick(filename,S,L,H,D,N,M,P,Ux):
	# open a file
	#fil = open(filename,"w+")
	elems = N*M*P # number of elements
	nods = (N+1)*(M+1)*(P+1) # number of nodes
	nodes = []
	# algorithm of generating nodes is based on
	# addition node by node in all horisontal layers
	# total number of horisontal layers is P+1
	#
	# Lets count steps - and therefore distances btw nodes
	step_x = 1.*L/P
	step_y = 1.*D/M
	step_z = 1.*H/N
	# loops btw layers and fill nodes array
	for k in range(P+1):
		x = S + step_x*k 
		#print "layer %d" % k
		for i in range(N+1):
			z = S + step_z*i
			for j in range(M+1):
				y = S + step_y*j
				#print " node(%d,%d)" % (i,j)
				p = Point(x,y,z)
				nodes.append(p)
	# loops btw layers and fill elements array
	for k in range(P):
		for i in range(N):
			for j in range(M):
				print "element %d %d %d" % (k,i,j)	
	# close a file
	#fil.close()
	 


export3dbrick("simple_brick.m",1,6,3,3,2,2,2,0.05)	 
