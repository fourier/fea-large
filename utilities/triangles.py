#/usr/bin/python

from base import Point, factorial, Element


# Class which represents triangle finite elements
# Triangle consist of 3 Points
# All 3 points should be enumerated counter-clockwise
# functions to be overloaded in derived classes, for example
# to use 6-nodes elements, are reset_tail(where nodes on edges 
# should be calculated) and N(i) - shape functions
# It is general requirement to use only 3 points as input.
# All other nodes should be calculated in derived classes,
# because the order of other nodes is important in shape functions
class Triangle3(Element):
	"Class of plane triangles FE with 3 nodes"
	def __init__(self,p1 = Point(),p2 = Point(),p3 = Point() ):
		self.p = []
		self.values = []
		self.deformed = []
		self.p.append(p1)
		self.p.append(p2)
		self.p.append(p3)
		self.deformed.append(p1)
		self.deformed.append(p2)
		self.deformed.append(p3)
	def reset_tail(self):
		tail_start = len(self.p) - 3
		if tail_start > 0:
			del self.p[3:len(self.p)]
	def reset(self,p1 = Point(),p2 = Point(), p3 = Point() ):
		self.p[0] = p1 
		self.p[1] = p2
		self.p[2] = p3
		self.deformed[0] = p1
		self.deformed[1] = p2
		self.deformed[2] = p3
		self.reset_tail()
	def from_L(self,L):
		result = Point()
		for i in range(len(L)):
			result.x += self.p[i].x*L[i]
			result.y += self.p[i].y*L[i]
		return x
	def a_coef(self,i):
		a = 0
		if i == 0:
			a = self.deformed[1].x*self.deformed[2].y-self.deformed[2].x*self.deformed[1].y
		elif i == 1:
			a = self.deformed[2].x*self.deformed[0].y-self.deformed[0].x*self.deformed[2].y
		elif i == 2:
			a = self.deformed[0].x*self.deformed[1].y-self.deformed[1].x*self.deformed[0].y
		else:
			print "Triangle3::a_coef:wrong argument i =", i
			return
		return a
	def b_coef(self,i):
		b = 0
		if i == 0:
			b = self.deformed[1].y - self.deformed[2].y
		elif i == 1:
			b = self.deformed[2].y - self.deformed[0].y
		elif i == 2:
			b = self.deformed[0].y - self.deformed[1].y
		else:
			print "Triangle3::b_coef:wrong argument i =", i
			return
		return b
	def c_coef(self,i):
		c = 0
		if i == 0:
			c = self.deformed[2].x - self.deformed[1].x
		elif i == 1:
			c = self.deformed[0].x - self.deformed[2].x
		elif i == 2:
			c = self.deformed[1].x - self.deformed[0].x
		else:
			print "Triangle3::c_coef:wrong argument i =", i
			return	
		return c
	def L( self,i,point = Point() ):
		if i >= len(self.p):
			print "Triangle3::L:wrong argument i =",i
		l = (self.a_coef(i)+self.b_coef(i)*point.x+self.c_coef(i)*point.y)/\
			(2.0*self.S_triangle())
		return l
	def N( self,i,point = Point() ):
		if i >= len(self.p):
			print "Triangle3::N:wrong argument i =",i
			return
		return self.L(i,point)
	def S_triangle(self):
		result = (self.p[1].x*self.p[2].y-self.p[2].x*self.p[1].y)+\
			(self.p[2].x*self.p[0].y-self.p[0].x*self.p[2].y)+\
			(self.p[0].x*self.p[1].y-self.p[1].x*self.p[0].y)
		return result/2.0
	def Integral_L_coords(self,m,n,k):
		numerator=factorial(m)*factorial(n)*factorial(k)*\
			2*self.S_triangle()
		denominator=factorial(m+n+k+2)
		return numerator/denominator



# Class for triangle finite elements
# with 6 nodes
# nodes 4,5,6 situated on edges of triangle
# with counter-clockwise order
# 4 in edge (1,2), 5 in edge (2,3) and 6 in (3,1)
class Triangle6(Triangle3):
	def __init__(self,p1 = Point(),p2 = Point(),p3 = Point() ):
		Triangle3.__init__(self,p1,p2,p3)
		self.reset_tail()
	def reset_tail(self):
		Triangle3.reset_tail(self)
		p4 = Point()
		p5 = Point()
		p6 = Point()
		p4.x = self.p[0].x+(self.p[1].x-self.p[0].x)/2.0
		p4.y = self.p[0].y+(self.p[1].y-self.p[0].y)/2.0
		p5.x = self.p[1].x+(self.p[2].x-self.p[1].x)/2.0
		p5.y = self.p[1].y+(self.p[2].y-self.p[1].y)/2.0
		p6.x = self.p[2].x+(self.p[0].x-self.p[2].x)/2.0
		p6.y = self.p[2].y+(self.p[0].y-self.p[2].y)/2.0
		self.p.append(p4)
		self.p.append(p5)
		self.p.append(p6)
	def N(self,i,point = Point()):
		if i >= len(self.p):
			print "Triangle6::N:Wrong argument i =", i
			return
		res = 0
		if i == 0 or i == 1 or i == 2:
			res = (2*self.L(i)-1)*self.L(i)
		elif i == 3:
			res = 4*self.L(0)*self.L(1)
		elif i == 4:
			res = 4*self.L(1)*self.L(2)
		else: # i == 5  
			res = 4*self.L(2)*self.L(0)
		return res


