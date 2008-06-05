#ifdef TRIANGULATION
#undef VOID
#define VOID void
extern "C"
{
	#include <triangle.h>
}
#endif

#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cassert>
#include <fstream>
#include <string>
#include <strstream>
#include <sstream>
#include <cassert>
#include <stdexcept>
#include <cmath>

#include "Geometrydata.h"

const double LGeometryData::s_delta_x = 0.4;
const double LGeometryData::s_delta_y = 0.4;

#ifdef USE_STD_EXCEPTIONS
#define THROW_RETURN(arg) throw(std::logic_error(arg));
#else
#define THROW_RETURN(arg) return false;
#endif

namespace 
{
	typedef LGeometryData::LNode Node;
	bool isNodesEqual(const Node& node1, const Node& node2)
	{
		return ( ( fabs( node1.x - node2.x ) < LGeometryData::s_delta_x ) &&
			( fabs( node1.y - node2.y ) < LGeometryData::s_delta_y ) );
	};

	struct IsNodesEqual : 
		public std::binary_function<Node,Node,bool>
	{
		inline 
		result_type operator()( const first_argument_type& node1,
								const second_argument_type& node2 ) const
		{
			return isNodesEqual(node1,node2);
		}
	};
	
};

namespace 
{
	typedef LGeometryData::LSegment Segment;
	typedef LGeometryData::LNode Node;
	typedef LGeometryData::LTriangle Triangle;
	bool isSegmentsEqual(const Segment& seg1,const Segment& seg2)
	{
		return ( ( seg1.n1 == seg2.n1 && seg1.n2 == seg2.n2 ) ||
			 ( seg1.n1 == seg2.n2 && seg1.n2 == seg2.n1 ) );
	};

	struct IsSegmentsEqual :
		public std::binary_function<Segment,Segment,bool>
	{
		
		inline 
		result_type operator() ( const first_argument_type& seg1,
								 const second_argument_type& seg2 ) const
		{
			return isSegmentsEqual(seg1,seg2);
		}
	};

	bool ifSegmentHasIndex(const Segment& seg, int index)
	{
		return ( seg.n1 == index || seg.n2 == index );
	}

	struct IfSegmentHasIndex :
		public std::binary_function<Segment,int,bool> 
	{
		inline
		result_type operator() ( const first_argument_type& seg, 
								 second_argument_type index ) const
		{
			return ifSegmentHasIndex( seg, index );
		}
	};

	bool ifTriangleHasIndex(const Triangle& tr, int index)
	{
		return ( tr.ind1 == index || tr.ind2 == index || tr.ind3 == index ||
				 tr.ind12 == index || tr.ind23 == index || tr.ind31 == index );
	};

	struct IfTriangleHasIndex :
		public std::binary_function<Triangle,int,bool>
	{
		inline
		result_type operator() ( const first_argument_type& tr, 
								 second_argument_type index ) const
		{
			return ifTriangleHasIndex( tr, index );
		};
	};
};

namespace 
{
	typedef LGeometryData::LBoundaryNode BNode;
	
	bool ifBoundaryNodeHasIndex(const BNode& nod, int index)
	{
		return nod.index == index;
	};

	struct IfBoundaryNodeHasIndex :
		public std::binary_function<BNode,int,bool>
	{
		inline
		result_type operator() ( const first_argument_type& nod, 
								 second_argument_type index ) const
		{
			return ifBoundaryNodeHasIndex( nod, index );
		};
	};
		
};

namespace 
{
	typedef LGeometryData::NodesVector NodesVector;
	typedef LGeometryData::SegmentsVector SegmentsVector;
	typedef LGeometryData::TrianglesVector TrianglesVector;
	typedef LGeometryData::BoundaryVector BoundaryVector;
	typedef LGeometryData::SymmetricVector SymmetricVector;
	bool isNodeOneLessThanTwoByX(const NodesVector& v, int ind1, int ind2)
	{
		return v[ind1].x < v[ind2].x;
	};

	struct IsNodeOneLessThanTwoByX : 
		public std::binary_function<int,int,bool>
	{
		IsNodeOneLessThanTwoByX(const NodesVector& v):v_(v){};
		
		result_type operator() ( first_argument_type ind1, second_argument_type ind2)
		{
			return isNodeOneLessThanTwoByX(v_,ind1,ind2);
		};

		const NodesVector& v_;
	};

	struct SegmentsDecreaseIndexesMoreThan : 
		public std::binary_function<int,int,void>
	{
	};
	struct TrianglesDecreaseIndexesMoreThan
	{
	};
	struct BoundaryDecreaseIndexesMoreThan
	{
	};
	struct SymmentricDecreaseIndexesMoreThan
	{
	};
};

std::ostream& operator << ( std::ostream& s, const LGeometryData& data)
{
	data.DropData(s);
	return s;
}

std::istream& operator >> ( std::istream& s, LGeometryData& data)
{
	data.ReadData(s);
	//	data.ReadFile2Data(s);
	return s;
}

LGeometryData::LGeometryData(void)
{
	m_vecNodesArray.reserve(1000);
	m_vecSymmetricArray.reserve(1000);
	m_eBoundaryCondition = ENoConditions;
}

LGeometryData::~LGeometryData(void)
{
}

bool LGeometryData::AddNode(double x, double y)
{
	LNode node(x,y);
	NodesConstIterator it = 
		std::find_if(m_vecNodesArray.begin(),
					 m_vecNodesArray.end(),
					 std::binder2nd<IsNodesEqual>(IsNodesEqual(),node));

	if ( it == m_vecNodesArray.end() )
		m_vecNodesArray.push_back(node);
	else 
		return false;
	return true;
}

bool LGeometryData::AddNode(const LNode& node)
{
	NodesConstIterator it = 
		std::find_if(m_vecNodesArray.begin(),
					 m_vecNodesArray.end(),
					 std::binder2nd<IsNodesEqual>(IsNodesEqual(),node));
	

	if ( it == m_vecNodesArray.end() )
		m_vecNodesArray.push_back(node);
	else 
		return false;
	return true;
}

bool LGeometryData::AddSegment(double x1,double y1,double x2,double y2)
{
	
	LNode nod1(x1,y1);
	LNode nod2(x2,y2);
	if ( isNodesEqual(nod1,nod2) )
		return false;

	int n1 = -1;
	int n2 = -1;
	
	for ( unsigned int i = 0; i < m_vecNodesArray.size(); i++ )
	{
		if ( isNodesEqual(LNode(x1,y1),m_vecNodesArray[i]) )
			n1 = i;
		if ( isNodesEqual(LNode(x2,y2),m_vecNodesArray[i]) )
			n2 = i;
	};

	if ( n1 == -1 || n2 == -1 ) 
		return false;

	LSegment seg(n1,n2);

	SegmentsConstIterator it = std::find_if(m_vecSegmentsArray.begin(),
		m_vecSegmentsArray.end(),
		std::binder2nd<IsSegmentsEqual>(IsSegmentsEqual(),seg));
	if ( it == m_vecSegmentsArray.end() )
		m_vecSegmentsArray.push_back(seg);
	else
		return false;
	
	return true;
}



bool LGeometryData::AddSegment(const LNode& nod1,const LNode& nod2)
{
	if ( isNodesEqual(nod1,nod2) )
		return false;

	int n1 = -1;
	int n2 = -1;
	
	for ( unsigned int i = 0; i < m_vecNodesArray.size(); i++ )
	{
		if ( isNodesEqual(nod1,m_vecNodesArray[i]) )
			n1 = i;
		if ( isNodesEqual(nod2,m_vecNodesArray[i]) )
			n2 = i;
	};

	if ( n1 == -1 || n2 == -1 ) 
		return false;

	LSegment seg(n1,n2);

	SegmentsConstIterator it = std::find_if(m_vecSegmentsArray.begin(),
		m_vecSegmentsArray.end(),
		std::binder2nd<IsSegmentsEqual>(IsSegmentsEqual(),seg));
	if ( it == m_vecSegmentsArray.end() )
		m_vecSegmentsArray.push_back(seg);
	else
		return false;
	
	return true;
}


LGeometryData::NodesConstIterator 
	LGeometryData::GetCorrespondingNode(double x,double y)
{
	LNode nod(x,y);
	return std::find_if(m_vecNodesArray.begin(),
		m_vecNodesArray.end(),
		std::binder2nd<IsNodesEqual>(IsNodesEqual(),nod));
}

LGeometryData::NodesConstIterator 
	LGeometryData::GetCorrespondingNode(const LNode& nod)
{
	return std::find_if(m_vecNodesArray.begin(),
		m_vecNodesArray.end(),
		std::binder2nd<IsNodesEqual>(IsNodesEqual(),nod));
}

#ifdef TRIANGULATION
bool LGeometryData::TriangulateArrays(double angle, double area, EDegree degree, NodesVector& nodes,
	SegmentsVector& segments, TrianglesVector& triangles,NodesVector& holes)
{
	if ( !nodes.size() || !segments.size() )
		return false;

	std::stringstream triswitches;
	
	triswitches << "pXQzq" << angle << "a" << area;
	if ( degree == ESixNodes )
		triswitches << "o2";	

	struct triangulateio in;
	struct triangulateio out;                                           
// Basic input initialization
	in.numberofpoints = (int)nodes.size();
	in.pointlist = new REAL[in.numberofpoints*2];
	in.pointattributelist = NULL;
	in.pointmarkerlist = new int [nodes.size()]; 
	in.numberofsegments = (int)segments.size();
	in.segmentlist = new int[in.numberofsegments*2];
	in.segmentmarkerlist = NULL;
    in.holelist = holes.size() ? new REAL[holes.size()*2] : NULL;
	in.numberofholes = int(holes.size());
	in.numberofpointattributes = 0;
	in.regionlist = NULL;
	in.numberofregions = 0;
	

// Basic output initialization	
	out.pointlist = (REAL *) NULL;            
	out.pointattributelist = (REAL *) NULL;
	out.pointmarkerlist = (int *) NULL; 
	out.trianglelist = (int *) NULL;
	out.triangleattributelist = (REAL *) NULL;
	out.neighborlist = (int *) NULL;
	out.segmentlist = (int *) NULL;
	out.segmentmarkerlist = (int *) NULL;
	out.edgelist = (int *) NULL;
	out.edgemarkerlist = (int *) NULL;
// Let's fill input data

	int i = 0;
	NodesConstIterator it1 = nodes.begin();
	
	for ( ; it1 != nodes.end(); it1++ )
	{
		in.pointlist[i++] = it1->x;
		in.pointlist[i++] = it1->y;
	}

	for ( i = 0; i < in.numberofpoints; i++ )
		in.pointmarkerlist[i] = nodes[i].bnd;

	SegmentsConstIterator it2 = segments.begin();
	
	for ( i = 0; it2 != segments.end(); it2++ )
	{
		in.segmentlist[i++] = it2->n1;
		in.segmentlist[i++] = it2->n2;
	}

	if ( holes.size() )
	{
		for ( it1 = holes.begin(),i = 0; it1 != holes.end(); it1++ )
		{
			in.holelist[i++] = it1->x;
			in.holelist[i++] = it1->y;
		}

	}

//	DropData2PolyFile("data.poly",nodes,segments);
	triangulate(const_cast<char*>(triswitches.str().c_str()), &in, &out, NULL);

	nodes.clear();
	segments.clear();
	triangles.clear();

	assert( nodes.size() == 0 );
	assert( segments.size() == 0 );
	assert( triangles.size() == 0 );

	for ( i = 0; i < out.numberofpoints*2; i+=2 )
	{
		LNode node(out.pointlist[i],
			out.pointlist[i+1],
			EBoundaryType(out.pointmarkerlist[i/2]));
		nodes.push_back(node);
	}

	for ( i = 0; i < out.numberofsegments*2; i+=2 )
	{
		LSegment seg(out.segmentlist[i],out.segmentlist[i+1]);
		segments.push_back(seg);
	}

	int nodesPerTriangle = degree == ESixNodes ? 6 : 3; 

	for ( i = 0; 
		  i < out.numberoftriangles*nodesPerTriangle; 
		  i += nodesPerTriangle )
	{
		LTriangle tng;
		tng.ind1 = out.trianglelist[i];
		tng.ind2 = out.trianglelist[i+1];
		tng.ind3 = out.trianglelist[i+2];
		if ( degree == ESixNodes )
		{
			tng.ind23 = out.trianglelist[i+3];
			tng.ind31 = out.trianglelist[i+4];
			tng.ind12 = out.trianglelist[i+5];
		}

		triangles.push_back(tng);
	};

	// Let's free occupied memory


	// in:
	delete [] in.pointlist;
	delete [] in.pointmarkerlist;
	delete [] in.segmentlist;
	if ( in.holelist )
		delete [] in.holelist;
	// out:
	trifree(out.pointlist);
	trifree(out.edgelist);
	trifree(out.segmentlist);
	trifree(out.pointmarkerlist);
	trifree(out.segmentmarkerlist);
	trifree(out.trianglelist);

	// all ok, return true
	return true;
}

bool LGeometryData::Triangulate( double angle, double area, EDegree degree )
{
	if ( !m_vecNodesArray.size() || !m_vecSegmentsArray.size() )
		return false;
	bool bResult = TriangulateArrays(angle, area, degree, m_vecNodesArray,
		m_vecSegmentsArray,m_vecTrianglesArray,m_vecHolesArray);

	SetBoundaryCondition(ENoConditions);

	return bResult;
}
#endif // TRIANGULATION

#ifdef _DEBUG
void LGeometryData::DropData2PolyFile( const std::string& fname,
	NodesVector& nodes,SegmentsVector& segments) 
{
	if ( !nodes.size() || !segments.size() )
		return;
	
	char* buf = new char [1024*1024*16];

	std::strstream stream(buf,1024*1024*16);
		    
	std::ofstream file;

	file.open( fname.c_str(), std::ios_base::out | std::ios_base::trunc );

	stream << "# A dump of the data from LGeometry class of Geometry2d project" << std::endl;

	stream<<(unsigned int)nodes.size()<<" "<<2<<" "<<0<<" "<<1<<std::endl;
	unsigned int i = 0;
	
	for ( ; i < nodes.size(); i++ )
	{
		stream<<i+1<<" "<<nodes[i].x<<" "<<nodes[i].y<<" "<<nodes[i].bnd<<std::endl;

	}
	stream << "# Now let's describe segments:" << std::endl;
	
	stream << (unsigned int)segments.size() << " " << 0 << std::endl;
	
	i = 0;
	for ( SegmentsConstIterator it = segments.begin(); it != segments.end(); it++ )
	{
		stream<<i<<" "<<it->n1+1<<" "<<it->n2+1<<std::endl;
		i++;	
	}
	stream << "# Number of holes:" << std::endl;
	stream << 0 << std::endl << std::ends;
	file << stream.str();
	
	file.close();

	delete [] buf;
}
#endif

bool LGeometryData::AddHole(double x, double y)
{
	LNode node;
	node.x = x;
	node.y = y;
	NodesConstIterator it = 
		std::find_if(m_vecHolesArray.begin(),
					 m_vecHolesArray.end(),
					 std::binder2nd<IsNodesEqual>(IsNodesEqual(),node));
	

	if ( it == m_vecHolesArray.end() )
		m_vecHolesArray.push_back(node);
	else 
		return false;
	return true;
}

bool LGeometryData::AddHole(const LNode& node)
{
	NodesConstIterator it = 
		std::find_if(m_vecHolesArray.begin(),
					 m_vecHolesArray.end(),
					 std::binder2nd<IsNodesEqual>(IsNodesEqual(),node));
	

	if ( it == m_vecHolesArray.end() )
		m_vecHolesArray.push_back(node);
	else 
		return false;
	return true;
}

bool LGeometryData::AddHole(NodesConstIterator it)
{
	LNode node;
	node.x = it->x;
	node.y = it->y;
	NodesConstIterator it1 = 
		std::find_if(m_vecHolesArray.begin(),
					 m_vecHolesArray.end(),
					 std::binder2nd<IsNodesEqual>(IsNodesEqual(),node));
	

	if ( it1 == m_vecHolesArray.end() )
		m_vecHolesArray.push_back(node);
	else 
		return false;
	return true;
}

void LGeometryData::SetBoundary(NodesConstIterator it, EBoundaryType bnd)
{
	NodesVector::iterator it1 = std::find_if(m_vecNodesArray.begin(),
					 m_vecNodesArray.end(),
					 std::binder2nd<IsNodesEqual>(IsNodesEqual(),*it));

	it1->bnd = bnd;

}

bool LGeometryData::DropData( std::ostream& stream ) const
{
	using std::endl;
	
	stream << "# A dump of the data from LGeometry class of Geometry2d project" << endl;
	stream << "# the '#' character means beginning of the comment, empty lines will be skipped" << endl;
	stream << "# Format: " << endl;
	stream << "# 1-st line - global geomtry settings" << endl;
	stream << "# 2-nd line - boundary settings" << endl;
	stream << "# Data sections: nodes, segments, holes,triangles,[symmetry][boundary]" << endl;
	stream << "# First meaningfull line should contain all global data parameters" << endl;
	stream << "# First line format: " << endl;
	stream << "# <Number of nodes> <number of segments> <number of holes> <number of triangles>" << endl << endl;

	typedef unsigned int uint;

	stream << (uint)m_vecNodesArray.size() << " " << (uint)m_vecSegmentsArray.size() << " ";
	stream << (uint)m_vecHolesArray.size() << " " << (uint)m_vecTrianglesArray.size() << endl << endl;

	stream << "# The second meaningfull line contains information about boundary conditions" << endl;
	stream << "# It has following format:" << endl;
	stream << "# <Boundary condition type> <# of nodes in symmetric condition> <# of bounary nodes>" << endl;
	stream << "# <Boundary condition type> is a string. It may has one of the following values:" << endl;
	stream << "# \"NoConditions\" - then no boundary conditions at all" << endl;
	stream << "# \"Stresses\" - then the stress boundary condition only" << endl;
	stream << "# \"Displacements\" - then the displacements boundary condition only" << endl;
	stream << "# \"SymAndStresses\" - then the symmetry and stress boundary condition given" << endl;
	stream << "# \"SymAndDisplacements\" - than the symmetry and displacements boundary condition given" << endl << endl;

	std::string condition;

	switch ( m_eBoundaryCondition )
	{
	case ENoConditions:
		condition = "NoConditions";
		break;
	case EStress:
		condition = "Stresses";
		break;
	case EDisplacement:
		condition = "Displacements";
		break;
	case ESymAndDisplacement:
		condition = "SymAndDisplacements";
		break;
	case ESymAndStress:
		condition = "SymAndStresses";
		break;
	default:
		throw(false);
		break;
	};
		
	stream << condition << " " << (uint)m_vecSymmetricArray.size() << " ";
	stream << (uint)m_vecBoundaryArray.size() << " " << endl << endl;

	stream << "# The following lines will contail information about nodes" << endl;
	stream << "# The format of node line is: " << endl;
	stream << "# <Zero-based index> <X> <Y> <boundary type>" << endl;
	stream << "# there X and Y - floating-point coordinates of the node" << endl;
	stream << "# and <boundary type> is \"0\" for inner and \"1\" for outer nodes" << endl << endl;

	uint i = 0;

	for ( ; i < (uint)m_vecNodesArray.size(); i++ )
	{
		stream << i << " " << m_vecNodesArray[i].x << " " << m_vecNodesArray[i].y;
		stream << " " << m_vecNodesArray[i].bnd << endl;
	}
	
	stream << endl;
	stream << "# The next couple of lines describes a segment structures" << endl;
	stream << "# Format: " << endl;
	stream << "# <First node number> <Second node number>" << endl << endl;
	
	SegmentsConstIterator it = m_vecSegmentsArray.begin();

	for ( ; it != m_vecSegmentsArray.end(); it ++ )
	{
		stream << it->n1 << " " << it->n2 << endl;
	};

	stream << endl;

	stream << "# The next contains information about holes" << endl;
	stream << "# Format: " << endl;
	stream << "# <zero-based hole index> <X> <Y>" << endl;
	stream << "# The holes format is similar to nodes format except of boundary markers" << endl << endl;

	for ( i = 0; i < m_vecHolesArray.size(); i++ )
	{
		stream << i << " " << m_vecHolesArray[i].x << " " << m_vecHolesArray[i].y << endl;
	}
	stream << endl;

	stream << "# Triangles section" << endl;
	stream << "# Format: " << endl;
	stream << "# <index1> <index2> <index3> [<index12> <index23> <index31>]" << endl;
	stream << "# where index1[2,3] is a index of corners of triangle in counterclockwise form" << endl;
	stream << "# In case of 6-nodes triangulation index12[23,31] is an indexes" << endl;
	stream << "# of a points in the middle of corresponding edge." << endl;
	stream << "# For example, index12 is an index of the point that lies on a edge" << endl;
	stream << "# witch vertex indexes is index1 and index2" << endl << endl;

	TrianglesConstIterator itr = m_vecTrianglesArray.begin(); 
	for ( ; itr != m_vecTrianglesArray.end(); itr++ )
	{
		stream << itr->ind1 << " " << itr->ind2 << " " << itr->ind3;
		if ( itr->ind12 != -1 &&
			 itr->ind23 != -1 &&
			 itr->ind31 != -1 )
			 stream << " " << itr->ind12 << " " << itr->ind23 << " " << itr->ind31;
		stream << endl;
	}

	stream << endl;

	if ( m_eBoundaryCondition == ESymAndDisplacement ||
		m_eBoundaryCondition == ESymAndStress && m_vecSymmetricArray.size())
	{
		stream << "# The next line is a very-long :) line of indexes of nodes," << endl;
		stream << "# that lies on a symmetry boundary" << endl <<endl;
		for ( i = 0; i < m_vecSymmetricArray.size() ; i++ )
			stream << m_vecSymmetricArray[i] << " ";
		stream << endl << endl;
	}

	if ( m_eBoundaryCondition != ENoConditions && m_vecBoundaryArray.size() )
	{
		stream << "# The following lines contains indexes of boundary nodes" << endl;
		stream << "# and boundary condition params" << endl;
		stream << "# Format: " << endl;
		if ( m_eBoundaryCondition == EStress || 
			m_eBoundaryCondition == ESymAndStress )
		{
			stream << "# <index> <S1> <S2> <dof type>" << endl;
			stream << "# Index - node index, S1 and S2 - stresses" << endl;
			stream << "# dof type - 1 - use 1st, 2 - use 2nd, 3 - use both" << endl;
		}
		else
		{
			stream << "# <index> <u1> <u2> <dof type>" << endl;
			stream << "# Index - node index, u1 and u2 - displacements" << endl;
			stream << "# dof type - 1 - use 1st, 2 - use 2nd, 3 - use both" << endl;
		}
		stream << endl;

		BoundaryConstIterator iter = m_vecBoundaryArray.begin(); 
		
		for ( ; iter != m_vecBoundaryArray.end(); iter++ )
		{
			stream << iter->index << " " << iter->Z1 << " " << iter->Z2 << " " << iter->doftype << endl;
		};
		stream << endl;
	}

	stream << "# End of data" ;

	return true;
}

bool LGeometryData::ReadFile2Data(std::istream& stream)
{
	int vertexcount = 0, trianglescount = 0;
	int number = 0;
	double x = 0, y = 0, min = 0;
	int v1 = 0, v2 = 0, v3 = 0, v4 = 0, v5 = 0, v6 = 0;
	m_vecBoundaryArray.clear();
	m_vecHolesArray.clear();
	m_vecNodesArray.clear();
	m_vecSegmentsArray.clear();
	m_vecSymmetricArray.clear();
	m_vecTrianglesArray.clear();
	m_eBoundaryCondition = ENoConditions;

	stream >> vertexcount >> trianglescount;
	for ( int i = 0; i < vertexcount; i++ )
	{
		stream >> number >> x >> y ;
		m_vecNodesArray.push_back(LNode(x,y));
	}

	for ( int i = 0; i < trianglescount; i++ )
	{
		stream >> v1 >> v2 >> v3 >> v4 >> v5 >> v6;
		m_vecTrianglesArray.push_back(LTriangle(v1,v2,v3,v4,v5,v6));
	}




    return true;
}


bool LGeometryData::RemoveAllContainingNode( NodesConstIterator it )
{
	// let's find index of iterator in vector
	int index = 0;
	for ( ; index < (int)m_vecNodesArray.size(); index++ )
		if ( it->x == m_vecNodesArray[index].x &&
			 it->y == m_vecNodesArray[index].y )
			break;

	// remove segments
	m_vecSegmentsArray.remove_if(std::binder2nd<IfSegmentHasIndex>(
		IfSegmentHasIndex(), index ) );
	// remove triangles
	m_vecTrianglesArray.remove_if(std::binder2nd<IfTriangleHasIndex>(
		IfTriangleHasIndex(), index ) );
	
	SetBoundaryCondition(ENoConditions);

	return true;
}

void LGeometryData::SetCoordinates(NodesConstIterator it, double x, double y)
{
	NodesVector::iterator it1 = std::find_if(m_vecNodesArray.begin(),
					 m_vecNodesArray.end(),
					 std::binder2nd<IsNodesEqual>(IsNodesEqual(),*it));

	if ( it1 != m_vecNodesArray.end() )
	{
		it1->x = x;
		it1->y = y;
	};
}

void LGeometryData::SetBoundaryCondition(EBoundaryCondition cond)
{
	m_eBoundaryCondition = cond;

	if ( m_eBoundaryCondition == ENoConditions )
	{
		m_vecBoundaryArray.clear();
		m_vecSymmetricArray.clear();
	};
}

bool LGeometryData::ChangeCondition(NodesConstIterator it, double Z1, double Z2,int doffixed)
{
	if ( m_eBoundaryCondition == ENoConditions )
		return false;
	int index = 0;
	for ( ; index < (int)m_vecNodesArray.size(); index++ )
		if ( isNodesEqual(*it,m_vecNodesArray[index]) )
			break;
	;
	if ( index == m_vecNodesArray.size() )
		return false;

	BoundaryIterator it1 = find_if( m_vecBoundaryArray.begin(),
		m_vecBoundaryArray.end(),
		std::binder2nd<IfBoundaryNodeHasIndex>(IfBoundaryNodeHasIndex(),
			index));

	if ( it1 != m_vecBoundaryArray.end() )
	{
		it1->Z1 = Z1;
		it1->Z2 = Z2;
		it1->doftype = doffixed;

	}
	else
	{
		LBoundaryNode node(index,Z1,Z2,doffixed);
		m_vecBoundaryArray.push_back(node);
	}
	return true;
}

LGeometryData::SegmentsConstIterator LGeometryData::GetCorrespondingSegment(
		NodesConstIterator it1, NodesConstIterator it2) const
{
	int index1 = 0;
	int index2 = 0;

	for ( ; index1 < (int)m_vecNodesArray.size(); index1++ )
		if ( isNodesEqual(*it1,m_vecNodesArray[index1]) )
			break;

	for ( ; index2 < (int)m_vecNodesArray.size(); index2++ )
		if ( isNodesEqual(*it2,m_vecNodesArray[index2]) )
			break;

	if ( index1 == m_vecNodesArray.size() )
		return m_vecSegmentsArray.end();
	if ( index2 == m_vecNodesArray.size() )
		return m_vecSegmentsArray.end();

	SegmentsConstIterator pos = m_vecSegmentsArray.begin();
	SegmentsConstIterator end = m_vecSegmentsArray.end();
		
	while ( pos != m_vecSegmentsArray.end() )
	{
		pos = find_if(pos,end,
				std::binder2nd<IfSegmentHasIndex>(IfSegmentHasIndex(),index1));
		if ( pos != m_vecSegmentsArray.end() )
			if ( ifSegmentHasIndex(*pos,index2) )
			{
				 break;
			}
			else
				pos++;
	} 
	
	return pos;
}

bool LGeometryData::IsNodesLiesOnOneLine( int index_begin, int index_end, int index_node)
{
	// let's check array boundaries

	using std::min;
	using std::max;

	if ( index_begin < 0 || index_end < 0 || index_node < 0 )
		return false;

	if ( index_begin >= (int)m_vecNodesArray.size() ||
		 index_end   >= (int)m_vecNodesArray.size() ||
		 index_node  >= (int)m_vecNodesArray.size() )
		 return false;

	const double precision = min(s_delta_x,s_delta_y)/100;

	if ( m_vecNodesArray[index_node].x >= 
			min(m_vecNodesArray[index_begin].x,m_vecNodesArray[index_end].x) && 
		 m_vecNodesArray[index_node].x <= 
			max(m_vecNodesArray[index_begin].x,m_vecNodesArray[index_end].x) &&
		m_vecNodesArray[index_node].y >= 
			min(m_vecNodesArray[index_begin].y,m_vecNodesArray[index_end].y) && 
		 m_vecNodesArray[index_node].y <= 
			max(m_vecNodesArray[index_begin].y,m_vecNodesArray[index_end].y) )
	{
		if ( fabs(Det(m_vecNodesArray[index_begin],m_vecNodesArray[index_end],
			     m_vecNodesArray[index_node])) < precision )
				 return true;
		else return false;
	}
	else 
		return false;

	return true;
}

double LGeometryData::Det(const LNode& begin, const LNode& end, const LNode& node)
{
	return (begin.x*node.y-begin.y*node.x)+(node.x*end.y-end.x*node.y)-
		(begin.x*end.y-begin.y*end.x);
}

bool LGeometryData::SetSymmetricArray(NodesConstIterator it1, NodesConstIterator it2 )
{
	if ( isNodesEqual(*it1, *it2 ) )
		return false;

	int index1 = 0;
	int index2 = 0;

	for ( ; index1 < (int)m_vecNodesArray.size(); index1++ )
		if ( isNodesEqual(*it1,m_vecNodesArray[index1]) )
			break;

	for ( ; index2 < (int)m_vecNodesArray.size(); index2++ )
		if ( isNodesEqual(*it2,m_vecNodesArray[index2]) )
			break;

	if ( index1 == m_vecNodesArray.size() )
		return false;
	if ( index2 == m_vecNodesArray.size() )
		return false;

	// now it's an corresponding indexes index1 and index2

	m_vecSymmetricArray.clear();
	
	m_vecSymmetricArray.push_back(index1);

	for ( int i = 0; i < (int)m_vecNodesArray.size(); i++ )
	{
		if ( i != index1 && i != index2 )
			if ( IsNodesLiesOnOneLine(index1,index2,i) )
				m_vecSymmetricArray.push_back(i);
	}
	
	m_vecSymmetricArray.push_back(index2);
	
	IsNodeOneLessThanTwoByX my_less( m_vecNodesArray );

	sort(m_vecSymmetricArray.begin(),
		 m_vecSymmetricArray.end(),
		 my_less);


	return true;
}

bool LGeometryData::AddToBoundaryArray( NodesConstIterator it1, 
				NodesConstIterator it2, double Z1,double Z2,int doffixed)
{
	if ( isNodesEqual(*it1, *it2 ) )
		return false;

	int index1 = 0;
	int index2 = 0;

	for ( ; index1 < (int)m_vecNodesArray.size(); index1++ )
		if ( isNodesEqual(*it1,m_vecNodesArray[index1]) )
			break;

	for ( ; index2 < (int)m_vecNodesArray.size(); index2++ )
		if ( isNodesEqual(*it2,m_vecNodesArray[index2]) )
			break;

	if ( index1 == m_vecNodesArray.size() )
		return false;
	if ( index2 == m_vecNodesArray.size() )
		return false;

	// now it's an corresponding indexes index1 and index2

	BoundaryConstIterator it;
	bool bFound = false;
	

	for ( it = m_vecBoundaryArray.begin();
		it != m_vecBoundaryArray.end(); it++)
		if ( it->index == index1 )
		{
			bFound = true;
			break;
		}

	if ( !bFound )
		m_vecBoundaryArray.push_back(LBoundaryNode(index1,Z1,Z2,doffixed));
	
	bFound = false;

	for ( it = m_vecBoundaryArray.begin();
		it != m_vecBoundaryArray.end(); it++)
		if ( it->index == index2 )
		{
			bFound = true;
			break;
		}

	if ( !bFound )
		m_vecBoundaryArray.push_back(LBoundaryNode(index2,Z1,Z2,doffixed));

	int siz = int(m_vecNodesArray.size());

	for ( int i = 0; i < siz; i++ )
	{
		if ( i != index1 && i != index2 )
			if ( IsNodesLiesOnOneLine(index1,index2,i) )
			{
				bFound = false;
	
				for (it = m_vecBoundaryArray.begin();
					 it != m_vecBoundaryArray.end(); it++)
					if ( it->index == i )
					{
						bFound = true;
						break;
					}

				if ( !bFound )
					m_vecBoundaryArray.push_back(LBoundaryNode(i,Z1,Z2,doffixed));
			}
	}
	
	return true;
}

bool LGeometryData::ReadData(std::istream& stream) throw (std::logic_error)
{
 	std::string str;

	enum EReadState
	{
		EFirstConfigLine,
		ESecondConfigLine,
		ENodesList,
		ESegmentsList,
		EHolesList,
		ETrianglesList,
		ESymmetryList,
		EBoundaryList
	};

	EReadState CurState = EFirstConfigLine;
	EReadState NextState = EFirstConfigLine;

	int nds = -1;
	int segms = -1;
	int trngls = -1;
	int hls = -1;

	EBoundaryCondition cond = EUknownCondition;
	int syms = -1;
	int bnds = -1;

	bool bComment = false;
    bool bEmpty = false;

	bool bThreeNodePerTriangle = true;

	int nCount = 0;

	// empty all arrays
	m_vecBoundaryArray.clear();
	m_vecHolesArray.clear();
	m_vecNodesArray.clear();
	m_vecSegmentsArray.clear();
	m_vecSymmetricArray.clear();
	m_vecTrianglesArray.clear();

	while ( std::getline(stream, str) )
	{
		nCount++;

		const char * ptr = str.c_str();
		if ( *ptr == '\0' )
			continue;

		bComment = bEmpty = false;
		do 
		{
			if ( !isspace(*ptr) )
			{
				if ( *ptr == '#' )
				{
					bComment = true;
					break;
				}
				else 
				{
					bComment = bEmpty = false;
					break;
				}
			}
			bEmpty = true;
		} while ( *++ptr );

		if ( bComment || bEmpty ) 
			continue;

		std::istringstream istr(str);

		bool bContinue = false;
		do 
		{
			bContinue = false;

			if ( CurState != NextState )
			{
				CurState = NextState;			
			}	
		
			switch ( CurState )
			{
			case EFirstConfigLine:
				{
					istr >> nds >> segms >> hls >> trngls;
					if ( istr.bad() || istr.fail() )
						THROW_RETURN("Error while reading first config line: wrong syntax");
					if ( nds < 0 || segms < 0 || hls < 0 || trngls < 0 ) 
						THROW_RETURN("Error while reading first config line: less than zero values");
					m_vecNodesArray.reserve(nds);

					NextState = ESecondConfigLine;	
				}
				continue;
			case ESecondConfigLine:
				{
					std::string strCond;
					istr >> strCond >> syms >> bnds ;
					if ( istr.bad() || istr.fail() )
						THROW_RETURN("Error while reading second config line: wrong syntax");

					if ( strCond == "NoConditions" )
						m_eBoundaryCondition = ENoConditions;
					else if ( strCond == "Stresses" )
						m_eBoundaryCondition = EStress;
					else if ( strCond == "Displacements" )
						m_eBoundaryCondition = EDisplacement;
					else if ( strCond == "SymAndStresses" )
						m_eBoundaryCondition = ESymAndStress;
					else if ( strCond == "SymAndDisplacements" )
						m_eBoundaryCondition = ESymAndDisplacement;
					else 
						THROW_RETURN("Error while reading second config line: invalid boundary type specifier");

					if ( syms < 0 && bnds < 0 )
						THROW_RETURN("Error while reading second config line: less than zero values");

					if ( m_eBoundaryCondition != ENoConditions )
					{
						if ( syms )
						{
							if ( m_eBoundaryCondition == ESymAndDisplacement ||
								m_eBoundaryCondition == ESymAndStress ) 
								m_vecSymmetricArray.reserve(syms);
							else
								THROW_RETURN("Error while reading second config line: \
											symmetry nodes given but not symmetry condition");
						}
					} else if ( syms || bnds )
						THROW_RETURN("Error while reading second config line: \
											No conditions but nodes given");

					NextState = ENodesList;
				}
				continue;
			case ENodesList:
				{
					if ( nds == 0 )
						return true;

					if ( nds == m_vecNodesArray.size() )
					{
						NextState = ESegmentsList;
						continue;
					}
					else 
						NextState = ENodesList;


					LNode node;
					int index;

					istr >> index >> node.x >> node.y >> *(int*)(&node.bnd);
					if ( istr.bad() || istr.fail() )
						THROW_RETURN("Error while reading nodes line: \
											wrong syntax");
					switch ( node.bnd )
					{
					case EInner:
					case EOuter:
					case ENotMarked:
						break;
					default:
						THROW_RETURN("Error while reading nodes line: \
											wrong boundary type");
					}
					
					m_vecNodesArray.push_back(node);

					NextState = nds == m_vecNodesArray.size() ? ESegmentsList:ENodesList;
				}
				continue;
			case ESegmentsList:
				{
					if ( !segms )
					{
						bContinue = true;
						NextState = EHolesList;
						continue;
					}
					if ( segms == m_vecSegmentsArray.size() )
					{
						NextState = EHolesList;
						continue;
					}
					else 
						NextState = ESegmentsList;
					
					LSegment seg;
				
					istr >> seg.n1 >> seg.n2;
					if ( istr.bad() || istr.fail() )
						THROW_RETURN("Error in segment line: \
									wrong syntax");

					if ( seg.n1 < 0 || seg.n2 < 0 )
						THROW_RETURN("Error in segment line: \
									less than zero values");

					if ( seg.n1 >= (int)m_vecNodesArray.size() ||
						seg.n2 >= (int)m_vecNodesArray.size() )
						return false;
					
					m_vecSegmentsArray.push_back(seg);
					
					NextState = segms == m_vecSegmentsArray.size() ? EHolesList : ESegmentsList;
				}
				continue;
			case EHolesList:
				{
					if ( !hls )
					{
						bContinue = true;
						NextState = ETrianglesList;
						continue;
					}				
					if ( hls == m_vecHolesArray.size() )
					{
						NextState = ETrianglesList;
						continue;
					}
					else 
						NextState = EHolesList;

					LNode node;
					int index; 
					
					istr >> index >> node.x >> node.y;
					if ( istr.bad() || istr.fail() )
						THROW_RETURN("Error in holes line: \
									wrong syntax");

					m_vecHolesArray.push_back(node);

					NextState = hls == m_vecHolesArray.size() ? ETrianglesList:EHolesList;
				}
				continue;
			case ETrianglesList:
				{
					if ( !trngls )
					{
						bContinue = true;
						NextState = ESymmetryList;
						continue;
					}
					if ( trngls == m_vecTrianglesArray.size() )
					{
						NextState = ESymmetryList;
						continue;
					}
					else
						NextState = ETrianglesList;
					
					LTriangle tng;

					std::vector<int> vec;
					vec.reserve(6);
					int index;
					
					while ( istr >> index )
						vec.push_back(index);

					if ( vec.size() != 3 && vec.size() != 6 )
						THROW_RETURN("Error in reading triangles line: \
									number of indexes is not 3 or 6");

					if ( m_vecTrianglesArray.size() == 0 )
						bThreeNodePerTriangle = ( vec.size() == 3 );

					if ( bThreeNodePerTriangle != ( vec.size() == 3 ) )
						THROW_RETURN("Error in reading triangles line: \
									different nodes count in different lines");

					tng.ind1 = vec[0];
					tng.ind2 = vec[1];
					tng.ind3 = vec[2];

					if ( !bThreeNodePerTriangle )
					{
						tng.ind12 = vec[3];
						tng.ind23 = vec[4];
						tng.ind31 = vec[5];
					}
		
					if ( tng.ind1 < 0 || tng.ind1 >= (int)m_vecNodesArray.size() )
						THROW_RETURN("Error in reading triangles line: \
									wrong index for 1st node");
					if ( tng.ind2 < 0 || tng.ind2 >= (int)m_vecNodesArray.size() )
						THROW_RETURN("Error in reading triangles line: \
									wrong index for 2nd node");
					if ( tng.ind3 < 0 || tng.ind3 >= (int)m_vecNodesArray.size() )
						THROW_RETURN("Error in reading triangles line: \
									wrong index for 3rd node");
					if ( !bThreeNodePerTriangle )
					{
						if ( tng.ind12 < 0 || tng.ind12 >= (int)m_vecNodesArray.size() )
							THROW_RETURN("Error in reading triangles line: \
									wrong index for 12 node");
						if ( tng.ind23 < 0 || tng.ind23 >= (int)m_vecNodesArray.size() )
							THROW_RETURN("Error in reading triangles line: \
									wrong index for 23 node");
						if ( tng.ind31 < 0 || tng.ind31 >= (int)m_vecNodesArray.size() )
							THROW_RETURN("Error in reading triangles line: \
									wrong index for 31 node");
					}

					m_vecTrianglesArray.push_back(tng);
				
					if ( trngls == m_vecTrianglesArray.size() )
					{
						NextState = ESymmetryList;
						continue;
					}
					else
						NextState = ETrianglesList;

					NextState = trngls == m_vecTrianglesArray.size() ? ESymmetryList : ETrianglesList;
				}
				continue;
			case ESymmetryList:
				{
					if ( syms == 0 )
					{
						NextState = EBoundaryList;
						bContinue = true;
						continue;
					}

					int index;

					for ( int i = 0; i < syms; i++ )
					{
						istr >> index;
						if ( istr.bad() || istr.fail() )
							THROW_RETURN("Error in reading symmetry line: \
										wrong syntax");
						if ( index < 0 || index >= (int)m_vecNodesArray.size() )
							THROW_RETURN("Error in reading symmetry line: \
										wrong index");
						m_vecSymmetricArray.push_back(index);
					}

					if ( syms == m_vecSymmetricArray.size() )
						NextState = EBoundaryList;
					else
						THROW_RETURN("Error in reading symmetry line: \
									number of readed indexes is not equal to specification");
				}
				continue;
			case EBoundaryList:
				{
					if ( bnds == m_vecBoundaryArray.size() )
						return true;
					else
						NextState = EBoundaryList;
					
					LBoundaryNode node;
					
					istr >> node.index >> node.Z1 >> node.Z2;
					if ( istr.bad() || istr.fail() )
						THROW_RETURN("Error in reading boundary line: \
									wrong syntax");
					if ( node.index < 0 || node.index >= (int)m_vecNodesArray.size() )
						THROW_RETURN("Error in reading boundary line: \
									wrong index");
					
					int doftype = 3;
					istr >> doftype;
					node.doftype = doftype;
					m_vecBoundaryArray.push_back(node);
			
				}
				continue;
			default:
				return false;

			};
					
		} while ( bContinue );
	};


	return true;
}


LGeometryData& LGeometryData:: operator = ( const LGeometryData& data )
{
	// empty all arrays
	m_vecBoundaryArray.clear();
	m_vecHolesArray.clear();
	m_vecNodesArray.clear();
	m_vecSegmentsArray.clear();
	m_vecSymmetricArray.clear();
	m_vecTrianglesArray.clear();

	m_eBoundaryCondition = data.m_eBoundaryCondition;

	std::copy(data.m_vecBoundaryArray.begin(),data.m_vecBoundaryArray.end(),
		std::back_inserter(m_vecBoundaryArray));
	std::copy(data.m_vecHolesArray.begin(),data.m_vecHolesArray.end(),
		std::back_inserter(m_vecHolesArray));
	std::copy(data.m_vecNodesArray.begin(),data.m_vecNodesArray.end(),
		std::back_inserter(m_vecNodesArray));
	std::copy(data.m_vecSegmentsArray.begin(),data.m_vecSegmentsArray.end(),
		std::back_inserter(m_vecSegmentsArray));
	std::copy(data.m_vecSymmetricArray.begin(),data.m_vecSymmetricArray.end(),
		std::back_inserter(m_vecSymmetricArray));
	std::copy(data.m_vecTrianglesArray.begin(),data.m_vecTrianglesArray.end(),
		std::back_inserter(m_vecTrianglesArray));

	return *this;
}

