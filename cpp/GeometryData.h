#pragma once

#include <vector>
#include <list>
#include <map>
#include <string>
#include <ostream>
#include <istream>
#include <stdexcept>

#pragma warning( disable : 4290 )


class LGeometryData
{
public:
	enum EBoundaryType 
	{
		EInner = 0,
		EOuter = 1,
		ENotMarked  = -1
	};

	enum EDegree 
	{
		EThreeNodes,
		ESixNodes,
	};

	enum EBoundaryCondition
	{
		ENoConditions,
		ESymAndStress,
		ESymAndDisplacement,
		EStress,
		EDisplacement,
		EUknownCondition
	};

	struct LNode 
	{
		double x; // X coord
		double y ; // Y coord 
		EBoundaryType bnd; // boundary index
		LNode(double x1=0 ,double y1=0,EBoundaryType bnd1=EInner)
			:x(x1),y(y1),bnd(bnd1){};
	};
	struct LSegment
	{
		int n1; // first point index 
		int n2; // second point index
		EBoundaryType bnd; // boundary index
		LSegment(int nn1=0,int nn2=0,EBoundaryType bnd1=EInner)
			:n1(nn1),n2(nn2),bnd(bnd1){};
	};

	
	struct LTriangle 
	{
		int ind1; // first vertex index
		int ind2; // second vertex index
		int ind3; // third vertex index
		int ind12;
		int ind23;
		int ind31;
		LTriangle(int i1=0,int i2=0,int i3=0,int i4=-1,int i5=-1,int i6=-1)
			: ind1(i1),ind2(i2),ind3(i3),ind12(i4),ind23(i5),ind31(i6)
		{
		};
	};

	struct LBoundaryNode 
	{
		enum EFixedDOF
		{
			EIgnore = 0,
			EFixedX = 1,
			EFixedY = 2,
			EFixedBoth = 3
		};
		int index; // index of boundary node
		double Z1; // acts like S_1 or u_1, depending on bool m_bBoundaryType
		double Z2; // acts like S_2 or u_2 
		int doftype;
		LBoundaryNode(int ind=-1,double z1=0,double z2=0,int dof=EFixedBoth)
			: index(ind),Z1(z1),Z2(z2),doftype(dof) {};
	};

	typedef std::vector<LNode> NodesVector;
	typedef std::list<LSegment> SegmentsVector;
	typedef std::list<LTriangle> TrianglesVector;
	typedef std::list<LBoundaryNode> BoundaryVector;
	typedef std::vector<int> SymmetricVector;
	
	typedef NodesVector::const_iterator NodesConstIterator;
	typedef SegmentsVector::const_iterator SegmentsConstIterator;
	typedef TrianglesVector::const_iterator TrianglesConstIterator;
	typedef BoundaryVector::const_iterator BoundaryConstIterator;
	typedef SymmetricVector::const_iterator SymmetricConstIterator;

	typedef NodesVector::iterator NodesIterator;
	typedef SegmentsVector::iterator SegmentsIterator;
	typedef TrianglesVector::iterator TrianglesIterator;
	typedef BoundaryVector::iterator BoundaryIterator;
	typedef SymmetricVector::iterator SymmetricIterator;

	// C'tor/D'tor
	LGeometryData(void);
	~LGeometryData(void);
	
	// operator == 
	LGeometryData& operator = ( const LGeometryData& );


	// Adders
	bool AddNode(double x, double y);
	bool AddNode(const LNode& node);

	bool AddSegment(double x1,double y1,double x2,double y2);
	bool AddSegment(const LNode& nod1, const LNode& nod2); 

	bool AddHole(double x, double y);
	bool AddHole(const LNode& node);
	bool AddHole(NodesConstIterator it);

	// Getters
	NodesConstIterator GetCorrespondingNode(double x,double y);
	NodesConstIterator GetCorrespondingNode(const LNode& nod);

	const NodesVector& GetNodesArray() const 
	{
		return m_vecNodesArray;
	}

	const NodesVector& GetHolesArray() const 
	{
		return m_vecHolesArray;
	}

	const SegmentsVector& GetSegmentsArray() const
	{
		return m_vecSegmentsArray;
	}

	const TrianglesVector& GetTrianglesArray() const
	{
		return m_vecTrianglesArray;
	}

	const BoundaryVector& GetBoundaryArray() const
	{
		return m_vecBoundaryArray;
	}
	
	const EBoundaryCondition& GetBoundaryCondition(void) const
	{
		return m_eBoundaryCondition;
	};

	const SymmetricVector& GetSymmetricArray(void) const
	{
		return m_vecSymmetricArray;
	};

	// Setters
	
	bool SetSymmetricArray( NodesConstIterator it1, NodesConstIterator it2 );

	bool AddToBoundaryArray( NodesConstIterator it1, NodesConstIterator it2, double Z1,double Z2,
		int doffixed = LBoundaryNode::EFixedBoth);
		
	SegmentsConstIterator GetCorrespondingSegment(NodesConstIterator it1,
		NodesConstIterator it2) const;

	void SetBoundaryCondition(EBoundaryCondition cond);
	
	void SetBoundary(NodesConstIterator it, EBoundaryType bnd);
	
	void SetCoordinates(NodesConstIterator it, double x, double y);

	bool RemoveAllContainingNode( NodesConstIterator it );

	bool ChangeCondition(NodesConstIterator it, double Z1, double Z2,int doffixed);
	
	// Misc 

	bool ReadFile2Data(std::istream& stream);	
	
#ifdef TRIANGULATION
	// Triangulation
	bool Triangulate( double angle, double area, EDegree degree );
#endif

	// Statics
	static const double s_delta_x;
	static const double s_delta_y;

	// Friends

	friend std::ostream& operator << ( std::ostream&, const LGeometryData&);
	friend std::istream& operator >> ( std::istream&, LGeometryData&);
protected:

	bool DropData(std::ostream& stream) const;

	bool ReadData(std::istream& stream) throw (std::logic_error);

#ifdef TRIANGULATION
	static bool TriangulateArrays(double angle, double area, EDegree degree, NodesVector& nodes, 
		SegmentsVector& segments, TrianglesVector& triangles,NodesVector& holes);
#endif

	static double Det(const LNode& begin, const LNode& end, const LNode& node);

	bool IsNodesLiesOnOneLine( int index_begin, int index_end, int index_node);

	#ifdef _DEBUG
	static void DropData2PolyFile( const std::string& fname,
		NodesVector& nodes,SegmentsVector& segments);
	#endif

	// Attributes
	NodesVector m_vecNodesArray;
	NodesVector m_vecHolesArray;
	SegmentsVector m_vecSegmentsArray;	
	TrianglesVector m_vecTrianglesArray;
	BoundaryVector m_vecBoundaryArray;
	SymmetricVector m_vecSymmetricArray;

	EBoundaryCondition m_eBoundaryCondition;

};


std::ostream& operator << ( std::ostream&, const LGeometryData&);
std::istream& operator >> ( std::istream& s, LGeometryData& data);
