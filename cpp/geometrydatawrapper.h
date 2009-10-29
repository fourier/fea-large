#ifndef __GEOMETRYDATAWRAPPER_H__
#define __GEOMETRYDATAWRAPPER_H__

// precompiled headers
#include "std.h"

// local includes
#include "dataloaderwrapper.h"
#include "elements.h"
#include "GeometryData.h"

template<typename ElementT/* = axisymmetric::triangle::element6*/>
class geometrydata_wrapper : 
	public data_loader_wrapper<LGeometryData,ElementT>
{
public:	
	// parent type definition
	typedef data_loader_wrapper<LGeometryData,ElementT> Parent;
	// type definitions from parent
	typedef typename Parent::DataLoader DataLoader;
	typedef typename Parent::Element Element;
	typedef typename Parent::ElementsArray ElementsArray;
	typedef typename Parent::DataPointer DataPointer;
	typedef typename Parent::Node Node;
	typedef typename Parent::BoundaryItem BoundaryItem;
	typedef typename Parent::PrescribedBoundaryArray PrescribedBoundaryArray;

	// nodes array definitions
	typedef LGeometryData::NodesVector NodesVector;
	typedef LGeometryData::NodesConstIterator NodesConstIterator;
	typedef LGeometryData::NodesIterator NodesIterator;
public:
	geometrydata_wrapper() :  Parent() {}
	virtual void load_elements(ElementsArray& array);
	virtual void prepare_index(index_matrix_interface& m);

protected:
	virtual void prepare_boundary();
	virtual void prepare_nodes();
protected:
	using Parent::data_;
	using Parent::ptr_array_;
	using Parent::prescribed_boundary_array_;
	using Parent::nodes_vector_;
};

namespace
{
template<typename ElementT>
class TriangleToElement
{
public:
	TriangleToElement( const LGeometryData::NodesVector& nodes ) : nodes_(nodes) 
	{
	}
	
	ElementT operator ( ) ( const LGeometryData::LTriangle& tri ) const 
	{
		ElementT el;
		// TODO: check carefully with orientation
		el.node(0).dof[0] = nodes_[tri.ind1].x;el.node(0).dof[1] = nodes_[tri.ind1].y;
		el.node(1).dof[0] = nodes_[tri.ind2].x;el.node(1).dof[1] = nodes_[tri.ind2].y;
		el.node(2).dof[0] = nodes_[tri.ind3].x;el.node(2).dof[1] = nodes_[tri.ind3].y;
		el.node(3).dof[0] = nodes_[tri.ind12].x;el.node(3).dof[1] = nodes_[tri.ind12].y;
		el.node(4).dof[0] = nodes_[tri.ind23].x;el.node(4).dof[1] = nodes_[tri.ind23].y;
		el.node(5).dof[0] = nodes_[tri.ind31].x;el.node(5).dof[1] = nodes_[tri.ind31].y;
		el.calculate_volume();
		return el;
	}
private:
	const LGeometryData::NodesVector& nodes_;
};

template<typename ElementT>
class LNodeToNode
{
public:
	typename ElementT::NodeT operator( ) ( const LGeometryData::LNode& lnode ) const
	{
		typename ElementT::NodeT node;
		node.dof[0] = lnode.x;
		node.dof[1] = lnode.y;
		return node;
	}
};

} // anonymous namespace


template<typename ElementT>
void geometrydata_wrapper<ElementT>::load_elements(typename geometrydata_wrapper<ElementT>::ElementsArray& array)
{
	const typename LGeometryData::NodesVector& nodes = data_->GetNodesArray();
	const typename LGeometryData::TrianglesVector& triangles = data_->GetTrianglesArray();

	array.clear();
	array.resize(triangles.size());
	std::transform(triangles.begin(),triangles.end(),array.begin(),TriangleToElement<ElementT>(nodes));
	ptr_array_ = &array;
}

template<typename ElementT>
void geometrydata_wrapper<ElementT>::prepare_boundary()
{
	const typename LGeometryData::BoundaryVector& boundary = data_->GetBoundaryArray();
	const typename DataLoader::SymmetricVector& sym = data_->GetSymmetricArray();
	
	prescribed_boundary_array_.clear();
	prescribed_boundary_array_.resize(boundary.size()+ sym.size());

	typename LGeometryData::BoundaryConstIterator it = boundary.begin();
	size_type i = 0;
	for ( ; it != boundary.end(); ++ it, ++ i )
	{
		switch (it->doftype)
		{
		case LGeometryData::LBoundaryNode::EIgnore:
			prescribed_boundary_array_[i] = BoundaryItem(it->index,Parent::EFree);
			break;
		case LGeometryData::LBoundaryNode::EFixedX:
			{
				Node node;
				node.dof[0] = it->Z1;
				prescribed_boundary_array_[i] = BoundaryItem(it->index,Parent::EPrescribedX,node);
			}
			break;
		case LGeometryData::LBoundaryNode::EFixedY:
			{
				Node node;
				node.dof[1] = it->Z2;
				prescribed_boundary_array_[i] = BoundaryItem(it->index,Parent::EPrescribedY,node);
			}
			break;
		case LGeometryData::LBoundaryNode::EFixedBoth:
			{
				Node node;
				node.dof[0] = it->Z1;
				node.dof[1] = it->Z2;
				prescribed_boundary_array_[i] = BoundaryItem(it->index,Parent::EPrescribedXY,node);
			}
			break;
		default:
			assert(false);
		};
	}
	Node empty_node;
	typename DataLoader::SymmetricConstIterator it1 = sym.begin();
	for ( ; it1 != sym.end(); ++ it1, ++ i )
		prescribed_boundary_array_[i] = BoundaryItem(*it1,Parent::EPrescribedX,empty_node);
}

template<typename ElementT>
void geometrydata_wrapper<ElementT>::prepare_index(index_matrix_interface& m)
{
	const LGeometryData::TrianglesVector& triangles = data_->GetTrianglesArray();
	LGeometryData::TrianglesConstIterator it = triangles.begin();
	size_type i = 0;
	for ( ; it != triangles.end(); ++ it, ++ i  )
	{
		m(i,0) = it->ind1*2;
		m(i,1) = it->ind1*2+1;

		m(i,2) = it->ind2*2;
		m(i,3) = it->ind2*2+1;

		m(i,4) = it->ind3*2;
		m(i,5) = it->ind3*2+1;

		m(i,6) = it->ind12*2;
		m(i,7) = it->ind12*2+1;

		m(i,8) = it->ind23*2;
		m(i,9) = it->ind23*2+1;

		m(i,10) = it->ind31*2;
		m(i,11) = it->ind31*2+1;
	}
}

template<typename ElementT>
void geometrydata_wrapper<ElementT>::prepare_nodes()
{
	nodes_vector_.clear();
	const LGeometryData::NodesVector& nodes = data_->GetNodesArray();
	nodes_vector_.resize(nodes.size());
	std::transform(nodes.begin(),nodes.end(),nodes_vector_.begin(),LNodeToNode<ElementT>());
}




#endif // __GEOMETRYDATAWRAPPER_H__
