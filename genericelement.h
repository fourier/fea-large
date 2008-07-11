#ifndef __GENERICELEMENT_H__
#define __GENERICELEMENT_H__

// precompiled headers
#include "std.h"

// local includes
#include "functions.h"
#include "gauss_coeffs.h"

//
// Nodes 
// 
// Parameters: number degrees on freedom :DoF
// type of stored values: T
//
//
template<size_type DoF,typename T = value_type>
struct Node
{
	typedef T value_type;
	Node() { memset(dof,0,sizeof(dof)); }
	value_type dof[DoF];
};

//
// serialization of nodes from input stream
// 
template<size_type DoF>
std::istream& operator >> (std::istream& s,Node<DoF>& n)
{
	for (size_type i = 0; i < DoF; ++ i)
	{
        s >> n.dof[i];
		s.get();
	}
	return s;
}


class non_copyconstructable
{
public:
	non_copyconstructable()
	{
	};
private:
	non_copyconstructable(const non_copyconstructable&)
	{
	};
};

// NodesN - number of nodes per element
// DoF - number of degrees of freedom per node
// VoigtSize - size of voigt vector of strain:
// can be 3 (for plane elements), 4(for axisymmetric elements)
// 6 for 3d elements
template<size_type NodesN,size_type DoF,size_type VoigtSize, typename GaussNodes = triangle::gauss_tri7p>
class element
{
public:
	enum
	{
		DofNumber = DoF,
		NodesNumber = NodesN,
		VoigtNumber = VoigtSize
	};
public:
	typedef element<NodesN,DoF,VoigtSize,GaussNodes> SelfT;
	typedef Node<DoF> NodeT;
	typedef typename NodeT::value_type value_type;
	typedef std::vector<NodeT> NodesArrayT;
	typedef GaussNodes GaussCoeffs;
	typedef boost::multi_array<value_type, 2> GaussForms;
	typedef boost::multi_array<value_type, 3> GaussFormsDerivatives;
public:
	// C'tor/D'tor 
	element() : volume_(0) 
	{ 
		presize_arrays();
	}
	
	// unsafe construction
	virtual ~element(){};
	// Getters/Setters
	NodeT& node(size_type i) { return nodes_[i]; }
	const NodeT& node(size_type i) const { return nodes_[i]; }
	const value_type volume() const { return volume_;}

	// form function
	// @param index - index of form function[0-(number_of_nodes-1)]
	// @param node - point to calculate form function
	virtual value_type form(size_type index, const NodeT& node) const = 0;

	// function for calculation derivatives of form functions
	// @param index - index of form function
	// @param dof - index of degree of freedom to integrate by
	virtual value_type dform(size_type index, size_type dof, const NodeT& node) const = 0;

	NodeT center() const
	{
		NodeT node;
		for ( size_type i = 0; i < NodesNumber; ++ i )
			for ( size_type j = 0; j < DoF; ++ j )
				node.dof[j] += nodes_[i].dof[j]/(value_type)NodesNumber; 
		return node;
	}

	// function for calculation of the gauss point
	// index - index of the Gauss point
	template<typename G>
	NodeT gauss_point(size_type index) const 
	{
// TODO: check why it can't be compiled
//		BOOST_STATIC_ASSERT(GaussCoeffs::NodesNumber == NodesNumber);
		assert(index >= 0 && index < G::GaussNumber);

		NodeT node;
		for ( size_type i = 0; i < G::NodesNumber; ++ i )
			for ( size_type j = 0; j < DoF; ++ j )
				node.dof[j] += nodes_[i].dof[j]*G::coeffs[index][i];
		
		return node;
	}

	const element& operator = (const element& el)
	{
		copy(el);
		return *this;
	}

	virtual void calculate_volume() = 0;
	virtual void calculate()
	{
		calculate_volume();
		init_gauss_nodes();
	}
protected:
	//
	void presize_arrays()
	{
		nodes_.reserve(NodesNumber);
		gauss_nodes_.reserve(GaussCoeffs::GaussNumber);

		boost::array<typename GaussForms::index,2> dims1 = 
			{{GaussNodes::GaussNumber,NodesNumber}};
		boost::array<typename GaussFormsDerivatives::index,3> dims2 = 
			{{GaussNodes::GaussNumber, NodesNumber, DofNumber}};       

		gauss_forms_.resize(dims1);
		gauss_derivatives_.resize(dims2);
	}
	void init_gauss_nodes()
	{
		gauss_nodes_.clear();
		for ( size_type i = 0; i < GaussCoeffs::GaussNumber; ++ i )
			gauss_nodes_.push_back(gauss_point<GaussCoeffs>(i));

		for ( size_type i = 0; i < GaussCoeffs::GaussNumber; ++ i )
		{
			for ( size_type j = 0; j < NodesNumber; ++ j )
			{
				gauss_forms_[i][j] = form(j,gauss_nodes_[i]);
				for ( size_type k = 0; k < DofNumber; ++ k )
					gauss_derivatives_[i][j][k] = dform(j,k,gauss_nodes_[i]);
			}
		}
	}
	void copy(const element& el)
	{
		for ( size_type i = 0; i < NodesNumber; ++ i )
			nodes_[i] = el.nodes_[i];
		calculate();
	}
protected:
	// nodes of FE
	NodesArrayT nodes_;
	// volume
	value_type volume_;
	// array of gauss nodes
	NodesArrayT gauss_nodes_;
	// values of shape functions in gauss nodes
	// access: gauss_shapes_[i][j]
	// where i - number of gauss node, j - number of form function
	GaussForms gauss_forms_;
	// derivatives of shape functions in gauss nodes
	// access: gauss_derivatives_[i][j][k] 
	// where i - number of gauss node,
	// j - number of form function
	// k - number of DoF
	// Example: DN1/dy in gauss point 2 of element
	// assuming what numeration from 0
	// gauss_derivatives_[1][0][1];
	// since y is a second DoF
	GaussFormsDerivatives gauss_derivatives_;
};


#endif // __GENERICELEMENT_H__
