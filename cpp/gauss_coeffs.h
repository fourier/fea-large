#ifndef __GAUSS_COEFFS_H__
#define __GAUSS_COEFFS_H__

// precompiled headers
#include "std.h"

// Template class describing Gauss nodes and weights for every element type
// GaussN - number of Gauss integration nodes
// NodesN - number of nodes in element (3 for triangle element)
template<size_type NodesN, size_type GaussN>
struct gauss_coeffs
{
	enum 
	{
		GaussNumber = GaussN,
		NodesNumber = NodesN
	};

	static const value_type coeffs[GaussN][NodesN];
	static const value_type weights[GaussN];
};

// Gauss coefficients collectios for triangle elements
namespace triangle
{
	// Linear order of integration:
	// 1 Gauss intergration point - center of the triangle 
	typedef gauss_coeffs<3,1> gauss_tri1p;
	// Quadratic order of integration:
	// 3 Gauss integration points - centers of every side of the triangle
	typedef gauss_coeffs<3,3> gauss_tri3p;
	// Cubic order of integration:
	// 4 Gauss integration points - center of the triangle, 3 nodes - corners of "small" triangle
	typedef gauss_coeffs<3,4> gauss_tri4p;
	// Quintic order of integration:
	// 7 Gauss integration points - center of the triangle, 3 nodes - corners of "small"
	// internal triangle, 3 nodes - centers of sides of "small" internal triangle
	typedef gauss_coeffs<3,7> gauss_tri7p;
}

#endif // __GAUSS_COEFFS_H__
