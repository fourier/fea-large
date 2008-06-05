#include "gauss_coeffs.h"
// enumeration of the particular Gauss integration points and weights

// Triangles section
// Taken from Zienkiewicz, vol 1, p222

using triangle::gauss_tri1p;
using triangle::gauss_tri3p;
using triangle::gauss_tri4p;
using triangle::gauss_tri7p;

// 1 integration point
template<>
const value_type gauss_tri1p::coeffs[gauss_tri1p::GaussNumber][gauss_tri1p::NodesNumber] = 
{{1/3.,1/3.,1/3.}};
template<>
const value_type gauss_tri1p::weights[gauss_tri1p::GaussNumber] = {1};

// 3 integration points
template<>
const value_type gauss_tri3p::coeffs[gauss_tri3p::GaussNumber][gauss_tri3p::NodesNumber] = 
{{0.5, 0.5, 0},{0, 0.5, 0.5},{0.5, 0, 0.5}};
template<>
const value_type gauss_tri3p::weights[gauss_tri3p::GaussNumber] = {1./3.,1./3.,1./3.};

// 4 integration points
template<>
const value_type gauss_tri4p::coeffs[gauss_tri4p::GaussNumber][gauss_tri4p::NodesNumber] = 
{{1/3.,1/3.,1/3.},{0.6, 0.2, 0.2}, {0.2, 0.6, 0.2},{0.2, 0.2, 0.6}};
template<>
const value_type gauss_tri4p::weights[gauss_tri4p::GaussNumber] = {-27/48.,25/48.,25/48.,25/48.};


// 7 integration points
const value_type a1 = 0.0597158717;
const value_type b1 = 0.4701420641;
const value_type a2 = 0.7974269853;
const value_type b2 = 0.1012865073;
const value_type w1 = 0.225;
const value_type w2 = 0.1323941527;
const value_type w3 = 0.1259391805;
template<>
const value_type gauss_tri7p::coeffs[gauss_tri7p::GaussNumber][gauss_tri7p::NodesNumber] = 
{{1/3.,1/3.,1/3.},
 {a1,b1,b1},
 {b1,a1,b1},
 {b1,b1,a1}, 
 {a2,b2,b2},
 {b2,a2,b2},
 {b2,b2,a2}};
template<>
const value_type gauss_tri7p::weights[gauss_tri7p::GaussNumber] = 
{w1,w2,w2,w2,w3,w3,w3};
