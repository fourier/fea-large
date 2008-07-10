// precompiled headers
#include "std.h"
// Local includes
#include "elements.h"
#include "gauss_coeffs.h"

/*
// test for element_plane_triangle_p6
void test_element_axisymmetric_triangle_p6()
{
	typedef axisymmetric::triangle::element6 testelement;
	typedef testelement::NodeT NodeT;
	NodeT nodes[6];
	
	// test geometry
	nodes[0].dof[0] = 1; nodes[0].dof[1] = 1;
	nodes[1].dof[0] = 4; nodes[1].dof[1] = 1;
	nodes[2].dof[0] = 1; nodes[2].dof[1] = 4;
	nodes[3].dof[0] = 2.5; nodes[3].dof[1] = 1;
	nodes[4].dof[0] = 2.5; nodes[4].dof[1] = 2.5;
	nodes[5].dof[0] = 1; nodes[5].dof[1] = 2.5;

	// test element
	testelement el(nodes);
	// test copyconstructor 
	testelement el1(el);
	
	// first of all - test local coordinates
	for ( size_type i = 0; i < 3; ++ i )
		for ( size_type j = 0; j < 3; ++ j )
			assert( el.local(i,nodes[j]) == kroneker_delta(i,j) );

	// second, test form functions
	for ( size_type i = 0; i < 6; ++ i )
	{
		for ( size_type j = 0; j < 6; ++ j )
		{
			assert( el.form(i,nodes[j]) == kroneker_delta(i,j) );
		}
	}

	// third, test derivatives of form functions
	// TODO: implement
	for ( size_type i = 0; i < 6; ++ i )
	{
		for ( size_type j = 0; j < 6; ++ j )
		{
			for ( size_type k = 0; k < 3; ++ k );
		}
	}


	// 4rth, test integration
	// For triangle elements analytical formulae:
	//                                        a!b!c!
	// Integral(L1^a * L2^b * L3^c dx dy) = ----------*2S
	//                                      (a+b+c+2)!
	// where S - volume of the triangle

	// Prerequisites: calculate integral of the simple function L1^a*L2^b*L3^c
	value_type S = el.volume();
	int a = 1, b = 1, c = 1;
	value_type integral = 2*S*factorial(a)*factorial(b)*factorial(c)/factorial(a+b+c+2);

	value_type value = 0;
	std::cout << "test 1 node integration:" << std::endl;
	for ( size_type i = 0; i < triangle::gauss_tri1p::GaussNumber; ++ i)
	{
		NodeT node = el.gauss_point<triangle::gauss_tri1p>(i);
		value_type fun = pow(el.local(0,node),a)*pow(el.local(1,node),b)*pow(el.local(2,node),c);
		value += S*fun*triangle::gauss_tri1p::weights[i];
	}
	std::cout << integral << " " << value << " " << fabs(integral-value) << std::endl;

	std::cout << "test 3 node integration:" << std::endl;
	value = 0;
	for ( size_type i = 0; i < triangle::gauss_tri3p::GaussNumber; ++ i)
	{
		NodeT node = el.gauss_point<triangle::gauss_tri3p>(i);
		value_type fun = pow(el.local(0,node),a)*pow(el.local(1,node),b)*pow(el.local(2,node),c);
		value += S*fun*triangle::gauss_tri3p::weights[i];
	}
	std::cout << integral << " " << value << " " << fabs(integral-value) << std::endl;

	std::cout << "test 4 node integration:" << std::endl;
	value = 0;
	for ( size_type i = 0; i < triangle::gauss_tri4p::GaussNumber; ++ i)
	{
		NodeT node = el.gauss_point<triangle::gauss_tri4p>(i);
		value_type fun = pow(el.local(0,node),a)*pow(el.local(1,node),b)*pow(el.local(2,node),c);
		value += S*fun*triangle::gauss_tri4p::weights[i];
	}
	std::cout << integral << " " << value << " " << fabs(integral-value) << std::endl;

	std::cout << "test 7 node integration:" << std::endl;
	value = 0;
	for ( size_type i = 0; i < triangle::gauss_tri7p::GaussNumber; ++ i)
	{
		NodeT node = el.gauss_point<triangle::gauss_tri7p>(i);
		value_type fun = pow(el.local(0,node),a)*pow(el.local(1,node),b)*pow(el.local(2,node),c);
		value += S*fun*triangle::gauss_tri7p::weights[i];
	}
	std::cout << integral << " " << value << " " << fabs(integral-value) << std::endl;


	std::cout <<"=======================================" << std::endl;
	for ( size_type i = 0; i < el.NodesNumber; ++ i)
		std::cout << el.node(i).dof[0] << " " << el.node(i).dof[1] << std::endl;
	std::cout <<"=======================================" << std::endl;
	value = 0;
	for ( size_type i = 0; i < triangle::gauss_tri3p::GaussNumber; ++ i)
	{
		NodeT node = el.gauss_point<triangle::gauss_tri3p>(i);
		std::cout << node.dof[0] << " " << node.dof[1] << " : weight = " << triangle::gauss_tri3p::weights[i];
		std::cout << " local1 = " << el.local(0,node) << " local2 = " << el.local(1,node) << " local3 = " << el.local(2,node);
		std::cout << std::endl;
	}

	assert(fabs(value-integral)<=1); // simple check

}
*/

void test_element_plane_triangle_p6()
{
	typedef plane::triangle::element6 testelement;
	typedef testelement::NodeT NodeT;
	NodeT nodes[6];
	
	// test geometry
	nodes[0].dof[0] = 1; nodes[0].dof[1] = 1;
	nodes[1].dof[0] = 4; nodes[1].dof[1] = 1;
	nodes[2].dof[0] = 1; nodes[2].dof[1] = 4;
	nodes[3].dof[0] = 2.5; nodes[3].dof[1] = 1;
	nodes[4].dof[0] = 2.5; nodes[4].dof[1] = 2.5;
	nodes[5].dof[0] = 1; nodes[5].dof[1] = 2.5;

	// test element
	testelement el(nodes);

	// test copyconstructor 
	testelement el1(el);
	//el1 = el;
	
	// first of all - test local coordinates
	for ( size_type i = 0; i < 3; ++ i )
		for ( size_type j = 0; j < 3; ++ j )
			assert( el1.local(i,nodes[j]) == kroneker_delta(i,j) );

	// second, test form functions
	for ( size_type i = 0; i < 6; ++ i )
	{
		for ( size_type j = 0; j < 6; ++ j )
		{
			assert( el1.form(i,nodes[j]) == kroneker_delta(i,j) );
		}
	}
	assert(testelement::VoigtNumber == 3);
}

