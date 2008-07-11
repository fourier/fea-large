#ifndef __ELEMENTS_H__
#define __ELEMENTS_H__

#include "genericelement.h"

namespace axisymmetric 
{
namespace triangle 
{
// 6-nodes triangle in axisymmetric case
// wherefore voigt number is 4:
// dim(Srr, Szz,Sff,Srz) = 4
class element6 : public ::element<6,2,4>
{
public:
	typedef ::element<6,2,4> Parent;
public:
	element6() : Parent()
	{
	}
	element6(const element6& el)
	{
		copy(el);
	}

	element6(const NodeT nodes[])
	{
		for ( size_type i = 0; i < NodesNumber; ++ i )
			nodes_[i] = nodes[i];
		calculate();
	}

	// local(L) coordinates
	// @param index - index of local function 0..2
	value_type local(size_type index, const NodeT& node) const
	{
		assert(index >= 0 && index < 3);
		value_type l = 0;
		
		const value_type r1 = nodes_[0].dof[0], z1 = nodes_[0].dof[1];
		const value_type r2 = nodes_[1].dof[0], z2 = nodes_[1].dof[1];
		const value_type r3 = nodes_[2].dof[0], z3 = nodes_[2].dof[1];

		value_type a[3];
		a[0] = r2*z3-r3*z2; // el(2,1)*el(3,2)-el(3,1)*el(2,2);
		a[1] = r3*z1-r1*z3; // el(3,1)*el(1,2)-el(1,1)*el(3,2);
		a[2] = r1*z2-r2*z1; // el(1,1)*el(2,2)-el(2,1)*el(1,2);
		
		value_type b[3];
		
		b[0] = z2-z3; // el(2,2)-el(3,2);
		b[1] = z3-z1; // el(3,2)-el(1,2);
		b[2] = z1-z2; // el(1,2)-el(2,2);

		value_type c[3];
		c[0] = r3-r2; // el(3,1)-el(2,1);
		c[1] = r1-r3; // el(1,1)-el(3,1);
		c[2] = r2-r1; // el(2,1)-el(1,1);

		l = (a[index]+b[index]*node.dof[0]+c[index]*node.dof[1])/(2.*volume_);;

		return l;
	}

	// form function
	// @param index - index of form function[0-(number_of_nodes-1)]
	// @param node - point to calculate form function
	value_type form(size_type index, const NodeT& node) const
	{
		assert(index >= 0 && index < 6);
		value_type f = 0;

		if ( index  < 3 )
			// f = (2*local(el,i,r,z)-1)*local(el,i,r,z);
			f = (2*local(index,node)-1)*local(index,node);
		else
			switch (index)
			{
			case 3:
				// f = 4*local(el,1,r,z)*local(el,2,r,z);
				f = 4*local(0,node)*local(1,node);
				break;
			case 4:
				// f = 4*local(el,2,r,z)*local(el,3,r,z);
				f = 4*local(1,node)*local(2,node);
				break;
			case 5:
				// f = 4*local(el,3,r,z)*local(el,1,r,z);
				f = 4*local(2,node)*local(0,node);
				break;
			}
		return f;
	};
	
	// function for calculation derivatives of form functions
	// @param index - index of form function
	// @param dof - index of degree of freedom to integrate by
	value_type dform(size_type index, size_type dof, const NodeT& node) const
	{
		assert(index >= 0 && index < 6);
		assert(dof >=0 && dof < 3);
		
		const value_type r1 = nodes_[0].dof[0], z1 = nodes_[0].dof[1];
		const value_type r2 = nodes_[1].dof[0], z2 = nodes_[1].dof[1];
		const value_type r3 = nodes_[2].dof[0], z3 = nodes_[2].dof[1];
		
		value_type f = 0;

		switch (dof)
		{
		case 0:
			switch (index)
			{
			case 0:
				f = 0.5*(z2-z3)*(4*local(0,node)-1)/volume_;
				break;
			case 1:
				f = 0.5*(z3-z1)*(4*local(1,node)-1)/volume_;
				break;
			case 2:
				f = 0.5*(z1-z2)*(4*local(2,node)-1)/volume_;
				break;
			case 3:
				f = 2*((z2-z3)*local(1,node)+(z3-z1)*local(0,node))/volume_;
				break;
			case 4:
				f = 2*((z3-z1)*local(2,node)+(z1-z2)*local(1,node))/volume_;
				break;
			case 5:
				f = 2*((z1-z2)*local(0,node)+(z2-z3)*local(2,node))/volume_;
				break;
			}
			break;
		case 1:
			switch (index)		
			{
			case 0:
				f = 0.5*(r3-r2)*(4*local(0,node)-1)/volume_;
				break;
			case 1:
				f = 0.5*(r1-r3)*(4*local(1,node)-1)/volume_;
				break;
			case 2:
				f = 0.5*(r2-r1)*(4*local(2,node)-1)/volume_;
				break;
			case 3:
				f = 2*((r3-r2)*local(1,node)+(r1-r3)*local(0,node))/volume_;
				break;			
			case 4:
				f = 2*((r1-r3)*local(2,node)+(r2-r1)*local(1,node))/volume_;
				break;
			case 5:
				f = 2*((r2-r1)*local(0,node)+(r3-r2)*local(2,node))/volume_;
				break;			
			}
			break;
		case 2:
			return 0; // f = 0
			break;
		}
		return f;
	}
	
	void calculate_volume()
	{
		const value_type r1 = nodes_[0].dof[0], z1 = nodes_[0].dof[1];
		const value_type r2 = nodes_[1].dof[0], z2 = nodes_[1].dof[1];
		const value_type r3 = nodes_[2].dof[0], z3 = nodes_[2].dof[1];
		Parent::volume_ = 0.5*((r2-r1)*z3+(r1-r3)*z2+(r3-r2)*z1);
	}
protected:
}; 

} // namespace triangle
} // namespace axisymmetric


namespace plane
{
namespace triangle 
{
// 6-nodes triangle in plane case
// wherefore voigt number is 3:
// dim(Sxx, Syy,Sxy) = 3
class element6 : public axisymmetric::triangle::element6
{
public:
	enum
	{
		VoigtNumber = 3
	};
	typedef axisymmetric::triangle::element6 Parent;
public:
	element6() : Parent()
	{
	}
	element6(const NodeT nodes[])
	{
		for ( size_type i = 0; i < NodesNumber; ++ i )
			nodes_[i] = nodes[i];
		calculate();
	}
	element6(const element6& el)
	{
		copy(el);
	}


protected:
}; 

// 3-node triangle in plane case
// wherefore voigt number is 3:
// dim(Sxx, Syy,Sxy) = 3
class element3 : public ::element<3,2,3>
{
public:
	typedef ::element<3,2,3> Parent;
public:

	element3(const element3& el)
	{
		copy(el);
	}

	// local(L) coordinates
	// @param index - index of local function 0..2
	value_type local(size_type index, const NodeT& node) const
	{
		assert(index >= 0 && index < 3);
		value_type l = 0;
		
		const value_type x1 = nodes_[0].dof[0], y1 = nodes_[0].dof[1];
		const value_type x2 = nodes_[1].dof[0], y2 = nodes_[1].dof[1];
		const value_type x3 = nodes_[2].dof[0], y3 = nodes_[2].dof[1];

		value_type a[3];
		a[0] = x2*y3-x3*y2; // el(2,1)*el(3,2)-el(3,1)*el(2,2);
		a[1] = x3*y1-x1*y3; // el(3,1)*el(1,2)-el(1,1)*el(3,2);
		a[2] = x1*y2-x2*y1; // el(1,1)*el(2,2)-el(2,1)*el(1,2);
		
		value_type b[3];
		
		b[0] = y2-y3; // el(2,2)-el(3,2);
		b[1] = y3-y1; // el(3,2)-el(1,2);
		b[2] = y1-y2; // el(1,2)-el(2,2);

		value_type c[3];
		c[0] = x3-x2; // el(3,1)-el(2,1);
		c[1] = x1-x3; // el(1,1)-el(3,1);
		c[2] = x2-x1; // el(2,1)-el(1,1);

		l = (a[index]+b[index]*node.dof[0]+c[index]*node.dof[1])/(2.*volume_);;

		return l;
	}

	// form function
	// @param index - index of form function[0-(number_of_nodes-1)]
	// @param node - point to calculate form function
	value_type form(size_type index, const NodeT& node) const
	{
		assert(index >= 0 && index < 3);
		value_type f = local(index,node);
		return f;
	};
	
	// function for calculation derivatives of form functions
	// @param index - index of form function
	// @param dof - index of degree of freedom to integrate by
	value_type dform(size_type index, size_type dof, const NodeT& node) const
	{
		assert(index >= 0 && index < 3);
		assert(dof >=0 && dof < MAX_DOF);
		
		const value_type x1 = nodes_[0].dof[0], y1 = nodes_[0].dof[1];
		const value_type x2 = nodes_[1].dof[0], y2 = nodes_[1].dof[1];
		const value_type x3 = nodes_[2].dof[0], y3 = nodes_[2].dof[1];
		
		value_type f = 0;

		switch (dof)
		{
		case 0:
			switch (index)
			{
			case 0:
				f = 0.5*(y2-y3)/volume_;
				break;
			case 1:
				f = 0.5*(y3-y1)/volume_;
				break;
			case 2:
				f = 0.5*(y1-y2)/volume_;
				break;
			}
			break;
		case 1:
			switch (index)		
			{
			case 0:
				f = 0.5*(x3-x2)/volume_;
				break;
			case 1:
				f = 0.5*(x1-x3)/volume_;
				break;
			case 2:
				f = 0.5*(x2-x1)/volume_;
				break;
			}
			break;
		default:
			assert(false);
			break;
		}
		return f;
	}
	
	void calculate_volume()
	{
		const value_type r1 = nodes_[0].dof[0], z1 = nodes_[0].dof[1];
		const value_type r2 = nodes_[1].dof[0], z2 = nodes_[1].dof[1];
		const value_type r3 = nodes_[2].dof[0], z3 = nodes_[2].dof[1];
		Parent::volume_ = 0.5*((r2-r1)*z3+(r1-r3)*z2+(r3-r2)*z1);
	}
protected:
}; 


} // namespace triangle
} // namespace plane


// Tests section
/*
void test_element_axisymmetric_triangle_p6();
*/
void test_element_plane_triangle_p6();



#endif // __ELEMENTS_H__
