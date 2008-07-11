#ifndef __JOBTYPE_H__
#define __JOBTYPE_H__

#include "global.h"

namespace job
{
	// typenames for static polymorphysm
	namespace plane_strain
	{
		enum 
		{
			Voigt = 3
		};
		template<typename ElementT>
		value_type b_matrix_proxy(const ElementT& el,size_type i,size_type j)
		{
			size_type J = j >> 1;
			if ((j&1) == 0)
			{
				switch(i)
				{
				case 0:
					return el.gauss_derivative(gauss,J,0);
				case 1: 
					return 0;
				case 2:
					return el.gauss_derivative(gauss,J,1);
				};
			}
			else
			{
				switch(i)
				{
				case 0:
					return 0;
				case 1: 
					return el.gauss_derivative(gauss,J,1);
				case 2:
					return el.gauss_derivative(gauss,J,0);
				};
			}
			assert(false); // error case
			return 0;	
		}
	}
	namespace axisymmetric
	{
		enum
		{
			Voigt = 4
		};
		template<typename ElementT>
		value_type b_matrix_proxy(const ElementT& el,size_type i,size_type j)
		{
			size_type J = j >> 1;
			if ((j&1) == 0)
			{
				switch(i)
				{
				case 0:
					return el.gauss_derivative(gauss,J,0);
				case 1: 
					return 0;
				case 2:
					return el.gauss_form(gauss,J)/el.gauss_node(gauss).dof[0];
				case 3:
					return 0.5*el.gauss_derivative(gauss,J,1);
				};
			}
			else
			{
				switch(i)
				{
				case 0:
					return 0;
				case 1: 
					return el.gauss_derivative(gauss,J,1);
				case 2:
					return 0;
				default: 
					return 0.5*el.gauss_derivative(gauss,J,0);
				};
			}
			return 0;
		}
	}
	namespace three_dimension
	{
		enum
		{
			Voigt = 6
		};
	}

	// enum for easy switching btw models
	enum type
	{
		PlaneStrain = job::plane_strain::Voigt,
		Axisymmetric = job::axisymmetric::Voigt,
		ThreeDimension = job::three_dimension::Voigt
	};

}

#endif // __JOBTYPE_H__
