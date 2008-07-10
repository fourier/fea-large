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
	};
	namespace axisymmetric
	{
		enum
		{
			Voigt = 4
		};
	};
	namespace three_dimension
	{
		enum
		{
			Voigt = 6
		};
	};

	// enum for easy switching btw models
	enum type
	{
		PlaneStrain = job::plane_strain::Voigt,
		Axisymmetric = job::axisymmetric::Voigt,
		ThreeDimension = job::three_dimension::Voigt
	};

};

#endif // __JOBTYPE_H__
