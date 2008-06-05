#ifndef __AXISYMMETRIC6P3_H__
#define __AXISYMMETRIC6P3_H__

// precompiled headers
#include "std.h"

// local includes
#include "genericsolver.h"
#include "geometrydatawrapper.h"

class axisymmetric6p3 : 
	public generic_solver<axisymmetric::triangle::element6,
													LGeometryData,
													triangle::gauss_tri7p>
{
public:
	typedef generic_solver<axisymmetric::triangle::element6,LGeometryData,triangle::gauss_tri7p> Parent;
	typedef geometrydata_wrapper<axisymmetric::triangle::element6> DataLoaderWrapper;
public:
	axisymmetric6p3 (Model* model) : Parent(model)
	{
		data_.reset(new DataLoaderWrapper);
		data_->reset(new DataLoader);
	}
protected:
	// overloaded functions
	virtual void linear_construct_local_matrix(size_type elem, Matrix& m) const;
	inline value_type b_matrix_proxy(size_type elem, size_type gauss, size_type i, size_type j) const;
	inline value_type b_vector_proxy(size_type el, size_type gauss, size_type i, size_type j) const; 
	virtual void deformation_gradient(size_type el,const Node& node, Matrix& grad) const;
	virtual void deformation_gradient(size_type el,const ElementsArray& elements,const Node& node, Matrix& grad) const;
	virtual void construct_ksigma(size_type el, Matrix& m) const;
	virtual void construct_residual_force(size_type el,Vector& v) const;
	virtual void postprocessing(const char* filename,const ElementsArray& elements) const;
private:
	void export_msh(const char* filename,const ElementsArray& elements) const;
};

inline value_type axisymmetric6p3::b_matrix_proxy(size_type el, size_type gauss, size_type i, size_type j) const
{
	size_type J = j >> 1;
	if ((j&1) == 0)
	{
		switch(i)
		{
		case 0:
			return gauss_derivatives_[el][gauss][J][1];
		case 1: 
			return 0;
		case 2:
			return gauss_derivatives_[el][gauss][J][0]/gauss_nodes_[el][gauss].dof[0];
		default: 
			return 0.5*gauss_derivatives_[el][gauss][J][2];
		};
	}
	else
	{
		switch(i)
		{
		case 0:
			return 0;
		case 1: 
			return gauss_derivatives_[el][gauss][J][2];
		case 2:
			return 0;
		default: 
			return 0.5*gauss_derivatives_[el][gauss][J][1];
		};
	}

	return 0;
}


inline value_type axisymmetric6p3::b_vector_proxy(size_type el, size_type gauss, size_type i, size_type j) const
{
	switch (i)
	{
	case 0:
		return gauss_derivatives_[el][gauss][j][1];
	case 1:
		return gauss_derivatives_[el][gauss][j][2];
	default:
		return gauss_derivatives_[el][gauss][j][0]/gauss_nodes_[el][gauss].dof[0];
	};
	return 0;
}


#endif // __AXISYMMETRIC6P3_H__
