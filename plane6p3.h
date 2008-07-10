#ifndef __PLANE6P3_H__
#define __PLANE6P3_H__

// precompiled headers
#include "std.h"

// local includes
#include "genericsolver.h"
#include "geometrydatawrapper.h"

class plane6p3 : 
	public generic_solver<plane::triangle::element6,
						  LGeometryData,
						  triangle::gauss_tri7p>
{
public:
	typedef generic_solver<plane::triangle::element6,LGeometryData,triangle::gauss_tri7p> Parent;
	typedef geometrydata_wrapper<plane::triangle::element6> DataLoaderWrapper;
	using Parent::Element;
	using Parent::ElementsArray;
	using Parent::DataLoader;
	using Parent::GaussNodes;
public:
	plane6p3 (Model* model) : Parent(model)
	{
		data_.reset(new DataLoaderWrapper);
		data_->reset(new DataLoader);
	}
protected:
	// overloaded functions
	virtual void postprocessing(const char* filename,const ElementsArray& elements) const;
	inline value_type b_matrix_proxy(size_type elem, size_type gauss, size_type i, size_type j) const;
	inline value_type b_vector_proxy(size_type el, size_type gauss, size_type i, size_type j) const;

private:
	void export_msh(const char* filename,const ElementsArray& elements) const;
	void get_local_displacements(size_type el,const ElementsArray& elements,Vector& u) const;
};

inline value_type plane6p3::b_matrix_proxy(size_type el, size_type gauss, size_type i, size_type j) const
{
	// TODO: implement plane6p3::b_matrix_proxy
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
			return gauss_derivatives_[el][gauss][J][2];
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
			return gauss_derivatives_[el][gauss][J][1];
		};
	}
	assert(false);
	return 0;
}


inline value_type plane6p3::b_vector_proxy(size_type el, size_type gauss, size_type i, size_type j) const
{
	// TODO: implement plane6p3::b_vector_proxy
	switch (i)
	{
	case 0:
		return gauss_derivatives_[el][gauss][j][1];
	case 1:
		return gauss_derivatives_[el][gauss][j][2];
	};
	assert(false); // only 2 values in vector b
	return 0;
}


void plane6p3::postprocessing(const char* filename,const ElementsArray& elements) const
{
	export_msh(filename,elements);
}
	

void plane6p3::export_msh(const char* filename,const ElementsArray& elements) const
{
	assert(initial_elements_.size() == elements.size());

	boost::array<Tensor4Rank::index, 4> shape = {{ MAX_DOF, MAX_DOF, MAX_DOF, MAX_DOF }};
	Tensor4Rank elasticity_tensor(shape);
	MATRIX(C,3,3);
	construct_elasticity_matrix(C,elasticity_tensor);


	std::ofstream f(filename, std::ios_base::out | std::ios_base::trunc );
	if ( f )
	{
		for ( size_type i = 0; i < initial_elements_.size(); ++ i )
		{
			size_type j = 0;
			size_type k = 0;
			for ( ; j < 3; ++ j )
			{
				for ( k = 0; k < 2; ++ k )
				{
					f << initial_elements_[i].node(j).dof[k] << ",";
				}
			}
			for ( j = 0; j < 3; ++ j )
			{
				for ( k = 0; k < 2; ++ k )
				{
					f << elements[i].node(j).dof[k] << ",";
				}
			}
			f << "0,";
			Node nod = initial_elements_[i].center();
//			MATRIX(F,3,3);
//			deformation_gradient(i,elements,nod,F);
//			MATRIX(S,3,3);
//			model_->stress_cauchy(F,S);
//			f << S(0,0) << "," << S(1,1) << "," << 0 << "," << S(0,1) << std::endl;		
			VECTOR(u,12);
			get_local_displacements(i,elements,u);
			std::cout << "u " << u << std::endl;
			MATRIX(B,3,12);
			for ( j = 0; j < 3; ++ j )
				for ( k = 0; k < 12; ++ k )
					B(j,k) = b_matrix_proxy(i,0,j,k);

			VECTOR(sigma,3);
			Vector strain = prod(B,u);
			std::cout << "strain " << strain << std::endl;
			sigma = prod(C,strain);
			f << sigma[0] << "," << sigma[1] << "," << "0" << "," << sigma[2] << std::endl;		
		}
		f.flush();
		f.close();
	}
}

void plane6p3::get_local_displacements(size_type el,const ElementsArray& elements,Vector& u) const
{
	for ( size_type i = 0; i < 6; ++ i )
	{
		for ( size_type j = 0; j < 2; ++ j )
		{
			size_type idx = i*2+j;
			u[idx] = elements[el].node(i).dof[j] - initial_elements_[el].node(i).dof[j];
		}
	}
}

#endif
