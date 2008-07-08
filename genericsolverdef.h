#ifndef __GENERICSOLVERDEF_H__
#define __GENERICSOLVERDEF_H__

// pragmas
#pragma warning( disable : 4267 )
#pragma once

// precompiled headers
#include "std.h"

// local includes
#include "application.h"
#include "GeometryData.h"
#include "gauss_coeffs.h"
#include "elements.h"
#include "functions.h"
#include "models.h"

// macros and constants

/////////////////////////////////////////////////////////////////////////////



template<typename E,typename D,typename G>
generic_solver<E,D,G>::generic_solver(Model* model) : 
	data_(NULL),
	msize_(0),
	model_(model)
{
}

template<typename E,typename D,typename G>
generic_solver<E,D,G>::~generic_solver()
{
}

template<typename E,typename D,typename G>
bool generic_solver<E,D,G>::load_file(const char* filename)
{
	if ( data_->load_file(filename) )
	{
		prepare_task();
		return true;
	}
	return false;
}


template<typename E,typename D,typename G>
void generic_solver<E,D,G>::prepare_task()
{
	// 1st: fill array of elements
	convert_geometry_to_elements(initial_elements_ );
	// prepare index matrix using data wrapper
	prepare_index();

	// 2nd: presize necessary arrays: global matrix and vector, multiarrays of 
	// gauss nodes and derivatives in nodes
	// presize global matrix and global right part vector and set appropriate msize_ variable
	presize_matrix();
	// presize gauss nodes multiarray and form functions derivatives in gauss nodes multiarray
	presize_gauss_nodes();

	// 3rd: fill gauss nodes multiarray and form functions derivatives in gauss nodes multiarray
	prepare_gauss_nodes(initial_elements_); // at this point elements_ == initial_elements_
	// ... and store initial gauss nodes and derivatives of form functions in gauss nodes
	initial_gauss_nodes_ = gauss_nodes_;
	initial_gauss_derivatives_ = gauss_derivatives_;
	// fill array of deformation gradients
	update_graddef_array(initial_elements_);
}



template<typename E,typename D,typename G>
void generic_solver<E,D,G>::presize_matrix()
{
	msize_ = data_->nodes_vector().size()*2;
	global_matrix_ = boost::numeric::ublas::zero_matrix<size_type>(msize_,msize_);
	global_vector_.resize(msize_);global_vector_.clear();
	global_solution_.resize(msize_);global_solution_.clear();
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::presize_gauss_nodes()
{
	// resize multiarray of derivatives of gauss nodes
	// sizes: 
	// [Number of elements][Number of Gauss nodes][number of form functions][number of DoF + 1]
	// needed to construct the deformation gradients at any point
	// So, for every gauss node:
	// | N1       N2       ...  Nk      |
	// | dN1/dx1  dN2/dx1  ...  dNk/dx1 | x gauss nodes per element
	// | dN1/dx2  dN2/dx2  ...  dNk/dx2 | 
	// | dN1/dx3  dN2/dx3  ...  dNk/dx3 | 
	// -> last indexes [0] - is form function for this node, >0 - derivative of this form function
	// Thats why last size is [number of DoF + 1]
	boost::array<typename MultiArray::index,4> dims = 
		{{initial_elements_.size(), GaussNodes::GaussNumber, Element::NodesNumber, Element::DofNumber + 1}};       
	
	gauss_derivatives_.resize(dims);
	initial_gauss_derivatives_.resize(dims);

	boost::array<typename GaussArray::index,2> dims1 = 
		{{initial_elements_.size(),GaussNodes::GaussNumber}};
	gauss_nodes_.resize(dims1);
	initial_gauss_nodes_.resize(dims1);
	graddefs_gauss_nodes_.resize(dims1);
}

template<typename E,typename D,typename G>
inline value_type generic_solver<E,D,G>::symm_tensor4rank_matrix_proxy(const Tensor4Rank& c_tensor,size_type i, size_type j) const
{
	// Mapping according to Voight notation
	// please take a look at constructor generic_solver
	size_type k = g_voigt_mapping.find(i)->second.first;
	size_type l = g_voigt_mapping.find(i)->second.second;
	size_type m = g_voigt_mapping.find(j)->second.first;
	size_type n = g_voigt_mapping.find(j)->second.second;
	return c_tensor[k][l][m][n];
}

template<typename E,typename D,typename G>
inline void generic_solver<E,D,G>::construct_elasticity_matrix(Matrix& m,const Tensor4Rank& c_tensor) const
{
/*
	const size_type voigt = Element::VoigtNumber;
	for ( size_type i = 0; i < voigt; ++ i )
	{
		for ( size_type j = 0; j < voigt; ++ j )
		{
			m(i,j) = symm_tensor4rank_matrix_proxy(c_tensor,i,j);
		}
	}
*/

	double nu = 0.3;
	double E = 1e9;
/*
	m(0,0) = 1;  m(0,1) = nu; 
	m(1,0) = nu; m(1,1) = 1;
							  m(2,2) = 0.5*(1-nu);
	m = m*E/(1-nu*nu);
*/
/*
	M(1,1) = E/(1-nu^2);
	M(2,2) = M(1,1);
	M(3,3) = 0.5*E/(1+nu);
	M(1,2) = nu*M(1,1);
	M(2,1) = M(1,2);

*/
	m(0,0) = E/(1-nu*nu);
	m(1,1) = m(0,0);
	m(2,2) = 0.5*E/(1+nu);
	m(0,1) = nu*m(0,0);
	m(1,0) = m(0,1);
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::prepare_gauss_nodes(const ElementsArray& elements)
{
	// fill multiarray with derivatives
	for ( size_type el = 0; el < elements.size(); ++ el )
	{
		const Element& elem = elements[el];
		for ( size_type gauss = 0; gauss < GaussNodes::GaussNumber; ++ gauss )
		{
			Node node = elem.template gauss_point<G>(gauss);
			std::cout << node.dof[0] << " " << node.dof[1] << std::endl;
			gauss_nodes_[el][gauss] = node;
			for ( size_type shape = 0; shape < Element::NodesNumber; ++ shape )
			{
				// first, put value of form function into the array
				gauss_derivatives_[el][gauss][shape][0] = elem.form(shape,node); 
				// next, put all derivatives of this form function into the array
				for ( size_type dof = 0; dof < Element::DofNumber; ++ dof )
					gauss_derivatives_[el][gauss][shape][dof+1] = elem.dform(shape,dof,node);
			}
		}	
	}
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::convert_geometry_to_elements(ElementsArray& elements)
{
	data_->load_elements(elements);
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::prepare_index()
{
	index_.clear();

	size_type rows = initial_elements_.size();
	size_type cols = Element::NodesNumber*Element::DofNumber;

	index_ = boost::numeric::ublas::zero_matrix<size_type>(rows,cols);
	
	index_matrix_proxy<IndexMatrix> m(index_);
	data_->prepare_index(m);
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::linear_construct_local_matrix(size_type el, Matrix& m) const
{
//	TODO: remove this
	linear_construct_local_matrix1(el,m);
	return;

	assert(m.size1() == m.size2());
	size_type size = m.size1();
	
	boost::array<Tensor4Rank::index, 4> shape = {{ MAX_DOF, MAX_DOF, MAX_DOF, MAX_DOF }};
	Tensor4Rank elasticity_tensor(shape);

	// first, iterate over gauss nodes
	for ( size_type gauss = 0; gauss < G::GaussNumber; ++ gauss )
	{
		// multiplication rule:
		// R = A*B*C;
		// R_il = A*B*C = A_ij * B_jk * C_kl
		// in our case 
		// K = B^T*C*B
		// so K_il = B_ji * C_jk * B_kl

		// First of all constuct [B] matrix
		const size_type voigt = Element::VoigtNumber;

		MATRIX(B,voigt,size);
		for ( size_type i = 0; i < voigt; ++ i )
			for ( size_type j = 0; j < size; ++ j )
				B(i,j) = b_matrix_proxy(el,gauss,i,j);

		// next construct C tensor in gauss node
		model_->construct_ctensor(elasticity_tensor,graddefs_gauss_nodes_[el][gauss]);
		// constuct symmetric matrix from tensor
		MATRIX(C,voigt,voigt);
		construct_elasticity_matrix(C,elasticity_tensor);

		// next calculate inversion of det(F)
#ifdef DETF
		const value_type inv_detF = 1.0/det3x3(graddefs_gauss_nodes_[el][gauss]);
		// next assemble B^T * C * B
		// optimization: move out multiplier
		const value_type multiplier = GaussNodes::weights[gauss]*inv_detF;
#else
		const value_type multiplier = GaussNodes::weights[gauss];	
#endif
		for ( size_type i = 0; i < size; ++ i )
		{
			for ( size_type j = 0; j < voigt; ++ j )
			{
				for ( size_type k = 0; k < voigt; ++ k )
				{
					for ( size_type l = 0; l < size; ++ l )
					{
						// K_il = B_ji * C_jk * B_kl
//						m(i,l) += B(j,i)*symm_tensor4rank_matrix_proxy(elasticity_tensor,j,k)*B(k,l)*multiplier;
						m(i,l) += B(j,i)*C(j,k)*B(k,l)*multiplier;
					}
				}
			}
		}
	}
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::linear_construct_local_matrix1(size_type el, Matrix& m) const
{
	boost::array<Tensor4Rank::index, 4> shape = {{ MAX_DOF, MAX_DOF, MAX_DOF, MAX_DOF }};
	Tensor4Rank elasticity_tensor(shape);
	size_type size = m.size1();
	// first, iterate over gauss nodes
	for ( size_type gauss = 0; gauss < G::GaussNumber; ++ gauss )
	{
		const size_type voigt = Element::VoigtNumber;

		MATRIX(B,voigt,size);
		for ( size_type i = 0; i < voigt; ++ i )
			for ( size_type j = 0; j < size; ++ j )
				B(i,j) = b_matrix_proxy(el,gauss,i,j);
		char fname[50];
		sprintf(fname,"Dump/B%d_%d.txt",el+1,gauss+1);
		Dump(B,fname);
		// constuct symmetric matrix from tensor
		MATRIX(C,voigt,voigt);
		construct_elasticity_matrix(C,elasticity_tensor);
		Dump(C,"Dump/C.txt");
		const value_type multiplier = GaussNodes::weights[gauss];	
		for ( size_type i = 0; i < size; ++ i )
		{
			for ( size_type j = 0; j < voigt; ++ j )
			{
				for ( size_type k = 0; k < voigt; ++ k )
				{
					for ( size_type l = 0; l < size; ++ l )
					{
						// K_il = B_ji * C_jk * B_kl
						m(i,l) += B(j,i)*C(j,k)*B(k,l)*multiplier;
					}
				}
			}
		}
	}


}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::distribute_local_in_global_matrix(size_type element_index, const Matrix& local)
{
	const size_type i = element_index; // element index
	for ( size_type k = 0; k < local.size1(); ++ k )
		for (size_type l = 0; l < local.size2(); ++ l )
			global_matrix_(index_(i,k),index_(i,l)) = global_matrix_(index_(i,k),index_(i,l)) + local(k,l);
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::distribute_local_in_global_vector(size_type el, const Vector& local)
{
	for ( size_type k = 0; k < local.size(); ++ k )
		global_vector_(index_(el,k)) += local[k];
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::linear_construct_global_matrix(const ElementsArray& elements)
{
	boost::array<Tensor4Rank::index, 4> shape = {{ MAX_DOF, MAX_DOF, MAX_DOF, MAX_DOF }};
	Tensor4Rank elasticity_tensor(shape);
	const size_type voigt = Element::VoigtNumber;
	MATRIX(C,voigt,voigt);
	construct_elasticity_matrix(C,elasticity_tensor);
	for ( int i = 0; i < voigt; ++ i )
	{
		for (int j = 0; j < voigt; ++ j )
		{
			std::cout << C(i,j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << "next" << std::endl;
	model_->construct_ctensor(elasticity_tensor,C);
	for ( int i = 0; i < voigt; ++ i )
	{
		for (int j = 0; j < voigt; ++ j )
		{
			std::cout << symm_tensor4rank_matrix_proxy(elasticity_tensor,i,j) << " ";
		}
		std::cout << std::endl;
	}	

	// iterate through all elements
	for ( size_type el = 0; el < elements.size(); ++ el  )
	{
		const Element& elem = elements[el];
		const value_type vol = elem.volume();
		const size_type size = Element::NodesNumber*Element::DofNumber;
		MATRIX(local,size,size);
		// prepare local matrix
		linear_construct_local_matrix(el,local);
		local = local*vol;
		char fname[50];
		sprintf(fname,"Dump/K%d.txt",el+1);
		Dump(local,fname);
		// distribute local matrix in global matrix
		distribute_local_in_global_matrix(el,local);
	}
	Dump(global_matrix_,"Dump/global.txt");
}


template<typename E,typename D,typename G>
void generic_solver<E,D,G>::apply_bc_prescribed(bool fixed_all)
{
	const typename DataLoaderWrapper::PrescribedBoundaryArray& boundary = data_->prescribed_boundary_array();
	// apply x(r for axisymmetric)-coordinate conditions
	typename DataLoaderWrapper::PrescribedBoundaryArray::const_iterator it = boundary.begin();
	for ( ; it != boundary.end(); ++ it )
	{
		if ( it->template get<1>() == DataLoaderWrapper::EPrescribedX || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedXZ  || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedXY  ||
			 it->template get<1>() == DataLoaderWrapper::EPrescribedXYZ )
		{
			size_type index = it->template get<0>()*Element::DofNumber; // r-coordinate
			Matrix::value_type tmp = global_matrix_(index,index);
			Matrix::value_type condition = fixed_all ? 0 : it->template get<2>().dof[0];
			for ( size_type j = 0; j < msize_; ++ j )
			{
				global_vector_[j] = global_vector_[j] - global_matrix_(index,j)*condition;
				global_matrix_(index,j) = 0;
				global_matrix_(j,index) = 0;
			}
			global_matrix_(index,index) = tmp;
			global_vector_[index] = tmp*condition;
		}
	}
	// apply y(z for axisymmetric)-coordinate conditions
	it = boundary.begin();
	for ( ; it != boundary.end(); ++ it )
	{
		if ( it->template get<1>() == DataLoaderWrapper::EPrescribedY || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedXY || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedYZ || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedXYZ )
		{
			size_type index = it->template get<0>()*Element::DofNumber + 1; // z-coordinate
			Matrix::value_type tmp = global_matrix_(index,index);
			Matrix::value_type condition = fixed_all ? 0 : it->template get<2>().dof[1];
			for ( size_type j = 0; j < msize_; ++ j )
			{
				global_vector_[j] = global_vector_[j] - global_matrix_(index,j)*condition;
				global_matrix_(index,j) = 0;
				global_matrix_(j,index) = 0;
			}
			global_matrix_(index,index) = tmp;
			global_vector_[index] = tmp*condition;
		}
	}
	// apply z(for 3d case)-coordinate conditions
	it = boundary.begin();
	for ( ; it != boundary.end(); ++ it )
	{
		if ( it->template get<1>() == DataLoaderWrapper::EPrescribedZ || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedXZ || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedYZ || 
			 it->template get<1>() == DataLoaderWrapper::EPrescribedXYZ )
		{
			size_type index = it->template get<0>()*Element::DofNumber + 2; // z-coordinate
			Matrix::value_type tmp = global_matrix_(index,index);
			Matrix::value_type condition = fixed_all ? 0 : it->template get<2>().dof[2];
			for ( size_type j = 0; j < msize_; ++ j )
			{
				global_vector_[j] = global_vector_[j] - global_matrix_(index,j)*condition;
				global_matrix_(index,j) = 0;
				global_matrix_(j,index) = 0;
			}
			global_matrix_(index,index) = tmp;
			global_vector_[index] = tmp*condition;
		}
	}
}



template<typename E,typename D,typename G>
void generic_solver<E,D,G>::construct_global_stiffness_matrix(const ElementsArray& elements)
{
	global_matrix_.clear();
	// iterate through all elements
	const size_type size = Element::NodesNumber*Element::DofNumber;
	for ( size_type el = 0; el < elements.size(); ++ el  )
	{
		const Element& elem = elements[el];
		const value_type vol = elem.volume();
		
		MATRIX(kconst,size,size);
		// prepare local constitutive matrix
		construct_kconstr(el,kconst);
		kconst = kconst*vol;
		
		// prepare matrix of initial stiffness
		MATRIX(ksigma,size,size);
		construct_ksigma(el,ksigma);
		ksigma = ksigma*vol;

		// distribute local matrixes in global matrix
		distribute_local_in_global_matrix(el,ksigma);
		distribute_local_in_global_matrix(el,kconst);
	}
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::construct_global_residual_force_vector(const ElementsArray& elements)
{
	global_vector_.clear();
	const size_type size = Element::NodesNumber*Element::DofNumber;
	for ( size_type el = 0; el < elements.size(); ++ el  )
	{
		const Element& elem = elements[el];
		const value_type vol = elem.volume();
		
		VECTOR(v,size);
		construct_residual_force(el,v);
		v = v*vol;
		distribute_local_in_global_vector(el,v);
	}
}


template<typename E,typename D,typename G>
void generic_solver<E,D,G>::update_with_prescribed_conditions(ElementsArray& elements)
{
	const typename DataLoaderWrapper::PrescribedBoundaryArray& boundary = data_->prescribed_boundary_array();
	// define temporary 'solution-like' vector with displacements from boundary conditions
	VECTOR(X,msize_);
	typename DataLoaderWrapper::PrescribedBoundaryArray::const_iterator it = boundary.begin();
	for ( ; it != boundary.end(); ++ it )
	{
		size_type index = it->template get<0>(); // index of node in nodes array
		size_type global_index = index*Element::DofNumber;
		switch (it->template get<1>())
		{
		case DataLoaderWrapper::EPrescribedX:
			X[global_index] = it->template get<2>().dof[0];
			break;
		case DataLoaderWrapper::EPrescribedY:
			X[global_index+1] = it->template get<2>().dof[1];
			break;
		case DataLoaderWrapper::EPrescribedZ:
			X[global_index+2] = it->template get<2>().dof[2];
			break;
		case DataLoaderWrapper::EPrescribedXY:
			X[global_index] = it->template get<2>().dof[0];
			X[global_index+1] = it->template get<2>().dof[1];
			break;
		case DataLoaderWrapper::EPrescribedXZ:
			X[global_index] = it->template get<2>().dof[0];
			X[global_index+2] = it->template get<2>().dof[2];
			break;
		case DataLoaderWrapper::EPrescribedYZ:
			X[global_index+1] = it->template get<2>().dof[1];
			X[global_index+2] = it->template get<2>().dof[2];
			break;
		case DataLoaderWrapper::EPrescribedXYZ:
			X[global_index] = it->template get<2>().dof[0];
			X[global_index+1] = it->template get<2>().dof[1];
			X[global_index+2] = it->template get<2>().dof[2];
			break;
		default:
			assert(false);
		}
	}
	update_with_solution(elements,X);
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::update_with_solution(ElementsArray& elements, const Vector& X) const
{
	for ( size_type i = 0; i < elements.size(); ++ i )
	{
		for ( size_type j = 0; j < Element::NodesNumber; ++ j )
		{
			// dof in global matrix	and vector
			for ( size_type k = 0; k < Element::DofNumber; ++ k )
			{
				size_type local_index = j*Element::DofNumber+k;
				size_type global_index = index_(i,local_index);
				elements[i].node(j).dof[k] += X[global_index];
			}
		}
		elements[i].calculate_volume();
	}
}


template<typename E,typename D,typename G>
void generic_solver<E,D,G>::update_graddef_array(ElementsArray& elements)
{
	for ( size_type el = 0; el < elements.size(); ++ el )
	{
		const Element& elem = elements[el];
		for ( size_type gauss = 0; gauss < GaussNodes::GaussNumber; ++ gauss )
		{
			MATRIX(F,MAX_DOF,MAX_DOF);
			deformation_gradient(el,elements,initial_gauss_nodes_[el][gauss],F);
			assert(det3x3(F) > 0);
			graddefs_gauss_nodes_[el][gauss] = F;
		}	
	}
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::construct_kconstr(size_type el, Matrix& m) const
{
	// TODO: carefully review construct kconstr
	linear_construct_local_matrix(el,m);
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::construct_ksigma(size_type el, Matrix& m) const
{
	using boost::numeric::ublas::prod;
	using boost::numeric::ublas::inner_prod;
	using boost::numeric::ublas::trans;

	assert(m.size1() == m.size2());
	const size_type size = m.size1();

	// first, iterate over gauss nodes.
	for ( size_type gauss = 0; gauss < GaussNodes::GaussNumber; ++ gauss )
	{
		MATRIX(Sigma,3,3);
		model_->stress_cauchy(graddefs_gauss_nodes_[el][gauss],Sigma);
		std::vector<Vector> B_vectors;
		B_vectors.resize(Element::NodesNumber);
		for ( size_type i = 0; i < Element::NodesNumber; ++ i )
		{
			B_vectors[i].resize(3);
			for ( size_type j = 0; j < 3; ++ j )
				B_vectors[i][j] = b_vector_proxy(el,gauss,j,i);
		}
		// optimization: move out multiplier
		// for axisymmetric task 2*PI*r
		const value_type multiplier = 2*PI*GaussNodes::weights[gauss];
		for ( size_type k = 0; k < Element::NodesNumber; ++ k )
		{
			for ( size_type l = 0; l < Element::NodesNumber; ++ l )
			{
				Vector tmp( prod(trans(B_vectors[k]),Sigma) );
				value_type value = inner_prod(tmp,B_vectors[l]);		
				const size_type offset_row = k*Element::DofNumber;
				const size_type offset_column = l*Element::DofNumber;
				value = value*multiplier;
				m(offset_row,offset_column) += value;	
				m(offset_row+1,offset_column+1) += value;
			}
		}
	}
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::construct_residual_force(size_type el,Vector& v) const
{
	using boost::numeric::ublas::prod;
	using boost::numeric::ublas::inner_prod;
	using boost::numeric::ublas::trans;

	size_type size = v.size();

	// first, iterate over gauss nodes
	for ( size_type gauss = 0; gauss < GaussNodes::GaussNumber; ++ gauss )
	{
		// First of all constuct [B] matrix
		MATRIX(B,Element::VoigtNumber,size);
			for ( size_type i = 0; i < Element::VoigtNumber; ++ i )
				for ( size_type j = 0; j < size; ++ j )
					B(i,j) = b_matrix_proxy(el,gauss,i,j);
		MATRIX(Sigma,Element::VoigtNumber,Element::VoigtNumber);
		model_->stress_cauchy(graddefs_gauss_nodes_[el][gauss],Sigma);
		VECTOR(S,Element::VoigtNumber);
		for ( size_type l = 0; l < Element::VoigtNumber; ++ l )
		{
			S[l] = Sigma(g_voigt_mapping.find(l)->second.first,g_voigt_mapping.find(l)->second.second);
		}
		v += GaussNodes::weights[gauss]*prod(trans(B),S);
	}
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::postprocessing(const char* filename,const ElementsArray& elements) const
{
	// TODO: add generic post processing
	assert(false);
}

template<typename E,typename D,typename G>
value_type generic_solver<E,D,G>::max_prescribed_bc(size_type dof_number) const
{
	value_type displ = 0;
	const typename DataLoaderWrapper::PrescribedBoundaryArray& boundary = data_->prescribed_boundary_array();
	// find x(r for axisymmetric)-coordinate conditions
	typename DataLoaderWrapper::PrescribedBoundaryArray::const_iterator it = boundary.begin();

	switch(dof_number)
	{
	case 0: // x(r):
		for ( ; it != boundary.end(); ++ it )
		{
			if ( it->template get<1>() == DataLoaderWrapper::EPrescribedX || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedXZ  || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedXY  ||
				it->template get<1>() == DataLoaderWrapper::EPrescribedXYZ )
			{
				Matrix::value_type condition = it->template get<2>().dof[0];
				if ( fabs(condition) > fabs(displ) )
					displ = condition;
			}
		}
		break;
	case 1:// y(z)
		for ( ; it != boundary.end(); ++ it )
		{
			if ( it->template get<1>() == DataLoaderWrapper::EPrescribedY || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedXY || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedYZ || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedXYZ )
			{
				Matrix::value_type condition = it->template get<2>().dof[1];
				if ( fabs(condition) > fabs(displ) )
					displ = condition;
			}
		}
		break;
	case 2: // z
		for ( ; it != boundary.end(); ++ it )
		{
			if ( it->template get<1>() == DataLoaderWrapper::EPrescribedZ || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedXZ || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedYZ || 
				it->template get<1>() == DataLoaderWrapper::EPrescribedXYZ )
			{
				Matrix::value_type condition = it->template get<2>().dof[2];
				if ( fabs(condition) > fabs(displ) )
					displ = condition;
			}
		}
		break;
	}
	return displ;
}


template<typename E,typename D,typename G>
void generic_solver<E,D,G>::deformation_gradient(size_type el,const Node& node, Matrix& F) const
{
	assert( F.size1() == F.size2() && F.size1() == MAX_DOF );
	F.clear();
	eye(F);
	
	const Element& elem = initial_elements_[el];
	for ( size_type i = 0; i < E::DofNumber; ++ i )
	{
		for ( size_type j = 0; j < E::DofNumber; ++ j )
		{
			value_type sum = 0;
			for ( size_type k = 0; k < E::NodesNumber; ++ k )
			{
				const value_type diff = elem.dform(k,j,node);

				const size_type dof_pos = k*E::DofNumber;
				const size_type index = index_(el,dof_pos+i);
				const value_type x = global_solution_[index];
				sum += diff*x;
			}
			F(i,j) += sum;
		}
	}
}

template<typename E,typename D,typename G>
void generic_solver<E,D,G>::deformation_gradient(size_type el,const ElementsArray& elements,const Node& node, Matrix& F) const
{
	assert( F.size1() == F.size2() && F.size1() == MAX_DOF );
	F.clear();
	eye(F);
	
	const Element& elem = initial_elements_[el];
	const Element& deformed = elements[el];
	
	for ( size_type i = 0; i < E::DofNumber; ++ i )
	{
		for ( size_type j = 0; j < E::DofNumber; ++ j )
		{
			value_type sum = 0;
			for ( size_type k = 0; k < E::NodesNumber; ++ k )
			{
				const value_type diff = elem.dform(k,j,node);
				const value_type x = deformed.node(k).dof[i] - elem.node(k).dof[i];
				sum += diff*x;
			}
			F(i,j) += sum;
		}
	}
}

#endif // __GENERICSOLVERDEF_H__


