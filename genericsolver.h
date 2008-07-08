#ifndef __GENERICSOLVER_H__
#define __GENERICSOLVER_H__

#pragma once

// precompiled headers
#include "std.h"

// local includes
#include "global.h"
#include "elements.h"
#include "dataloaderwrapper.h"
#include "models.h"
#include "LSLAE.h"


// using namespaces
namespace ublas = boost::numeric::ublas;

template<typename ElementT,typename DataLoaderT,typename GaussNodesT>
class generic_solver_base
{
public:
	// template typdefs
	typedef ElementT Element;
	typedef typename ElementT::NodeT Node;
	typedef std::vector<Element> ElementsArray;
	typedef DataLoaderT DataLoader;
	typedef GaussNodesT GaussNodes;
	typedef constitutive_equation_base<Matrix> Model;
	typedef std::auto_ptr<Model> ModelPtr;
	// internal  typedefs
	typedef boost::multi_array<Node,2> GaussArray;
	typedef boost::multi_array<value_type, 4> MultiArray;
	typedef boost::multi_array<Matrix,2> GradDefGaussArray;
	typedef data_loader_wrapper<DataLoader,Element,ElementsArray> DataLoaderWrapper;
public:
	virtual Matrix& global_matrix() = 0;
	virtual const Matrix& global_matrix() const = 0;
	virtual Vector& global_vector() = 0;
	virtual const Vector& global_vector() const = 0; 

	virtual const ElementsArray& elements() const = 0;
	virtual const Model& model() const = 0;

	// apply prescribed boundary conditions to the global matrix and global vector
	// fixed_all - make all prescribed displacements equal to 0
	// fixed_all must be set for nonlinear task 
	// used together with update_with_prescribed_conditions
	// virtual to make it possible to update boundary derived class
	virtual void apply_bc_prescribed(bool fixed_all = false) = 0;
	// distribute matrix from local matrix into the global stiffness matrix using index matrix
	virtual void distribute_local_in_global_matrix(size_type element_index, const Matrix& local) = 0;
	// distribute vector from local matrix into the global stiffness matrix using index matrix
	virtual void distribute_local_in_global_vector(size_type element_index, const Vector& local) = 0;
	// update deformation gradient in gauss nodes array with 
	// new calculated elements
	virtual void update_graddef_array(ElementsArray& elements) = 0;
	virtual void prepare_gauss_nodes(const ElementsArray& elements) = 0;
	virtual void update_with_prescribed_conditions(ElementsArray& elements) = 0;
	virtual void update_with_solution(ElementsArray& elements, const Vector& X) const = 0;

	// linear solution steps
	virtual void linear_construct_global_matrix(const ElementsArray& elements) = 0;
	virtual void construct_global_stiffness_matrix(const ElementsArray& elements) = 0;
	virtual void construct_global_residual_force_vector(const ElementsArray& elements) = 0;


	// matrixes
	virtual void construct_kconstr(size_type el, Matrix& m) const = 0;
	virtual void construct_ksigma(size_type el, Matrix& m) const = 0;

	virtual void construct_residual_force(size_type el,Vector& v) const = 0;

	virtual void postprocessing(const char* filename,const ElementsArray& elements) const = 0;

	virtual void deformation_gradient(size_type el,const ElementsArray& elements,const Node& node, Matrix& grad) const = 0;

	// returns maximum by absolute value prescribed displacement
	virtual value_type max_prescribed_bc(size_type dof_number) const = 0;

	// static
	// static solution of SLAE
	static Vector solve(const Matrix& m, const Vector& v)
	{
		VECTOR(x,v.size());
		try 
		{
//			x = LSLAE::gauss_solve_ex(m,v);
//			x = LSLAE::lu_substitute_solve(m,v);
			x = LSLAE::cg_solve(m,v,10000,1e-15);
			Vector residual = prod(m,x) - v;
			value_type norm_residual = boost::numeric::ublas::norm_2(residual);
			std::cout << "||A*x-b|| residual: " << norm_residual << std::endl;; 
		}
		catch ( std::logic_error& e )
		{
			std::cout << e.what() << std::endl;
			return x;
		}
		return x;
	}
};



template<typename ElementT,typename DataLoaderT,typename GaussNodesT>
class generic_solver : public generic_solver_base<ElementT, DataLoaderT, GaussNodesT>
{
public: // typedefs
	typedef generic_solver_base<ElementT, DataLoaderT, GaussNodesT> Parent;
	// template typdefs
	typedef typename Parent::Element Element;
	typedef typename Parent::Node Node;
	typedef typename Parent::ElementsArray ElementsArray;
	typedef typename Parent::DataLoader DataLoader;
	typedef typename Parent::GaussNodes GaussNodes;
	typedef typename Parent::Model Model;
	typedef typename Parent::ModelPtr ModelPtr;
	// internal  typedefs
	typedef typename Parent::GaussArray GaussArray;
	typedef typename Parent::MultiArray MultiArray;
	typedef typename Parent::GradDefGaussArray GradDefGaussArray;
	typedef typename Parent::DataLoaderWrapper DataLoaderWrapper;

public:
	template<typename T>
	struct index_matrix_proxy : public index_matrix_interface
	{
		index_matrix_proxy(T& m) : index_(m) {}
		size_type& operator () (size_type i, size_type j)
		{
			return index_(i,j);
		}
		T& index_;
	};
public: // enums
	enum TState
	{
		EFileLoaded,
		EUnknown
	};
public:
	generic_solver(Model* model);
	virtual ~generic_solver();
	virtual bool load_file(const char* filename);
public:
	virtual Matrix& global_matrix() { return global_matrix_; };
	virtual const Matrix& global_matrix() const { return global_matrix_; };
	virtual Vector& global_vector() { return global_vector_; };
	virtual const Vector& global_vector() const { return global_vector_; };
	virtual const Model& model() const { return *model_;};

	// apply prescribed boundary conditions to the global matrix and global vector
	// fixed_all - make all prescribed displacements equal to 0
	// fixed_all must be set for nonlinear task 
	// used together with update_with_prescribed_conditions
	// virtual to make it possible to update boundary derived class
	virtual void apply_bc_prescribed(bool fixed_all = false); 
	// distribute matrix from local matrix into the global stiffness matrix using index matrix
	virtual void distribute_local_in_global_matrix(size_type element_index, const Matrix& local);
	// distribute vector from local matrix into the global stiffness matrix using index matrix
	virtual void distribute_local_in_global_vector(size_type element_index, const Vector& local);
	// condition derivative solution steps
	// update deformation gradient in gauss nodes array with 
	// new calculated elements
	void update_graddef_array(ElementsArray& elements);
	// prepare gauss nodes array and gauss derivatives array
	// using given elements vector, to use in Newton iterations
	void prepare_gauss_nodes(const ElementsArray& elements);

	void update_with_prescribed_conditions(ElementsArray& elements);
	void update_with_solution(ElementsArray& elements, const Vector& X) const;

	void linear_construct_global_matrix(const ElementsArray& elements);
	void construct_global_stiffness_matrix(const ElementsArray& elements);
	void construct_global_residual_force_vector(const ElementsArray& elements);



	virtual void construct_kconstr(size_type el, Matrix& m) const;
	virtual void construct_ksigma(size_type el, Matrix& m) const;
	virtual void construct_residual_force(size_type el,Vector& v) const;

	virtual const ElementsArray& elements() const { return initial_elements_; };

	// returns maximum by absolute value prescribed displacement
	virtual value_type max_prescribed_bc(size_type dof_number) const;

	// perform post solution steps
protected:
	// all steps in one
	void prepare_task();
	// presize global matrix and vector
	void presize_matrix();
	// presize gauss nodes array and derivatives array
	void presize_gauss_nodes();
	// load elements from geometry data loader wrapper
	void convert_geometry_to_elements(ElementsArray& elements);
	// virtual steps - can be impemented in derived classes
	// prepare index matrix
	virtual void prepare_index(); // virtual to make it possible to update index in derived class
	virtual void postprocessing(const char* filename,const ElementsArray& elements) const;

	// local matrix may be constructed differently
	virtual void linear_construct_local_matrix(size_type elem, Matrix& m) const;
	void linear_construct_local_matrix1(size_type elem, Matrix& m) const;

	// auxulary functions
	// function which returns component of matrix with indexes (i,j) from 
	// symmetric tensor of 4th rank c_tensor using Voigt notation
	inline value_type symm_tensor4rank_matrix_proxy(const Tensor4Rank& c_tensor,size_type i, size_type j) const;
	// constructor for elasticity matrix
	inline void construct_elasticity_matrix(Matrix& m,const Tensor4Rank& c_tensor) const;
	// element of B matrix(matrix of derivatives of shape functions)
	// must be implemented in derived class
	// elem - number of element
	// gauss - number of gauss node
	// i,j - row and column indexes in matrix
	virtual value_type b_matrix_proxy(size_type elem, size_type gauss, size_type i, size_type j) const = 0;
	// vector of derivatives of shape functions
	// may be implemented in derived class
	// el - number of element
	// gauss - number of gauss node
	// i - row in vector
	// j - number of form function
	virtual value_type b_vector_proxy(size_type el, size_type gauss, size_type i, size_type j) const = 0;
	// deformation gradient in element with index el, in point 'node'	
	// node is the node in not deformed configuration
	virtual void deformation_gradient(size_type el,const Node& node, Matrix& grad) const;
	virtual void deformation_gradient(size_type el,const ElementsArray& elements,const Node& node, Matrix& grad) const;

protected:
	// auto pointer to a geometry data loaded from a file
	std::auto_ptr<DataLoaderWrapper> data_;
	// array of elements constructed from geometry data and updated with solution
//	ElementsArray elements_;
	// array of initial elements constructed from geometry data
	ElementsArray initial_elements_;
	// matrix of indexes for disctribution dof in global matrix
	IndexMatrix index_;
	// global matrix itself
	Matrix global_matrix_;
	// global vector of right part
	Vector global_vector_;
	// global vector of solution
	Vector global_solution_;

	// multidimensional array of derivatives in Gauss nodes of elements
	// number of Gauss nodes depends on the value of scheme_
	// access: 
	// gauss_derivatives_[Number of an element][Number of Gauss nodes][number of form function][form,derivative]
	// needed to construct the deformation gradients at any point
	// So, for every gauss node:
	// | N1       N2       ...  Nk      |
	// | dN1/dx1  dN2/dx1  ...  dNk/dx1 | x gauss nodes per element
	// | dN1/dx2  dN2/dx2  ...  dNk/dx2 | 
	// | dN1/dx3  dN2/dx3  ...  dNk/dx3 | 
	// 
	MultiArray gauss_derivatives_;
	// multiarray of initial gauss derivatives
	MultiArray initial_gauss_derivatives_;
	// multiarray of Gauss nodes per element 
	GaussArray gauss_nodes_;
	// multiarray of initial Gauss nodes per element 
	GaussArray initial_gauss_nodes_;
	// multiarray of deformation gradients in Gauss nodes per element 
	GradDefGaussArray graddefs_gauss_nodes_;
	// size of global matrix. Convinient to store it in an external variable
	size_type msize_;
	// Material model
	ModelPtr model_;
};

#include "genericsolverdef.h"

#endif // __GENERICSOLVER_H__
