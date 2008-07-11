#ifndef __MODELS_H__
#define __MODELS_H__

// precompiled headers
#include "std.h"

// local includes
#include "functions.h"
#include "jobtype.h"
#include "global.h"

// namespace for various voigt mappings
namespace mapping
{

// 
// Mapping from symmetric tensor of 4th rank to symmetric matrix indicies
// 
value_type tensor4rank_matrix_proxy(const Tensor4Rank& c_tensor,size_type i, size_type j) 
{
	// Mapping according to Voight notation
	// please take a look at constructor generic_solver
	size_type k = g_voigt_mapping.find(i)->second.first;
	size_type l = g_voigt_mapping.find(i)->second.second;
	size_type m = g_voigt_mapping.find(j)->second.first;
	size_type n = g_voigt_mapping.find(j)->second.second;
	return c_tensor[k][l][m][n];
}

// 
// Convert matrix to vector using Voigt notation
// 
template<int>
void vector_from_matrix(const Matrix& m,Vector& v);

template<>
void vector_from_matrix<job::PlaneStrain>(const Matrix& m,Vector& v)
{
	v.clear();
	v.resize(job::plane_strain::Voigt);
	v[0] = m(0,0);
	v[1] = m(1,1);
	v[2] = m(0,1);
}
template<>
void vector_from_matrix<job::Axisymmetric>(const Matrix& m,Vector& v)
{
	v.clear();
	v.resize(job::axisymmetric::Voigt);
	v[0] = m(0,0);
	v[1] = m(1,1);
	v[2] = m(2,2);
	v[3] = m(0,1);
}
template<>
void vector_from_matrix<job::ThreeDimension>(const Matrix& m,Vector& v)
{
	v.clear();
	v.resize(job::three_dimension::Voigt);
	for ( size_type i = 0; i < g_voigt_mapping.size() ; ++ i )
	{
		size_type k = g_voigt_mapping.find(i)->second.first;
		size_type l = g_voigt_mapping.find(i)->second.second;
		v[i] = m(k,l);
	}
}

// 
// Convert vector to matrix using Voigt notation
// 

template<int>
void vector_from_matrix<job::ThreeDimension>(const Vector& v, Matrix& m);


} // end namespace mapping

// this namespace will include major kinematics operations
//
namespace kinematics
{
	using boost::numeric::ublas::trans;
	using boost::numeric::ublas::prod;
	using boost::numeric::ublas::identity_matrix;

	template<int kind,typename MatrixT = boost::numeric::ublas::matrix<value_type> >
	struct invariant;

	template<typename MatrixT>
	struct invariant<1,MatrixT>
	{
		typedef typename MatrixT::value_type value_type;
		static inline value_type get(const MatrixT& matrix) 
		{
			// I1(Q) = E**Q
			return matrix(0,0)+matrix(1,1)+matrix(2,2);
		}
	};
	template<typename MatrixT>
	struct invariant<2,MatrixT>
	{
		typedef typename MatrixT::value_type value_type;
		static inline value_type get(const MatrixT& matrix) 
		{
			// I2(Q) = 0.5 [ I1^2(Q)-I1(Q^2) ]
			value_type i1 = invariant<1>(matrix);
			value_type i2 = invariant<1>(boost::numeric::ublas::prod(matrix,matrix));
			return 0.5*(i1*i1-i2);
		}
	};
	template<typename MatrixT>
	struct invariant<3,MatrixT>
	{
		typedef typename MatrixT::value_type value_type;
		static inline value_type get(const MatrixT& matrix) 
		{
			// I3(Q) = det(Q)
			return det3x3(matrix);
		}
	};
	namespace measures
	{
		namespace right
		{
			template<typename MatrixT>
			inline void cauchy_green(const MatrixT& F, MatrixT& RCG)
			{
				RCG = prod(trans(F),F);
			}
		};
		namespace left
		{

		};
	};
	namespace strains
	{
		namespace right
		{
			template<typename MatrixT>
			inline void cauchy_green(const MatrixT& F, MatrixT& RCG)
			{
				RCG = 0.5*(prod(trans(F),F) - identity_matrix<value_type>(3));
			}
		};
		namespace left
		{

		};
	};
};



template<typename MatrixT>
class constitutive_equation_base
{
public:
	constitutive_equation_base(value_type const1, value_type const2) : const1_(const1), const2_(const2){}
	virtual void stress_cauchy(const MatrixT& F,MatrixT& S) const = 0;
	virtual void stress_cauchy_small(const MatrixT& F,MatrixT& S,job::type type) const = 0;
	virtual void construct_ctensor(Tensor4Rank& c_tensor,const MatrixT& F) const = 0;
	virtual void elasticity_matrix(Matrix& C,job::type type) const = 0;
	

	virtual ~constitutive_equation_base(){}
protected:
	value_type const1_; // material constant1
	value_type const2_; // material constant2
};

template<typename MatrixT>
class model_A5 : public constitutive_equation_base<MatrixT>
{
public:
	typedef constitutive_equation_base<MatrixT> Parent;
public:
	model_A5 (value_type const1, value_type const2) : constitutive_equation_base<MatrixT>(const1,const2){}
	virtual void stress_cauchy(const MatrixT& F,MatrixT& S) const
	{
		assert(F.size1() == F.size2() && F.size1() == 3);
		assert(S.size1() == S.size2() && F.size1() == 3);
		MatrixT C(boost::numeric::ublas::template zero_matrix<value_type>(3,3));
		kinematics::strains::right::cauchy_green(F,C);
#ifdef DETF
		const typename MatrixT::value_type inv_detF = 1.0/det3x3(F);
		S = inv_detF*(Parent::const1_*kinematics::template invariant<1>::get(C)*boost::numeric::ublas::template identity_matrix<value_type>(3)+
			2*Parent::const2_*C);
#else
		S = (Parent::const1_*kinematics::template invariant<1>::get(C)*boost::numeric::ublas::template identity_matrix<value_type>(3)+
			2*Parent::const2_*C);
#endif
	}

	virtual void stress_cauchy_small(const MatrixT& F,MatrixT& S,job::type type) const
	{
		MATRIX(M,type,type);
		VECTOR(V,type);
		MatrixT C(boost::numeric::ublas::template zero_matrix<value_type>(3,3));
		kinematics::strains::right::cauchy_green(F,C);

	}

	virtual void construct_ctensor(Tensor4Rank& c_tensor,const MatrixT& F) const
	{
		//C_{IJKL} = lambda*delta(IJ)*delta(KL) + 2mu*delta(IK)*delta(JL);
		for ( size_type i = 0; i < MAX_DOF; ++ i )
			for ( size_type j = 0; j < MAX_DOF; ++ j )
				for ( size_type k = 0; k < MAX_DOF; ++ k )
					for ( size_type l = 0; l < MAX_DOF; ++ l )
						c_tensor[i][j][k][l] = 0;
		for ( size_type i = 0; i < MAX_DOF; ++ i )
			for ( size_type j = 0; j < MAX_DOF; ++ j )
				for ( size_type k = 0; k < MAX_DOF; ++ k )
					for ( size_type l = 0; l < MAX_DOF; ++ l )
					{
						c_tensor[i][j][k][l] = 0;
						for ( size_type I = 0; I < MAX_DOF; ++ I )
						for ( size_type J = 0; J < MAX_DOF; ++ J )
						for ( size_type K = 0; K < MAX_DOF; ++ K )
						for ( size_type L = 0; L < MAX_DOF; ++ L )
							c_tensor[i][j][k][l] += 
							(Parent::const1_*kroneker_delta(i,j)*kroneker_delta(k,l)+ \
//								2*Parent::const2_*(kroneker_delta(i,k)*kroneker_delta(j,l));
								2*Parent::const2_*(kroneker_delta(i,k)*kroneker_delta(j,l)+kroneker_delta(i,l)*kroneker_delta(j,k)))* \
									F(i,I)*F(j,J)*F(k,K)*F(l,L);
					}
	};
	virtual void elasticity_matrix(Matrix& C,job::type type) const
	{
		construct_elasticity_matrix<job>(C);
	}

protected:
	template<int>
	void construct_elasticity_matrix(Matrix& C) const;
	template<>
	void construct_elasticity_matrix<job::PlaneStrain>(Matrix& C) const
	{
		value_type l = Parent::const1_;
		value_type m = Parent::const1_;
		C.clear();
		C.resize(job::plane_strain::Voigt,job::plane_strain::Voigt);
		C(0,0) = l+2*m;	C(0,1) = l;		
		C(1,0) = l;		C(1,1) = l+2*m;
										C(2,2) = 2*m;
	}
	template<>
	void construct_elasticity_matrix<job::Axisymmetric>(Matrix& C) const
	{
		value_type l = Parent::const1_;
		value_type m = Parent::const1_;	
		C.clear();
		C.resize(job::axisymmetric::Voigt,job::plane_strain::Voigt);
		C(0,0) = l+2*m;	C(0,1) = l;		
		C(1,0) = l;		C(1,1) = l+2*m;
										C(2,2) = 2*m;
														C(3,3) = 2*m;
	}
	template<>
	void construct_elasticity_matrix<job::ThreeDimension>(Matrix& C) const
	{
		C.clear();
		C.resize(job::three_dimension::Voigt,job::plane_strain::Voigt);
	}

};

template<typename MatrixT>
class model_elastic : public constitutive_equation_base<MatrixT>
{
public:
	typedef constitutive_equation_base<MatrixT> Parent;
public:
	model_elastic (value_type E, value_type nu) : constitutive_equation_base<MatrixT>(E, nu){}
	virtual void stress_cauchy(const MatrixT& F,MatrixT& S) const
	{
		assert(F.size1() == F.size2() && F.size1() == 3);
		assert(S.size1() == S.size2() && F.size1() == 3);
		MatrixT C(boost::numeric::ublas::template zero_matrix<value_type>(3,3));
		kinematics::strains::right::cauchy_green(F,C);
		const value_type E = Parent::const1_;
		const value_type nu = Parent::const2_;
		const value_type lambda = E*nu/((1+nu)*(1-2*nu));
		const value_type mu = E/(2*(1+nu));
#ifdef DETF
		const typename MatrixT::value_type inv_detF = 1.0/det3x3(F);
		S = inv_detF*(lambda*kinematics::template invariant<1>::get(C)*boost::numeric::ublas::template identity_matrix<value_type>(3)+
			2*mu*C);
#else
		S = lambda*kinematics::template invariant<1>::get(C)*boost::numeric::ublas::template identity_matrix<value_type>(3)+
			2*mu*C;
#endif
	}
	virtual void construct_ctensor(Tensor4Rank& c_tensor,const MatrixT& F) const
	{
		//C_{IJKL} = lambda*delta(IJ)*delta(KL) + 2mu*delta(IK)*delta(JL);
		const value_type E = Parent::const1_;
		const value_type nu = Parent::const2_;
		const value_type lambda = E*nu/((1+nu)*(1-2*nu));
		const value_type mu = E/(2*(1+nu));
		const value_type K = E/(3*(1-2*nu));
	
		for ( size_type i = 0; i < MAX_DOF; ++ i )
			for ( size_type j = 0; j < MAX_DOF; ++ j )
				for ( size_type k = 0; k < MAX_DOF; ++ k )
					for ( size_type l = 0; l < MAX_DOF; ++ l )
					{
						const value_type Iijkl = 0.5*(kroneker_delta(i,k)*kroneker_delta(j,l)+kroneker_delta(i,l)*kroneker_delta(j,k));
						const value_type delta = kroneker_delta(i,j)*kroneker_delta(k,l);
						c_tensor[i][j][k][l] = 
//							lambda*kroneker_delta(i,j)*kroneker_delta(k,l)+2*mu*kroneker_delta(i,k)*kroneker_delta(j,l);
//									lambda*kroneker_delta(i,j)*kroneker_delta(k,l)+
//									2*mu*(kroneker_delta(i,k)*kroneker_delta(j,l)+kroneker_delta(i,l)*kroneker_delta(j,k));
									K*delta+2*mu*(Iijkl-1/3.*delta);
					}
	};
};



#endif // __MODELS_H__
