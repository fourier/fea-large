#ifndef __MODELS_H__
#define __MODELS_H__

// precompiled headers
#include "std.h"

// local includes
#include "functions.h"

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
	virtual void update_ctensor(Tensor4Rank& c_tensor,const MatrixT& F) const = 0;
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
	virtual void update_ctensor(Tensor4Rank& c_tensor,const MatrixT& F) const
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
	virtual void update_ctensor(Tensor4Rank& c_tensor,const MatrixT& F) const
	{
		//C_{IJKL} = lambda*delta(IJ)*delta(KL) + 2mu*delta(IK)*delta(JL);
		const value_type E = Parent::const1_;
		const value_type nu = Parent::const2_;
		const value_type lambda = E*nu/((1+nu)*(1-2*nu));
		const value_type mu = E/(2*(1+nu));
	
		for ( size_type i = 0; i < MAX_DOF; ++ i )
			for ( size_type j = 0; j < MAX_DOF; ++ j )
				for ( size_type k = 0; k < MAX_DOF; ++ k )
					for ( size_type l = 0; l < MAX_DOF; ++ l )
						c_tensor[i][j][k][l] = 
							lambda*kroneker_delta(i,j)*kroneker_delta(k,l)+2*mu*kroneker_delta(i,k)*kroneker_delta(j,l);
	};
};



#endif // __MODELS_H__
