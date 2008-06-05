#ifndef __LSLAE_H__
#define __LSLAE_H__

// precompiled headers
#include "std.h"

#pragma warning( disable : 4290 )

//using namespace boost::numeric::ublas;

#define ROUND_LESS_THAN_1E_MINUS_16

namespace LSLAE
{
	template<class MatrixT, class VectorT>
	VectorT gauss_solve_swap(const MatrixT& a, const VectorT& b) throw (std::logic_error)
	{
		typedef typename MatrixT::size_type size_type;
		typedef typename MatrixT::value_type value_type;
		if ( a.size1() != b.size() )
			throw std::logic_error("Wrong sizes");
		if ( a.size1() != a.size2() )
			throw std::logic_error("Not square matrix");

		size_type n = a.size2();

		VectorT result( b );

		MatrixT mtx( a );

		size_type i,j,k;
		
		value_type div=mtx( 0, 0 );

		size_type big = 0;

		for ( i = 0 ; i < n; ++ i )
		{
			big = i;
			for ( j = i; j < n; ++ j )
			{
				if ( fabs( mtx( j, i ) ) > fabs( mtx( big, i ) ) )
					big = j;
			}

			if ( big != i )
			{
/*				std::swap_ranges( boost::numeric::ublas::row( mtx, big ).begin(),
					boost::numeric::ublas::row( mtx, big ).end(),
					boost::numeric::ublas::row( mtx, i ).begin() );*/
					std::swap_ranges( 
							boost::numeric::ublas::matrix_row<MatrixT>( mtx, big ).begin(),
					boost::numeric::ublas::matrix_row<MatrixT>( mtx, big ).end(),
					boost::numeric::ublas::matrix_row<MatrixT>( mtx, i ).begin() );

				std::swap( result[i], result[big] );
			}

			if ( mtx ( i, i ) == 0.0 ) 
			{
				throw std::logic_error("Singular matrix");
			}

			div = mtx( i, i );
			
			mtx( i, i ) = 1;
			
			for ( j = i + 1; j < n ; ++ j )
				mtx( i, j ) = mtx( i, j ) / div;

			result(i) /= div;
			
			for ( k = i + 1; k < n ; ++ k )
			{ 
				if ( mtx( k, i ) != 0.0 )
				{
					div = mtx( k , i );
					for ( j = i; j < n; ++ j )
					{
						mtx( k, j ) = mtx( k, j ) - mtx( i, j ) * div;
					}
					result[k] = result[k] - result[i]*div;
				}
			}
		}
#ifdef DO_NOT_USE_BOOST_SOLVE
		for ( i = n - 1; i >= 0; -- i )
		{
			for ( k = i - 1; k >= 0; -- k )
			{
				if ( mtx( k, i ) != 0.0 )
				{
					div = mtx( k, i );

					for ( j = n - 1; j >= i; -- j )
					{
						mtx( k, j ) = mtx( k, j ) - mtx( i, j )*div;
					}
					result[k] = result[k] - result[i] * div; 
				}
			}
		}
#else
		boost::numeric::ublas::inplace_solve(mtx,
			result,
			boost::numeric::ublas::upper_tag());
#endif

#ifdef ROUND_LESS_THAN_1E_MINUS_16
		typename VectorT::iterator it = result.begin();
		for ( ;
			it != result.end(); 
			++ it )
		{
			if ( fabs(*it) < 1e-15 )
				*it = 0.0;
		}
#endif
		return result;
	}

	template<class MatrixT, class VectorT>
	VectorT gauss_solve(const MatrixT& a, const VectorT& b) throw (std::logic_error)
	{
		BOOST_STATIC_ASSERT( sizeof(MatrixT::size_type) == sizeof(VectorT::size_type) );
		typedef typename MatrixT::size_type size_type;
		typedef typename MatrixT::value_type value_type;

		if ( a.size1() != b.size() )
			throw std::logic_error("Wrong sizes");
		if ( a.size1() != a.size2() )
			throw std::logic_error("Not square matrix");


		size_type i,j,k;
		size_type n = a.size2();

		VectorT result( b );

		MatrixT mtx( a );

		value_type div=mtx( 0, 0 );

		size_type big = 0;

		std::vector<size_type> indexes(n);
		typename std::vector<size_type>::reverse_iterator rit,rit1;
		typename std::vector<size_type>::iterator it,it1;
		
		i = 0;

		for ( it = indexes.begin(); it != indexes.end(); ++ it )
			*it = i++ ;

		i = 0;

		for ( it = indexes.begin(); it != indexes.end(); ++ it )
		{
			big = *it;
			
			for ( it1 = it; it1 != indexes.end(); ++ it1 )
			{
				if ( fabs( mtx( *it1, i ) ) > fabs( mtx( big, i ) ) )
					big = *it1;
			}

			if ( big != *it )
				std::swap( indexes[i], indexes[big] );

			if ( mtx ( *it, i ) == 0.0 ) 
			{
				throw std::logic_error("Singular matrix");
			}

			div = mtx( *it, i );
			
			mtx( *it, i ) = 1;
			
			for ( j = i + 1; j < n ; ++ j )
				mtx( *it, j ) = mtx( *it, j ) / div;
			
			result(*it) /= div;
			
			for ( it1 = boost::next(it); it1 != indexes.end(); ++ it1 )
			{ 
				k = *it1;
			
				if ( mtx( k, i ) != 0 )
				{
					div = mtx( k , i );
					for ( j = i; j < n; ++ j )
					{
						mtx( k, j ) = mtx( k, j ) - mtx( *it, j ) * div;
					}
					result[k] = result[k] - result[*it]*div;
				}
			}
			i++;
		}

		i = n - 1;

		for ( rit = indexes.rbegin(); rit != indexes.rend(); ++ rit )
		{
			for ( rit1 = boost::next(rit); rit1 != indexes.rend(); ++ rit1 )
			{
				k = *rit1;
				if ( mtx( k, i ) != 0.0 )
				{
					div = mtx( k, i );

					for ( j = n - 1; j >= i; -- j )
					{
						mtx( k, j ) = mtx( k, j ) - mtx( *rit, j )*div;
					}
					result[k] = result[k] - result[*rit] * div; 
				}
			}
			i--;
		}

#ifdef ROUND_LESS_THAN_1E_MINUS_16
		for ( typename VectorT::iterator it = result.begin();
			it != result.end(); 
			++ it )
		{
			if ( fabs(*it) < 1e-15 )
				*it = 0.0;
		}
#endif
		VectorT solution(n);
		for ( i = 0; i < n; ++ i )
			solution[i] = result[indexes[i]];

		return solution;

	}
	
	template<class MatrixT, class VectorT>
	VectorT gauss_solve_ex(const MatrixT& a, const VectorT& b) throw (std::logic_error)
	{
		BOOST_STATIC_ASSERT( sizeof(typename MatrixT::size_type) == sizeof(typename VectorT::size_type) );
	
		typedef typename MatrixT::value_type value_type;
		typedef typename MatrixT::size_type size_type;


		if ( a.size1() != b.size() )
			throw std::logic_error("Wrong sizes");
		if ( a.size1() != a.size2() )
			throw std::logic_error("Not square matrix");

		size_type n = a.size2();

		VectorT result( b );

		MatrixT mtx( a );

		value_type div=mtx( size_type(), size_type() );

		std::vector<size_type> indexes_rows(n);
		std::vector<size_type> indexes_cols(n);
		typename std::vector<size_type>::reverse_iterator rit,ritc,rit1,rit2;
		typename std::vector<size_type>::iterator it,itc,it1,it2;
		
		const value_type null = value_type();

		size_type i = size_type();

		for ( it = indexes_rows.begin(),it1 = indexes_cols.begin();
			  it != indexes_rows.end() && it1 != indexes_cols.end(); ++ it, ++ it1, ++ i )
			  *it = i,*it1 = i;
		
		for ( it = indexes_rows.begin() , itc = indexes_cols.begin(); it != indexes_rows.end(); ++ it, ++ itc )
		{
			typename std::vector<size_type>::iterator big_row = it,big_col = itc;

			for ( it1 = it; it1 != indexes_rows.end(); ++ it1 )
			{
				if ( fabs( mtx( *it1, *itc ) ) > fabs( mtx( *big_row, *itc ) ) )
					big_row = it1;
			}
			for ( it1 = itc; it1 != indexes_cols.end(); ++ it1 )
			{
				if ( fabs( mtx( *it, *it1 ) ) > fabs( mtx( *it, *big_col ) ) )
					big_col = it1;
			}

			if ( fabs( mtx( *big_row, *itc ) ) > fabs( mtx( *it, *big_col ) ) )
			{
				if ( big_row != it )
					std::iter_swap( it, big_row );
			}
			else
			{
				if ( big_col != itc )
					std::iter_swap( itc, big_col );
			}

			if ( mtx ( *it, *itc ) == null ) 
			{
				std::cout << "Row " << *it << "is empty" << std::endl;
				throw std::logic_error("Singular matrix");
			}

			div = mtx( *it, *itc );
			
			mtx( *it, *itc ) = 1;
			
			for ( it1 = boost::next(itc); it1 != indexes_cols.end(); ++ it1 )
				mtx( *it, *it1 ) = mtx( *it, *it1 ) / div;
			
			result(*it) /= div;
			
			for ( it1 = boost::next(it); it1 != indexes_rows.end(); ++ it1 )
			{ 
				if ( mtx( *it1, *itc ) != null )
				{
					div = mtx( *it1 , *itc );
					for ( it2 = itc; it2 != indexes_cols.end(); ++ it2 )
					{
						mtx( *it1, *it2 ) = mtx( *it1, *it2 ) - mtx( *it, *it2 ) * div;
					}
					result[*it1] = result[*it1] - result[*it]*div;
				}
			}
		}

		for ( rit = indexes_rows.rbegin(), ritc = indexes_cols.rbegin(); rit != indexes_rows.rend(); ++ rit, ++ ritc )
		{
			for ( rit1 = boost::next(rit); rit1 != indexes_rows.rend(); ++ rit1 )
			{
				if ( mtx( *rit1, *ritc ) != null )
				{
					div = mtx( *rit1, *ritc );

					for ( rit2 = indexes_cols.rbegin(); rit2 != boost::next(ritc)/*rit2 >= i*/; ++ rit2 )
					{
						mtx( *rit1, *rit2 ) = mtx( *rit1, *rit2 ) - mtx( *rit, *rit2 )*div;
					}
					result[*rit1] = result[*rit1] - result[*rit] * div; 
				}
			}
		}

#ifdef ROUND_LESS_THAN_1E_MINUS_16
		for ( typename VectorT::iterator it = result.begin();
			it != result.end(); 
			++ it )
		{
			if ( fabs(*it) < 1e-15 )
				*it = null;
		}
#endif
		VectorT solution(n);
		for ( i = 0; i < n; ++ i )
			solution[i] = result[indexes_rows[i]];
		for ( i = 0; i < n; ++ i )
			result[i] = solution[indexes_cols[i]];

		return result;
	}

	template<class MatrixT, class VectorT>
	VectorT lu_substitute_solve(const MatrixT& a, const VectorT& b) throw (std::logic_error)
	{
		BOOST_STATIC_ASSERT( sizeof(typename MatrixT::size_type) == sizeof(typename VectorT::size_type) );
	
		typedef typename MatrixT::value_type value_type;
		typedef typename MatrixT::size_type size_type;


		if ( a.size1() != b.size() )
			throw std::logic_error("Wrong sizes");
		if ( a.size1() != a.size2() )
			throw std::logic_error("Not square matrix");

		size_type n = a.size2();

		VectorT result(n);

		MatrixT m(a);
		Matrix m_inv(n,n);
		m_inv.clear();

		boost::numeric::ublas::permutation_matrix<size_type> pm(n);
		boost::numeric::ublas::lu_factorize(m,pm);
		boost::numeric::ublas::lu_substitute(m, pm, m_inv);
		// compute w = A^-1*b
		result = boost::numeric::ublas::prod(m_inv,b);
		return result;
	}
/*
	template<class MatrixT, class VectorT>
	VectorT pcg_solve(const MatrixT& A, const MatrixT& M, const VectorT& b) throw (std::logic_error)
	{
		BOOST_STATIC_ASSERT( sizeof(typename MatrixT::size_type) == sizeof(typename VectorT::size_type) );
	
		typedef typename MatrixT::value_type value_type;
		typedef typename MatrixT::size_type size_type;

		if ( A.size1() != b.size() )
			throw std::logic_error("Wrong sizes");
		if ( A.size1() != A.size2() )
			throw std::logic_error("Not square matrix");

		size_type n = A.size2();

		VectorT result(n);

		MatrixT m(A);
	}
*/	
	template<class MatrixT, class VectorT>
	VectorT cg_solve(const MatrixT& A, const VectorT& b, typename MatrixT::size_type max_iter,typename MatrixT::value_type tol) throw (std::logic_error)
	{
		using boost::numeric::ublas::prod;
		using boost::numeric::ublas::norm_2;
		using boost::numeric::ublas::inner_prod;
	
		BOOST_STATIC_ASSERT( sizeof(typename MatrixT::size_type) == sizeof(typename VectorT::size_type) );
	
		typedef typename MatrixT::value_type value_type;
		typedef typename MatrixT::size_type size_type;

		if ( A.size1() != b.size() )
			throw std::logic_error("Wrong sizes");
		if ( A.size1() != A.size2() )
			throw std::logic_error("Not square matrix");

		size_type n = A.size2();

		VectorT x(n);
		x.clear();

		size_type iter;
		value_type alpha, beta, residn;
		
		VectorT resid = b - prod(A,x); // residual
		VectorT d = resid;             // search direction
		VectorT resid_old;
		VectorT temp;
  
		// CG loop
		for ( iter = 0; iter < max_iter; iter ++ )
		{
			temp = prod(A,d);
			alpha = inner_prod(resid,resid)/inner_prod(d,temp);
			x += (d*alpha);
			resid_old = resid;
			resid -= (temp*alpha);
			residn = norm_2(resid);

			if ( residn < tol )
				break;

			beta = inner_prod(resid,resid)/inner_prod(resid_old,resid_old);
			d = resid + d*beta;
		}
		return x;	
	}

};


#endif
