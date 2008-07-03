#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

// precompiled header
#include "std.h"

template<typename T>
T kroneker_delta(T i, T j)
{
	return i == j? 1: 0;
}

template<typename T>
T factorial(T n)
{
	return n > 0 ? n*factorial(n-1) : 1;
}

template<typename SquareMatrix>
inline void eye(SquareMatrix& matrix)
{
	assert(matrix.size1() == matrix.size2());
	typename SquareMatrix::size_type size = matrix.size1();
	for ( typename SquareMatrix::size_type i = 0; i < size; ++ i )
		matrix(i,i) = 1;
}

template<typename T>
inline void Dump(const T& arg, const char* fname)
{
	FILE* file = fopen(fname,"w+");
	if ( file )
	{
		for ( typename T::size_type i = 0; i < arg.size1(); ++i )
		{
			for ( typename T::size_type j = 0; j < arg.size2(); ++j )
				fprintf(file,"%.8f ",arg(i,j));
			fprintf(file,"\n");
		}
		fclose(file);
	}
}

template<>
inline void Dump<Vector>(const Vector& arg, const char* fname)
{
	FILE* file = fopen(fname,"w+");
	if ( file )
	{
		for ( Vector::size_type i = 0; i < arg.size(); ++i )
			fprintf(file,"%.8f\n",arg[i]);
		fclose(file);
	}
}

// determinant of matrix 3x3
template<typename MatrixT>
inline typename MatrixT::value_type det3x3(const MatrixT& m)
{
	return 
		m(0,0)*(m(1,1)*m(2,2)-m(1,2)*m(2,1)) - 
		m(0,1)*(m(1,0)*m(2,2)-m(1,2)*m(2,0)) + 
		m(0,2)*(m(1,0)*m(2,1)-m(1,1)*m(2,0));
}

// inverse matrix 3x3
template<typename MatrixT>
inline bool inv3x3(const MatrixT& a, MatrixT& inv)
{
	typedef typename MatrixT::value_type value_type;
	// calculate determinant
	value_type det = det3x3(a);
	if ( det == 0 )
		return false;
	// calculate components
	// first row
	inv(0,0) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))/det;
	inv(0,1) = (a(0,2)*a(2,1)-a(0,1)*a(2,2))/det;
	inv(0,2) = (a(0,1)*a(1,2)-a(0,2)*a(1,1))/det;
	// second row
	inv(1,0) = (a(1,2)*a(2,0)-a(1,0)*a(2,2))/det;
	inv(1,1) = (a(0,0)*a(2,2)-a(0,2)*a(2,0))/det;
	inv(1,2) = (a(0,2)*a(1,0)-a(0,0)*a(1,2))/det;
	// third row
	inv(2,0) = (a(1,0)*a(2,1)-a(1,1)*a(2,0))/det;
	inv(2,1) = (a(0,1)*a(2,0)-a(0,0)*a(2,1))/det;
	inv(2,2) = (a(0,0)*a(1,1)-a(0,1)*a(1,0))/det;
	return true;
};

#endif // __FUNCTIONS_H__
