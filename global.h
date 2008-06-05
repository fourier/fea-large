#ifndef __GLOBAL_H__
#define __GLOBAL_H__

// type definitions
// boost typedefs

typedef double value_type;
typedef boost::numeric::ublas::matrix<value_type> Matrix;
typedef boost::numeric::ublas::vector<value_type> Vector;
typedef Matrix::size_type size_type;
typedef boost::numeric::ublas::matrix<size_type> IndexMatrix;
typedef boost::multi_array<value_type, 4> Tensor4Rank;

#undef MATRIX
#define MATRIX(name,rows,cols) Matrix name(rows,cols);name.clear();
#undef VECTOR
#define VECTOR(name,rows) Vector name(rows);name.clear();

const value_type PI = 3.14159265358979323846;
const size_type MAX_DOF = 3;

#define DETF

#endif // __GLOBAL_H__
