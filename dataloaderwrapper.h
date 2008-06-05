#ifndef __DATALOADERWRAPPER_H__
#define __DATALOADERWRAPPER_H__

// precompiled headers
#include "std.h"

// local includes

struct index_matrix_interface
{
	virtual size_type& operator () (size_type i, size_type j) = 0;
};

template<typename DataLoaderT,typename ElementT, typename ElementsArrayT = std::vector<ElementT> >
class data_loader_wrapper
{ 
public:
	enum TPrescribedBoundaryType
	{
		EFree = 0, // free;
		EPrescribedX = 1, // : x prescribed (axisymmetric case: r, same the symmetry condition)
		EPrescribedY = 2, // : y prescribed (axisymmetric case: z)
		EPrescribedXY = 3, // : x, y prescribed (axisymmetric case: r,z)
		EPrescribedZ = 4, // : z prescribed
		EPrescribedXZ = 5, // : x, z prescribed
		EPrescribedYZ = 6, // : y, z prescribed
		EPrescribedXYZ = 7  // : x, y, z prescribed.
	};
public:
	typedef DataLoaderT DataLoader;
	typedef ElementT Element;
	typedef ElementsArrayT ElementsArray;
	typedef std::auto_ptr<DataLoader> DataPointer;
	typedef typename Element::NodeT Node;
	typedef std::vector<Node> NodesVector;
	typedef boost::tuple<size_type,size_type,Node> BoundaryItem;
	typedef std::vector<BoundaryItem> PrescribedBoundaryArray;
public:
	// constructor taking ownership
	data_loader_wrapper() : data_(NULL), ptr_array_(NULL)
	{
	};
	
	void reset(DataLoaderT* data)
	{
		data_.reset(data);
	}

	bool load_file(const char* filename)
	{
		std::ifstream f(filename);
		if (f.is_open())
		{
			f >> *data_;
			f.close();
			prepare_boundary();
			prepare_nodes();
			ptr_array_ = NULL;
			return true;
		}
		return false;
	}

	// empty and fill elements array with data from DataLoader
	virtual void load_elements(ElementsArray& array) = 0;

	virtual void prepare_index(index_matrix_interface& m) = 0;

	const PrescribedBoundaryArray& prescribed_boundary_array() const { return prescribed_boundary_array_; }
	const NodesVector& nodes_vector() const { return nodes_vector_;}
	NodesVector& nodes_vector() { return nodes_vector_;}
	
protected:
	virtual void prepare_boundary() = 0;
	virtual void prepare_nodes() = 0;

protected:
	// pointer to a data loader, holding ownership
	DataPointer data_;
	// pointer to elements array, NOT holding ownership
	ElementsArray* ptr_array_;
	PrescribedBoundaryArray prescribed_boundary_array_;
	NodesVector nodes_vector_;
};


#endif // __DATALOADERWRAPPER_H__
