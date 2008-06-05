#ifndef __METHODS_H__
#define __METHODS_H__

#include "genericsolver.h"

template<typename ElementT,typename DataLoaderT,typename GaussNodesT>
class method
{
public:
	typedef generic_solver_base<ElementT,DataLoaderT,GaussNodesT> generic_solver_interface;
public:
	method(generic_solver_interface& solver) : solver_(solver){}
protected:
	generic_solver_interface& solver_;

};

#endif // __METHODS_H__
