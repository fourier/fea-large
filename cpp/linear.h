#ifndef __LINEAR_H__
#define __LINEAR_H__

#include "methods.h"

template<typename ElementT,typename DataLoaderT,typename GaussNodesT>
class linear_method : public method<ElementT,DataLoaderT,GaussNodesT>
{
public:
	typedef method<ElementT,DataLoaderT,GaussNodesT> Parent;
	typedef typename Parent::generic_solver_interface generic_solver_interface;
	typedef typename Parent::generic_solver_interface::ElementsArray ElementsArray;
public:
	linear_method(generic_solver_interface& solver) : Parent(solver){}
	void solve(const char* filename);
protected:
	using Parent::solver_;
};

template<typename E,typename D,typename G>
void linear_method<E,D,G>::solve(const char* filename)
{
	ElementsArray elements = solver_.elements();
	solver_.linear_construct_global_matrix(elements);
	solver_.apply_bc_prescribed();
	Dump(solver_.global_matrix(),"Dump/global_applied.txt");
	Dump(solver_.global_vector(),"Dump/global_vector.txt");
	Vector global_solution = generic_solver_interface::solve(solver_.global_matrix(),solver_.global_vector());
	Dump(global_solution,"solution.txt");
	solver_.update_with_solution(elements,global_solution );
	solver_.postprocessing(filename,elements);
}


#endif // __LINEAR_H__
