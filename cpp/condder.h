#ifndef __CONDDER_H__
#define __CONDDER_H__

#include "methods.h"

template<typename ElementT,typename DataLoaderT,typename GaussNodesT>
class condition_derivative_method : public method<ElementT,DataLoaderT,GaussNodesT>
{
public:
	typedef method<ElementT,DataLoaderT,GaussNodesT> Parent;
	typedef typename Parent::generic_solver_interface generic_solver_interface;
	typedef typename Parent::generic_solver_interface::ElementsArray ElementsArray;
	typedef typename ElementT::NodeT Node;
public:
	condition_derivative_method(generic_solver_interface& solver) : Parent(solver){}
	void solve(const char* filename,
			   size_type max_loads, 
			   size_type max_newton, 
			   size_type max_linesearch,
			   value_type min_tolerance,
			   value_type max_tolerance);
protected:
	value_type line_search(const ElementsArray& elements,
						   const Vector& solution,
						   size_type max_iter) const;
	value_type average_stress(const ElementsArray& elements, size_type i,size_type j);
	value_type overall_size(size_type dof_number);
protected:
	using Parent::solver_;
};

template<typename E,typename D,typename G>
void condition_derivative_method<E,D,G>::solve(const char* filename,
								 size_type max_loads, 
								 size_type max_newton, 
								 size_type max_linesearch,
								 value_type min_tolerance,
								 value_type max_tolerance)
{
	ElementsArray elements = solver_.elements();
	
	value_type max_displ = solver_.max_prescribed_bc(1);
    value_type current_displ = 0;
	value_type size = overall_size(1);
	value_type step = max_displ/size;

	for ( size_type loading_step = 0; loading_step < max_loads; ++ loading_step )
	{
		value_type tolerance = max_tolerance;
		current_displ += step;
		FILE* fil = fopen("stress22.txt","a");
		solver_.update_with_prescribed_conditions(elements);
		for ( size_type newton_step = 0; newton_step < max_newton; ++ newton_step )
		{
			solver_.prepare_gauss_nodes(elements);
			solver_.update_graddef_array(elements);

			solver_.construct_global_stiffness_matrix(elements);
			solver_.construct_global_residual_force_vector(elements);
			solver_.apply_bc_prescribed(true);
			Vector solution = generic_solver_interface::solve(solver_.global_matrix(),-solver_.global_vector());
			tolerance = fabs(boost::numeric::ublas::inner_prod(solver_.global_vector(),solution));
			std::cout << "newton iteration " << newton_step << " tolerance " << tolerance << std::endl;
			solver_.update_with_solution(elements,solution);
			if ( tolerance < min_tolerance )
				break;
			value_type step = line_search(elements,solution,max_linesearch);
			std::cout << "step = " << step << std::endl;
			if ( step != 0.0 )
				solver_.update_with_solution(elements,step*solution);
		}
		std::cout << "loading iteration " << loading_step << " finished with tolerance " << tolerance <<  std::endl;
		char filename[50];
		sprintf(filename,"iteration%d.msh",loading_step);
		solver_.postprocessing(filename,elements);
		value_type avg_stress = average_stress(elements,1,1);
		fprintf(fil,"%f %f\n",1+current_displ,avg_stress);
		fclose(fil);
		if ( tolerance > max_tolerance )
			break;
	}

	std::cout << "Calculation stopped." << std::endl;
}

template<typename E,typename D,typename G>
value_type condition_derivative_method<E,D,G>::line_search(const ElementsArray& elements,const Vector& solution,size_type max_iter) const
{
	value_type a = 0.0;
	value_type b = 1;
	value_type c = 0.5;

	value_type R0 = 0, R1 = 0, Reta = 0;
	value_type eta = 0;

	const value_type rho = 0.5;

	// 1: R(0)
	ElementsArray elems = elements;
	solver_.update_with_solution(elems,a*solution);
	// update gauss nodes and derivatives
	solver_.prepare_gauss_nodes(elems);
	// update deformation gradients array
	solver_.update_graddef_array(elems);
	// constuct residual vector based on 
	solver_.construct_global_residual_force_vector(elems);
	R0 = (boost::numeric::ublas::inner_prod(solver_.global_vector(),solution));

	// 2: R(1)
	elems = elements;
	solver_.update_with_solution(elems,b*solution);
	// update gauss nodes and derivatives
	solver_.prepare_gauss_nodes(elems);
	// update deformation gradients array
	solver_.update_graddef_array(elems);
	// constuct residual vector based on 
	solver_.construct_global_residual_force_vector(elems);
	R1 = (boost::numeric::ublas::inner_prod(solver_.global_vector(),solution));


	Reta = R1;
	
	typedef std::pair<value_type,value_type> pair;
	std::list<pair> pair_list;
	pair_list.push_back(pair(0,fabs(R0)));

	for ( size_type i = 0; i < max_iter; ++ i )
	{
		// calculate eta
		value_type alpha = R0/Reta;
		if ( alpha < 0 )
			eta = alpha/2. + sqrt(pow(alpha/2.,2)-alpha);
		else
			eta = alpha/2.;

		// test eta
		elems = elements;
		solver_.update_with_solution(elems,eta*solution);
		// update gauss nodes and derivatives
		solver_.prepare_gauss_nodes(elems);
		// update deformation gradients array
		solver_.update_graddef_array(elems);
		// constuct residual vector based on 
		solver_.construct_global_residual_force_vector(elems);

		Reta = (boost::numeric::ublas::inner_prod(solver_.global_vector(),solution));
		pair_list.push_back(pair(eta,fabs(Reta)));

		if ( fabs(Reta) < rho*fabs(R0) )
			break;
	}
	pair_list.push_back(pair(1,fabs(R1)));
	std::list<pair>::const_iterator it = pair_list.begin();

	for ( ; it != pair_list.end(); ++ it )
	{
		if ( it->second < fabs(Reta) )
		{
			Reta = it->second;
			eta = it->first;
		}
	}

	return eta;
}

template<typename E,typename D,typename G>
value_type condition_derivative_method<E,D,G>::average_stress(const ElementsArray& elements, size_type i,size_type j)
{
	value_type avg_stress = 0;
	for ( size_type el = 0; el < elements.size(); ++ el )
	{
		Node nod = solver_.elements()[el].center();
		MATRIX(F,3,3);
		solver_.deformation_gradient(el,elements,nod,F);
		MATRIX(S,3,3);
		solver_.model().stress_cauchy(F,S);
		avg_stress += S(i,j);
	}	
	avg_stress /= elements.size();
	return avg_stress;
}

template<typename E,typename D,typename G>
value_type condition_derivative_method<E,D,G>::overall_size(size_type dof_number)
{
	value_type min = solver_.elements()[0].node(0).dof[dof_number];
	value_type max = solver_.elements()[0].node(0).dof[dof_number];

	for ( size_type el = 1; el < solver_.elements().size(); ++ el )
	{
		for ( size_type i = 0; i < E::NodesNumber; ++ i )
		{
			value_type value = solver_.elements()[el].node(i).dof[dof_number];
			if ( value < min )
				min = value;
			if ( value > max )
				max = value;
		}
	}
	return max - min;
}

#endif // __CONDDER_H__
