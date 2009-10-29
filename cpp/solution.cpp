#include "solution.h"
#include "interfaces.h"
#include "models.h"
#include "linear.h"
#include "axisymmetric6p3.h"
#include "plane6p3.h"
#include "condder.h"

void solver_command::fill_default_values()
{
	variables_.clear();
	string_variables_.clear();
	variables_["loads"] = 100;
	variables_["newton"] = 50;
	variables_["tolerance"] = 1e-7;
	variables_["searches"] = 20;
	variables_["max_tolerance"] = 10;
	variables_["mu"] = 100;
	variables_["lambda"] = 100;
	variables_["nu"] = 0.3;
	variables_["E"] = 1E9;
	string_variables_["model"] = "elastic";
	string_variables_["input"] = "matlab_geometry.node";
	string_variables_["output"] = "solution.msh";
}


solver_command::solver_command(config_interface& config,std::ostream& output_stream) : 
	config_(config), output_stream_(output_stream)
{
	fill_default_values();
	variables_map::const_iterator it = variables_.begin();
	for ( ; it != variables_.end(); ++ it )
		config.set_variable(it->first.c_str(), it->second);
	string_variables_map::const_iterator it1 = string_variables_.begin();
	for ( ; it1 != string_variables_.end(); ++ it1 )
		config.set_variable_string(it1->first.c_str(), it1->second.c_str());
}

void solver_command::get_values()
{
	variables_map::iterator it = variables_.begin();
	for ( ; it != variables_.end(); ++ it )
	{
		double value = config_.get_variable(it->first.c_str());
		it->second = value;
		
	}
	string_variables_map::iterator it1 = string_variables_.begin();
	for ( ; it1 != string_variables_.end(); ++ it1 )
	{
		std::string value = config_.get_variable_string(it1->first.c_str());
		it1->second = value;
	}
}

void solver_command::out_values()
{
	output_stream_ << "Calculation parameters: " << std::endl;
	
	variables_map::const_iterator it = variables_.begin();
	string_variables_map::const_iterator it1 = string_variables_.begin();
	
	for ( ; it1 != string_variables_.end(); ++ it1 )
		output_stream_ << it1->first.c_str() << " is " << it1->second.c_str() << std::endl;
	for ( ; it != variables_.end(); ++ it )
		output_stream_ << it->first.c_str() << " is " << it->second << std::endl;

}

void solver_command::execute(solver_command::TMethod method)
{
//	axisymmetric6p3 solver(new model_elastic<Matrix>(1e+9,0.4));

	get_values();

	size_type loads = (size_type)variables_["loads"];
	size_type newton = (size_type)variables_["newton"];
	value_type tolerance = variables_["tolerance"];
	size_type searches = (size_type)variables_["searches"];
	value_type max_tolerance = variables_["max_tolerance"];
	value_type mu = variables_["mu"] ;
	value_type lambda = variables_["lambda"];
	value_type nu = variables_["nu"];
	value_type E = variables_["E"];
	std::string input = string_variables_["input"];
	std::string output = string_variables_["output"];
	std::string model_type = string_variables_["model"];

	if ( input.length() == 0 )
	{
		output_stream_ << "Error: 'input' variable not set. Cannot continue" << std::endl;
		return;
	}

//	axisymmetric6p3 task(new model_A5<Matrix>(alpha,mu));
//	typedef generic_solver<plane::triangle::element6,LGeometryData,triangle::gauss_tri7p> solver;
	
	constitutive_equation_base<Matrix>* model = NULL;
	if ( model_type == "A5" )
		model = new model_A5<Matrix>(lambda,mu);
	else if (model_type == "elastic" )
		model = new model_elastic<Matrix>(E,nu);
	else 
	{
		std::cout << "Material model \"" << model_type << "\" not supported" << std::endl; 
		return;
	}
	plane6p3 task(model);
	
	if ( task.load_file(input.c_str()) )
	{
		out_values();
		switch(method)
		{
		case ELinear:
			{
				output_stream_ << "Starting linear solution process..." << std::endl;

//				linear_method<axisymmetric6p3::Element,axisymmetric6p3::DataLoader,axisymmetric6p3::GaussNodes> 
//					linear(task);

				linear_method<plane6p3::Element,plane6p3::DataLoader,plane6p3::GaussNodes> 
					linear(task);
				linear.solve(output.c_str());
	
			}
			break;
		case ENonlinear:
			{

				output_stream_ << "Starting nonlinear solution process..." << std::endl;
				// 1st: parse output path to get filenames w/o extension
				// boost::filesystem::path path(output.c_str());

//				condition_derivative_method<axisymmetric6p3::Element,axisymmetric6p3::DataLoader,
//					axisymmetric6p3::GaussNodes> 
//					connder(task);

				condition_derivative_method<plane6p3::Element,plane6p3::DataLoader,plane6p3::GaussNodes> 
					connder(task);
				connder.solve(output.c_str(),loads,newton,searches,tolerance,max_tolerance);
			}
			break;
		}
	}
	else
	{
		output_stream_ << "Error: Could not load " << input.c_str() << ": Cannot continue" << std::endl;
	}
}



