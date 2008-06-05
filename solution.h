#ifndef __SOLUTION_H__
#define __SOLUTION_H__

#include "std.h"

class config_interface;

class solver_command
{
public:
	enum TMethod
	{
		ELinear,
		ENonlinear
	};

	solver_command(config_interface& config,std::ostream& output_stream = std::cout);
	void execute(TMethod method);

protected:
	void fill_default_values();
	void get_values();
	void out_values();
protected: 
	typedef std::map<std::string,double> variables_map;
	typedef std::map<std::string,std::string> string_variables_map;

protected: 
	variables_map variables_;
	string_variables_map string_variables_;
	config_interface& config_;
	std::ostream& output_stream_;
};

#endif // __SOLUTION_H__
