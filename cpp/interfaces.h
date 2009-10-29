#ifndef __INTERFACES_H__
#define __INTERFACES_H__

#include "std.h"

class config_interface
{
public:
	virtual	bool set_variable(const char* var_name,double var_value) = 0;
	virtual	bool set_variable_string(const char* var_name,const char* var_value) = 0;
	virtual	double get_variable(const char* var_name) = 0;
	virtual	std::string get_variable_string(const char* var_name) = 0;

	virtual ~config_interface(){}
};

#endif // __INTERFACES_H__
