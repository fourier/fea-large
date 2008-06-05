#ifndef __CONSOLE_H__
#define __CONSOLE_H__

#include "std.h"
#include "interfaces.h"

class lua_wrapper
{
public:
	lua_wrapper();
	virtual ~lua_wrapper();
	// eval lua commands in 'buffer'
	// returns true if buffer evaluated successfully
	bool eval_buffer(const std::string& buffer);// throw std::logic_error;
	template<typename T> 
	T get_variable(const char*  var_name);

	template<typename T>
	bool set_variable(const char*  var_name,T var_value)
	{
		std::stringstream stream;
		stream << var_name << " = \"" << var_value << "\"";
		return eval_buffer(stream.str());
	}

	std::string last_error();
protected:
	lua_State* lua_state_;
};


class console : public config_interface
{
public:
	console();
	virtual ~console();
	virtual void run(std::istream& input);
protected:
	enum TParse
	{
		EExit,
		ELua,
		EStartLinear,
		EStartNonlinear,
		EHelp,
		EUnknown	
	};
	void out_intro();
	void out_exit();
	void out_error();
	void out_help();
	void out_prompt();
	TParse process_input(const std::string& input);

	bool set_variable(const char* var_name, double var_value);
	bool set_variable_string(const char* var_name,const char* var_value);
	double get_variable(const char* var_name);
	std::string get_variable_string(const char* var_name);
protected:
	std::auto_ptr<lua_wrapper> lua_;	
};


#endif
