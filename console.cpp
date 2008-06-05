#include "console.h"
#include "solution.h"

template<>
double lua_wrapper::get_variable<double>(const char*  var_name)
{
	lua_getglobal(lua_state_, var_name);
	double value = lua_tonumber(lua_state_, -1);
	return value;
}

template<>
std::string lua_wrapper::get_variable<std::string>(const char*  var_name)
{
	lua_getglobal(lua_state_, var_name);
	std::string value = lua_tostring(lua_state_, -1);
	return value;
}


lua_wrapper::lua_wrapper()
{
	lua_state_ = lua_open();
	assert(lua_state_);
	luaopen_base(lua_state_); // base
#ifdef ALL_LUA_LIBRARIES
	luaopen_table(L); // table
	luaopen_io(L); // I/O
	luaopen_string(L); // string
	luaopen_math(L); // math
#endif
}

lua_wrapper::~lua_wrapper()
{
	lua_close(lua_state_);
}

bool lua_wrapper::eval_buffer(const std::string& buffer)// throw std::logic_error
{
	int error = luaL_loadbuffer(lua_state_, buffer.c_str(), buffer.length(), "line") ||
		lua_pcall(lua_state_, 0, 0, 0);
		return error == 0;
}

std::string lua_wrapper::last_error()
{
	std::string error;
	error	= lua_tostring(lua_state_, -1);
	lua_pop(lua_state_, 1);
	if ( error ==  "" )
		error = "An error occured while processing statement";
	return error;
}

console::console() : lua_(new lua_wrapper)
{
}

console::~console()
{

}


void console::run(std::istream& input)
{
	// read here: http://www.ibm.com/developerworks/ru/library/l-lua/
	// and here: http://www.itc.ua/node/20951/
	// and here: http://lua-users.org/wiki/TutorialDirectory

	bool do_exit = false;
	out_intro();
	solver_command command(*this);
	std::string line;
	do
	{
		out_prompt();
		std::getline(input,line);
		switch (process_input(line))
		{
			case EExit:
				do_exit = true;
				break;
			case EHelp:
				out_help();
				break;
			case EStartLinear:
				command.execute(solver_command::ELinear);
				break;
			case EStartNonlinear:
				command.execute(solver_command::ENonlinear);
				break;
			case ELua:
				if ( !lua_->eval_buffer(line) )
					out_error();
				break;
		}
	} while ( !do_exit );
	out_exit();
}

void console::out_intro()
{
	std::cout << std::endl << "FEA finite strain program" << std::endl;
	std::cout << "=========================" << std::endl;
	std::cout << std::endl << "Interactive session" << std::endl;
	std::cout << std::endl << "Try \"help\" to display short help" << std::endl;
	std::cout << std::endl;
}

void console::out_exit()
{
	std::cout << std::endl << "Thank you!" << std::endl;
}

void console::out_error()
{
	std::cout << std::endl;
	std::cout << lua_->last_error() << std::endl << std::endl;
}


void console::out_help()
{	
	std::cout << "Short help:" << std::endl;
	std::cout << "===========" << std::endl;
	std::cout << "Calculation control paramers could be entered in a following way:" << std::endl;
	std::cout << "variable = value" << std::endl;
	std::cout << std::endl << "Meaningfull parameters" << std::endl;	
	std::cout << "'input' - contains file with input geometry" << std::endl;
	std::cout << "'output' - contains file for output, or folder";
	std::cout << " in case of nonlinear solution" << std::endl;
	std::cout << "'loads' - contains number of load iterations";
	std::cout << " for nonlinear solution" << std::endl;
	std::cout << "'tolerance' - contains tolerance for nonlinear solution" << std::endl;
	std::cout << "'newton' - contains a maximum number of Newton iterations" << std::endl;
	std::cout << "'searches' - maximum number of line searches" << std::endl;
	std::cout << "'alpha' - material constant Alpha(default 100)" << std::endl;
	std::cout << "'mu' - material constant Mu(default 100)" << std::endl;
	std::cout << "Commands:" << std::endl;
	std::cout << "help - display this help" << std::endl;
	std::cout << "exit(quit) - exit program" << std::endl;
	std::cout << "\"start linear\" - start linear solution" << std::endl;
	std::cout << "\"start nonlinear\" - start solution with finite strains" << std::endl;
	// TODO: add more commands and descriptions
	std::cout << "All other commands are commands on a Lua language" << std::endl;
	std::cout << std::endl;
}

void console::out_prompt()
{
	std::cout << "FEA: ";
}

console::TParse console::process_input(const std::string& input)
{
	if ( input == "help" )
		return EHelp;
	if ( input == "exit" || input == "quit" )
		return EExit;
	if ( input ==  "start linear" ) 
		return EStartLinear;
	if ( input ==  "start nonlinear" )
		return EStartNonlinear;
	return ELua;
}

bool console::set_variable(const char* var_name,double var_value)
{
	bool result = true;
	if ( !lua_->set_variable(var_name,var_value) )
	{
		out_error();
		result = false;
	}
	return result;
}

bool console::set_variable_string(const char* var_name,const char* var_value)
{
	bool result = true;
	if ( !lua_->set_variable(var_name,var_value) )
	{
		out_error();
		result = false;
	}
	return result;
}

double console::get_variable(const char* var_name)
{
	return lua_->get_variable<double>(var_name);
}

std::string console::get_variable_string(const char* var_name)
{
	return lua_->get_variable<std::string>(var_name);
}
