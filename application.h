#ifndef __APPLICATION_H__
#define __APPLICATION_H__

// precompiled headers
#include "std.h"

class application
{
public:
	application();
	virtual ~application();
	void run(int argc,const char* argv[]);
protected:
	void parse_cmdline(int argc,const char* argv[]);
protected:
	std::string filename_;	
	bool do_batch_;
};

#endif // __APPLICATION_H__

