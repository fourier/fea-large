#include "application.h"
#include "console.h"
#include "elements.h"

application::application() : do_batch_(false)
{
}

application::~application()
{
}


void application::run(int argc,const char* argv[])
{
	test_element_plane_triangle_p6();
	parse_cmdline(argc,argv);
	if ( do_batch_ ) // process files
	{
	}
	else // run interactive shell
	{
		console con;
		con.run(std::cin);
	}
}


void application::parse_cmdline(int argc,const char* argv[])
{
/*
	if ( argc != 1 ) // do batch processing
	{
		filename_ = argv[1];
	}
	else // call console
	{

	}
*/
}
