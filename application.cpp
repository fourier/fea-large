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
		std::list<std::string>::const_iterator it = filenames_.begin();
		for ( ; it != filenames_.end(); ++ it )
		{
			std::ifstream file(it->c_str());
			if ( file )
			{
				console con;
				con.run(file,true);
				file.close();
			}
		}
	}
	else // run interactive shell
	{
		console con;
		con.run(std::cin,false);
	}
}


void application::parse_cmdline(int argc,const char* argv[])
{

	if ( argc > 1 ) // do batch processing
	{
		do_batch_ = true;
		for ( int i = 1; i < argc; ++ i )
			filenames_.push_back(argv[i]);
	}
}
