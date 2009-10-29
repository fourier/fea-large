#include "application.h"
#include "console.h"
#include "elements.h"
#include "global.h"

VoigtMapping g_voigt_mapping;

application::application() : do_batch_(false)
{
}

application::~application()
{
}

void application::initialize_globals()
{
	// construct Voigt mapping
	// Mapping according to Voigt notation
	// 
	// 11->1 22->2 33->3
	// 23->4 13->5 12->6
	//
	// or, with zero-based index:
	//
	// 00->0 11->1 22->2
	// 12->3 02->4 01->5

	g_voigt_mapping[0] = std::pair<size_type,size_type>(0,0);
	g_voigt_mapping[1] = std::pair<size_type,size_type>(1,1);
	g_voigt_mapping[2] = std::pair<size_type,size_type>(2,2);

	g_voigt_mapping[3] = std::pair<size_type,size_type>(1,2);
	g_voigt_mapping[4] = std::pair<size_type,size_type>(0,2);
	g_voigt_mapping[5] = std::pair<size_type,size_type>(0,1);
}

void application::run(int argc,const char* argv[])
{
	initialize_globals();
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
