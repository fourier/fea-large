//
// Main project file. Defines an entry point to the application
//


// pragmas
#pragma warning( disable : 4267 )

// precompiled headers
#include "std.h"

// Local includes
#include "GeometryData.h"
#include "gauss_coeffs.h"
#include "elements.h"
#include "dataloaderwrapper.h"
#include "geometrydatawrapper.h"
#include "application.h"


application app;

int main(int argc,const char* argv[])
{
//	test_element_axisymmetric_triangle_p6();
	app.run(argc,argv);
	return 0;
}
