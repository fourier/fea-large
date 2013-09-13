#include <float.h>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include <list>
#include <ctype.h>
#include <cassert>
#include <map>

#ifdef WIN32
	#define isnan(arg) _isnan(arg)
#endif


struct Utility
{
	// strains
	static double strain_eng(double k)
	{
		// The engineering strain is defined by: 
		// strain_eng = (L-L0) / L0
		return k - 1; // k = L/L0;
	}
	static double strain_true(double k)
	{
		// The true strain is defined by: 
		// strain_true = ln( L/L0 ) = ln(1 + strain_eng)
		return log(k);
	}
	static double k_from_strain_eng(double strain)
	{
		return strain+1;
	}
	static double k_from_strain_true(double strain)
	{
		return exp(strain);
	}
	
	static
	double stress_true_from_stress_eng( double stress, double k)
	{
		// stress_true = (L/L0) * stress_eng 
		return k*stress;
	}
	static
	double stress_eng_from_stress_true( double stress, double k)
	{
		return stress/k;
	}
private: 
	Utility(){};
};

class Approximation
{
	enum TModelClass
	{
		A,
		B
	};
	// Experimental data
	struct ExperimentalDataT
	{
		double t; // time
		double k; // k = k1
		double s; // s = stress_eng = sigma_1 - Cauchy stress 
	};
	typedef std::list<ExperimentalDataT> ExperimentalDataListT;
public:
	Approximation() : engineering_(false), model_n_(3)
	{
	}
	~Approximation()
	{
	}

	static void usage();
	bool parse_cmdline(int argc, const char* argv[]);
	bool load_file();
	// prepare experimental data for calculations
	// data_ should contain first=k, second=pressure_Tzz
	bool prepare_experimental_data();

	void main_loop(int argc, const char* argv[]);

	double Txx(double k, double param1, double param2);
	double approximation_error(const ExperimentalDataListT& data,
									double param1, double param2);
	static double func(double x1, double x2);
	
	void output_for_gnuplot(double compression_param1, double compression_param2,
									double tension_param1,double tension_param2);
	void calculate_actual_params(const ExperimentalDataListT& data, 
									double& param1,double& param2);
	void enumerative_calculation_params(const ExperimentalDataListT& data,
									double& param1,double& param2);
	void search_calculation_params(const ExperimentalDataListT& data, 
									double& param1,double& param2);

	bool hook_jivs(const ExperimentalDataListT& data,
									double& param1,double& param2);

	bool set_model(const char* model_name);

private:
	std::string data_filename_;
	bool engineering_;
	int model_n_;
	TModelClass model_class_;
	ExperimentalDataListT data_compression_;
	ExperimentalDataListT data_tension_;
};

void Approximation::usage()
{
	using std::cout;
	using std::endl;

	cout << "Program for calculation approxiamtion of uniaxial data in case of nonlinear deformations." << endl;
	cout << "This program would calculate material parameters for any of semilinear models of solids" << endl;
	cout << "(except models A3,B3 and Cn/Dn) using provided data-file with experimental data" << endl;
	cout << "Usage: uniaxial model [-e] filename" << endl;
	cout << "model is the name of model to analyze: A1,A2,A4,A5,B1,B2,B4,B5" << endl;
	cout << "filename is the filename with experimental uniaxial test data." << endl;
	cout << "Such files can be obtained from polymerfem.com site for example" << endl;
	cout << "Inside file where are sequence of lines with space-separated 3 values: " << endl;
	cout << "time strain stress" << endl;
	cout << "Lines started with '#' or '%' will be skipped" << endl;
	cout << "Optional key '-e' mean that data in file is in engineering form: time eng_strain eng_stress" << endl;
	cout << "W/o this key data analyzed as true strain and true stress" << endl;
	cout << "Example of usage: uniaxial A1 -e experimental.txt" << endl;
	cout << "This means calculate parameters for model A1 using experimental data" << endl;
	cout << "from text file 'experimental.txt' where stored engineering strains" << endl;
	cout << "and engineering stresses" << endl;
}

bool Approximation::set_model(const char* model_name)
{
	if ( strlen(model_name) !=  2 )
		return false;
	
	int n = 0;
	char T = model_name[0];
	TModelClass model_class;

	if ( !isdigit(model_name[1]) )
		return false;
	n = model_name[1] - '0';
	if ( T ==  'A' || T ==  'a' )
		model_class = A;
	else if ( T ==  'B' || T ==  'b' )
		model_class = B;
	else 
		return false;
	
	if (!( n == 1 || n ==  2 || n ==  4 || n == 5 ))
		return false;
	
	model_class_ = model_class;
	model_n_ = n;
	return true;
}

bool Approximation::parse_cmdline(int argc, const char* argv[])
{
	if ( argc != 3 && argc != 4 )
		return false;
	
	// check 1st param. it should be 2-characters long , first - letter and second - digit
	if ( !set_model(argv[1]) )
		return false;
	if ( argc == 4 ) // '-e' key specified
	{
		if ( strcmp(argv[2], "-e" ) )
			return false;
		engineering_ = true;
		data_filename_ = argv[3];
	}
	else
		data_filename_ = argv[2];	 

		return true;
}

void Approximation::main_loop(int argc, const char* argv[])
{
	if ( parse_cmdline(argc,argv) )
	{
		if ( load_file() )
		{
			if ( prepare_experimental_data() )
			{
//				double param1 = -2, param2 = 1;
				double compression_param1 = 2.1, compression_param2 = 2.1;
				calculate_actual_params(data_compression_,
												compression_param1,compression_param2);
				double tension_param1 = 2.1, tension_param2 = 2.1;
				calculate_actual_params(data_tension_,
												tension_param1,tension_param2);
				output_for_gnuplot(compression_param1,compression_param2,
												tension_param1,tension_param2);
			}
		}
	}
}

void Approximation::enumerative_calculation_params(
								const ExperimentalDataListT& data,
								double& param1, double& param2)
{
	int steps = 400;
	double param1_start = -20;	
	double param1_end = 20;
	double param2_start = -20;
	double param2_end = 20;

	double step_param1 = (param1_end - param1_start)/steps;
	double step_param2 = (param2_end - param2_start)/steps;

	double best_param1 = param1;
	double best_param2 = param2;
	double min_error = 1e15;
	for ( param1 = param1_start; param1 < param1_end; param1 += step_param1 )
	{
		for ( param2 = param2_start; param2 < param2_end; param2 += step_param2 )
		{
			double err = approximation_error(data,param1,param2);
			if ( err < min_error )
			{
				min_error = err;
				best_param1 = param1;
				best_param2 = param2;
			}
//			std::cout << "param1 = " << param1 << " param2 = " <<param2;
//			std::cout << " error = " << err << std::endl;
		}
	}
	std::cout << "Best params: param1 = " << best_param1;
	std::cout << " param2 = " << best_param2;
	std::cout << std::endl << "Error value = " << min_error;
	std::cout << std::endl;
	param1 = best_param1;
	param2 = best_param2;
}

void Approximation::search_calculation_params(
								const ExperimentalDataListT& data,
								double& param1,double& param2)
{
	hook_jivs(data,param1,param2);	
}


void Approximation::calculate_actual_params(
								const ExperimentalDataListT& data,
								double& param1,double& param2)
{
//	enumerative_calculation_params(data, param1,param2);
	search_calculation_params(data,param1,param2);
}

void Approximation::output_for_gnuplot(
								double compression_param1, double compression_param2,
								double tension_param1, double tension_param2)
{
	FILE* etension_file = fopen("experimental_tension.txt","w");
	FILE* ecompression_file = fopen("experimental_compression.txt","w");
	FILE* tension_file = fopen("calculated_tension.txt","w");
	FILE* compression_file = fopen("calculated_compression.txt","w");
	if ( !etension_file || !ecompression_file || !tension_file || !compression_file)
	{
		fclose(etension_file);
		fclose(ecompression_file);
		fclose(tension_file);
		fclose(compression_file);
		return;
	}

	ExperimentalDataListT::const_iterator it = data_compression_.begin();
	for (; it != data_compression_.end(); ++ it )
	{
		fprintf(ecompression_file,"%f %f\n",it->k,it->s);
		fprintf(compression_file,"%f %f\n",it->k,Txx(it->k,compression_param1,
														compression_param2));
	}
	it = data_tension_.begin();
	for (; it != data_tension_.end(); ++ it )
	{
		fprintf(etension_file,"%f %f\n",it->k,it->s);
		fprintf(tension_file,"%f %f\n",it->k,Txx(it->k,tension_param1,tension_param2));
	}
	fclose(etension_file);
	fclose(ecompression_file);
	fclose(compression_file);
	fclose(tension_file);
}

bool Approximation::load_file()
{
	std::ifstream file(data_filename_.c_str());
	if ( file )
	{
		bool compression = true;
		std::string line;
		ExperimentalDataListT datalist;
		bool extract_next_datalist = false;
		while ( !std::getline(file,line).eof() )
		{
			if ( line.length() )
			{
				// 1st check for comment. Comment symbols: '%' and '#'
				if ( line[0] ==  '%' || line[0] == '#' )
					continue;
				std::stringstream str(line);
				double val1,val2,val3;
				str >> val1 >> val2 >> val3;
				#ifdef LOGGING
				std::cout << "time = " << val1 << " strain = ";
				std::cout << val2 << " stress = " << val3 << std::endl;
				#endif
				// check for compression or tension
				if ( datalist.size() == 1 ) // first value is 0,0,0 almost always
				{
					// check if second value is less than first => compression
					if ( val2 == datalist.begin()->k )
					{
						datalist.begin()->t = val1;
						datalist.begin()->s = val3;
						continue;
					}
					compression = val2 < datalist.begin()->k;
					std::cout << "Experimental data for ";
					if ( compression ) 
						std::cout << "compression found" << std::endl;
					else 
						std::cout << "tension found" << std::endl;
				}
				if ( datalist.size() > 1 ) 
				// try to extract data before hysteresis to avoid loop
				{
					if ( val2 == datalist.rbegin()->k )
					{
						datalist.rbegin()->t = val1;
						datalist.rbegin()->s = val3;
						continue;
					}
					if ( compression ) 
					{
						// for compression hysteresis will occure when strain next value
						// would be less than previous
						if ( val2 > datalist.rbegin()->k )
						{
							data_compression_ = datalist;
							extract_next_datalist = true;
						}
					}
					else
					{ // for tension
						if ( val2 < datalist.rbegin()->k )
						{
							data_tension_ = datalist;
							extract_next_datalist = true;
						}
					}
				}
				ExperimentalDataT data;
				data.t = val1;
				data.k = val2;
				data.s = val3;
				if ( extract_next_datalist )
				{
					if ( compression )
						data_tension_.push_back(data);
					else
						data_compression_.push_back(data);
				}
				else
				{
					datalist.push_back(data);
				}
			}
		}

		return true;
	}
	return false;
}
	
bool Approximation::prepare_experimental_data()
{
	ExperimentalDataListT::iterator it = data_compression_.begin();
	for (; it !=  data_compression_.end(); ++ it )
	{
		double k,s;
		if ( engineering_ ) // engineering strain/stress
		{
			k = Utility::k_from_strain_eng(it->k);
			s = it->s; // from MPa to Pa
		}
		else // true strain/stress
		{
			k = Utility::k_from_strain_true(it->k);
			s = Utility::stress_eng_from_stress_true(it->s,k);
		}
		it->k = k;
		it->s = s;		
	}
	it = data_tension_.begin();
	for (; it !=  data_tension_.end(); ++ it )
	{
		double k,s;
		if ( engineering_ ) // engineering strain/stress
		{
			k = Utility::k_from_strain_eng(it->k);
			s = it->s; // from MPa to Pa
		}
		else // true strain/stress
		{
			k = Utility::k_from_strain_true(it->k);
			s = Utility::stress_eng_from_stress_true(it->s,k);
		}
		it->k = k;
		it->s = s;		
	}
	return true;
}

double Approximation::Txx(double k, 
	double param1, double param2)
{
	double result = 0;
	double p0 = 1.0,C = 0.0;
	if ( model_class_ ==  A )
	{
		double lambda = param1;
		double mu = param2;
		double M = 0,N = 0;
		switch(model_n_)
		{
		case 1:
			M = 0.5*(3.0-1.0/pow(k,2.0)-2*k)*(1.0/pow(k,2.0)-k);
			N = 1./pow(k,2)*(1-1./pow(k,2))-k*(1-k);
			break;
		case 2:
			M = 0.5*(3-1./k-2*sqrt(k))*(1./k-sqrt(k));
			N = 1./k*(1-1./k)-sqrt(k)*(1-sqrt(k));
			break;
		case 4:
			M = 0.5*(k+2./sqrt(k)-3)*(k-1./sqrt(k));
			N = 1./sqrt(k)*(1-1./sqrt(k))-k*(1-k);
			break;
		case 5:
			M = 0.5*(pow(k,2)+2./k-3)*(pow(k,2)-1./k);
			N = 1./k*(1-1./k)-pow(k,2)*(1-pow(k,2));
			break;
		default: 
			return 0;
		};
		result = lambda*M+mu*N;
	}
	else	
	{
		double mu = param1;
		double beta = param2;
		double Z = 0,H = 0;
		switch(model_n_)
		{
		case 1:
			Z = 2*(k-1./pow(k,2));
			H = pow(k,2)-1./k;
			break;
		case 2:
			Z = sqrt(k)-1./k;
			H = k - 1./sqrt(k);
			break;
		case 4:
			Z = k - 1./sqrt(k);
			H = 1./k - sqrt(k);
			break;
		case 5:
			Z = 2*(pow(k,2)-1./k);
			H = 1./pow(k,2)-k;
			break;
		default: 
				return 0;

		};
		result = mu*(1+beta)*Z-mu*(1-beta)*H;			
	}
	return result;
}

double Approximation::approximation_error(
								const ExperimentalDataListT& data,
								double param1, double param2)
{
	double result = 0;
	int N = data.size() -1; // skip first
	ExperimentalDataListT::const_iterator it = data.begin();
	it++;
	for ( ; it != data.end(); ++ it )
	{
		double k = it->k;
		double sigma_calc = Txx(k,param1,param2);
		double sigma_experiment = it->s;
//		std::cout << "experimental = " << pressure_experiment << " calc = ";
//		std::cout << pressure_calc << std::endl;
		result += pow(1.0 - sigma_calc/sigma_experiment,2);
	}
	return sqrt(result/N);
}

double Approximation::func(double x1, double x2)
{
	return 6*pow(x1,2)-4*x1*x2+3*pow(x2,2)+4*sqrt(5)*(x1+2*x2)+22;
}
#define F(data,arg1,arg2) approximation_error(data,arg1,arg2)
//#define F(arg1,arg2) func(arg1,arg2)
bool Approximation::hook_jivs(
								const ExperimentalDataListT& data,
								double& param1,double& param2)
{
	int iter = 0;
	int stagnation_iter = 0;
	static const double gold_ratio=(sqrt(5)-1.0)/2.0;
	const int max_iter = 20000;
	const int max_stagnation_iter = 1000;
	const double eps = 0.00001;
	const double gamma = 2.0;
	double a = 1.5;
	double p = eps;
	double q = 1;
	double alpha,beta;
	double b1 = 1.,b2 = 1.;
	double _x1 = param1, _x2 = param2;
	std::cout << std::endl;
	do
	{
		iter ++;
		const double x1 = _x1, x2 = _x2;
		double f = F(data,_x1,_x2);
		double f1_plus = F(data,_x1+b1,_x2);
		double f1_minus = F(data,_x1-b1,_x2);
#ifdef LOGGING
		std::cout << "x1_="<<_x1<<" x2_="<<_x2<<std::endl;
		std::cout << "f = " << f<< " f1_plus = " << f1_plus << " f1_minus = " << f1_minus<< std::endl;
#endif
		if ( f1_plus < f && f1_plus <= f1_minus )
		{
#ifdef LOGGING
			std::cout << "f1_plus < f && f1_plus < f_minus" << std::endl;
#endif
			_x1 = _x1+b1;
#ifdef LOGGING
			std::cout << "_x1="<<_x1<<std::endl;
#endif
		}
		else if ( f1_minus < f && f1_minus < f1_plus )
		{
#ifdef LOGGING
			std::cout <<"f1_minus < f && f1_minus < f1_plus" << std::endl;	
#endif
			_x1 = _x1-b1;
#ifdef LOGGING
			std::cout << "_x1="<<_x1<<std::endl;
#endif
		}
		double f2_plus = F(data,_x1,_x2+b2);
		double f2_minus = F(data,_x1,_x2-b2);
		f = F(data,_x1,_x2);
#ifdef LOGGING
		std::cout << "f = " << f<< " f2_plus = " << f2_plus << " f2_minus = " << f2_minus<< std::endl;
#endif
		if ( f2_plus < f && f2_plus <= f2_minus ) 
		{
#ifdef LOGGING
				std::cout<<"f2_plus < f && f2_plus < f2_minus" << std::endl;
#endif
				_x2 = _x2+b2;
#ifdef LOGGING
				std::cout <<"_x2="<<_x2<<std::endl;
#endif
		}
		else if ( f2_minus < f && f2_minus < f2_plus )
		{
#ifdef LOGGING
			std::cout<<"f2_minus < f && f2_minus < f2_plus" << std::endl;
#endif
			_x2 = _x2-b2;
#ifdef LOGGING
			std::cout <<"_x2="<<_x2<<std::endl;
#endif
		}
		if ( _x1 == x1 && _x2 == x2 )
		{
			b1 = b1/gamma;
			b2 = b2/gamma;
			_x1 = x1;
			_x2 = x2;
#ifdef LOGGING
			std::cout << "stagnation" << std::endl;
#endif
			stagnation_iter ++;
			if ( stagnation_iter > max_stagnation_iter )
			{
				double err = sqrt(pow(x1-_x1,2)+pow(x2-_x2,2));
				std::cout << "err = " << err << std::endl;
				param1 = _x1;
				param2 = _x2;
				break;
			}
			continue;
		}
		double err = sqrt(pow(x1-_x1,2)+pow(x2-_x2,2));
		if ( err < eps || iter > max_iter )
		{
			std::cout << "err = " << err << std::endl;
			param1 = _x1;
			param2 = _x2;
			break;
		}
		while (fabs(p-q)>eps)
		{
			alpha=gold_ratio*p+(1.0-gold_ratio)*q;
			beta=(1.0-gold_ratio)*p+gold_ratio*q;
			if (F(data,x1+(_x1-x1)*alpha,x2+(_x2-x2)*alpha) > 
					F(data,x1+(_x1-x1)*beta,x2+(_x2-x2)*beta))
			{
				p=alpha;
				alpha=beta;
				beta=(1-gold_ratio)*p+gold_ratio*q;
			}
			else
			{
				q=beta;
				beta=alpha;
				alpha=gold_ratio*p+(1-gold_ratio)*q;
			}
		}
		a=(p+q)/2.0;
		_x1 = x1 + a*(_x1 - x1);
		_x2 = x2 + a*(_x2 - x2);
#ifdef LOGGING		
		std::cout << "iteration " << iter << " finished with values: " << std::endl;
		std::cout << "x1 = " << _x1 << " x2 = " << _x2 << std::endl;
#endif	
	} while(true);
	std::cout << std::endl;
	std::cout << "Number of iterations = " << iter << std::endl;
	std::cout << "x1 = " << param1 << " x2 = " << param2;
	std::cout << " f = " << F(data,param1,param2) << std::endl;
	return true;
}

int main(int argc, const char* argv[])
{
	Approximation appr;
	Approximation::usage();
	appr.main_loop(argc, argv);
	return 0;
}
