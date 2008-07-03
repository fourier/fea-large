// Precomiled headers
#include "std.h"

// local includes
#include "axisymmetric6p3.h"
#include "functions.h"


void axisymmetric6p3::linear_construct_local_matrix(size_type el, Matrix& m) const
{
	assert(m.size1() == m.size2());
	const size_type size = m.size1();

	boost::array<Tensor4Rank::index, 4> shape = {{ MAX_DOF, MAX_DOF, MAX_DOF, MAX_DOF }};
	Tensor4Rank elasticity_tensor(shape);
	// iterate over gauss nodes
	for ( size_type gauss = 0; gauss < GaussNodes::GaussNumber; ++ gauss )
	{
		// integration: Kil=2*PI*Kil*r

		// First of all constuct [B] matrix
		MATRIX(B,Element::VoigtNumber,size);
		for ( size_type i = 0; i < Element::VoigtNumber; ++ i )
			for ( size_type j = 0; j < size; ++ j )
				B(i,j) = b_matrix_proxy(el,gauss,i,j);
		// next construct C tensor in gauss node
		model_->construct_ctensor(elasticity_tensor,graddefs_gauss_nodes_[el][gauss]);
		MATRIX(C,Element::VoigtNumber,Element::VoigtNumber);
		construct_elasticity_matrix(C,elasticity_tensor);

		// next calculate inversion of det(F)
#ifdef DETF
		const value_type inv_detF = 1.0/det3x3(graddefs_gauss_nodes_[el][gauss]);
		// next assemble B^T * C * B
		// optimization: move out multiplier
		const value_type multiplier = 2*PI*GaussNodes::weights[gauss]*gauss_nodes_[el][gauss].dof[0]*inv_detF;
#else
		const value_type multiplier = 2*PI*GaussNodes::weights[gauss]*gauss_nodes_[el][gauss].dof[0];	
#endif
		for ( size_type i = 0; i < size; ++ i )
		{
			for ( size_type j = 0; j < Element::VoigtNumber; ++ j )
			{
				for ( size_type k = 0; k < Element::VoigtNumber; ++ k )
				{
					for ( size_type l = 0; l < size; ++ l )
					{
//						m(i,l) += B(j,i)*symm_tensor4rank_matrix_proxy(elasticity_tensor,j,k)*B(k,l)*multiplier;
						m(i,l) += B(j,i)*C(j,k)*B(k,l)*multiplier;
					}
				}
			}
		}
	}
}

void axisymmetric6p3::deformation_gradient(size_type el,const Node& node, Matrix& F) const
{
	assert( F.size1() == F.size2() && F.size1() == MAX_DOF );
	F.clear();	
	eye(F);
	const Element& elem = initial_elements_[el];
	for ( size_type k = 0; k < Element::NodesNumber; ++ k ) 
	{
		const size_type dof_pos = k*Element::DofNumber;
		const size_type index_r = index_(el,dof_pos);
		const size_type index_z = index_(el,dof_pos+1);
		const value_type r = global_solution_[index_r];
		const value_type z = global_solution_[index_z];
		const value_type dr = elem.dform(k,0,node);
		const value_type dz = elem.dform(k,1,node);
		const value_type f = elem.form(k,node)/node.dof[0];
		F(0,0) += dr*r;
		F(1,1) += dz*z;
		F(2,2) += f*r;
		F(0,1) += dz*r;
		F(1,0) += dr*z;
	}
}

void axisymmetric6p3::deformation_gradient(size_type el,const ElementsArray& elements,const Node& node, Matrix& F) const
{
	assert( F.size1() == F.size2() && F.size1() == MAX_DOF );
	F.clear();	
	eye(F);
	const Element& elem = initial_elements_[el];
	const Element& deformed = elements[el];
	for ( size_type k = 0; k < Element::NodesNumber; ++ k ) 
	{
		const value_type r = deformed.node(k).dof[0] - elem.node(k).dof[0];
		const value_type z = deformed.node(k).dof[1] - elem.node(k).dof[1];
		const value_type dr = elem.dform(k,0,node);
		const value_type dz = elem.dform(k,1,node);
		const value_type f = elem.form(k,node)/node.dof[0];
		F(0,0) += dr*r;
		F(1,1) += dz*z;
		F(2,2) += f*r;
		F(0,1) += dz*r;
		F(1,0) += dr*z;
	}
}

void axisymmetric6p3::construct_ksigma(size_type el, Matrix& m) const
{
	using boost::numeric::ublas::prod;
	using boost::numeric::ublas::inner_prod;
	using boost::numeric::ublas::trans;

	assert(m.size1() == m.size2());
	const size_type size = m.size1();

	// first, iterate over gauss nodes.
	for ( size_type gauss = 0; gauss < GaussNodes::GaussNumber; ++ gauss )
	{
		MATRIX(Sigma,3,3);
		model_->stress_cauchy(graddefs_gauss_nodes_[el][gauss],Sigma);
		std::vector<Vector> B_vectors;
		B_vectors.resize(Element::NodesNumber);
		for ( size_type i = 0; i < Element::NodesNumber; ++ i )
		{
			B_vectors[i].resize(3);
			for ( size_type j = 0; j < 3; ++ j )
				B_vectors[i][j] = b_vector_proxy(el,gauss,j,i);
		}
		// optimization: move out multiplier
		// for axisymmetric task 2*PI*r
		const value_type multiplier = 2*PI*GaussNodes::weights[gauss]*gauss_nodes_[el][gauss].dof[0];
		for ( size_type k = 0; k < Element::NodesNumber; ++ k )
		{
			for ( size_type l = 0; l < Element::NodesNumber; ++ l )
			{
				Vector tmp( prod(trans(B_vectors[k]),Sigma) );
				value_type value = inner_prod(tmp,B_vectors[l]);		
				const size_type offset_row = k*Element::DofNumber;
				const size_type offset_column = l*Element::DofNumber;
				value = value*multiplier;
				m(offset_row,offset_column) += value;	
				m(offset_row+1,offset_column+1) += value;
			}
		}
	}
}

void axisymmetric6p3::construct_residual_force(size_type el,Vector& v) const
{
	using boost::numeric::ublas::prod;
	using boost::numeric::ublas::inner_prod;
	using boost::numeric::ublas::trans;

	size_type size = v.size();

	// first, iterate over gauss nodes
	for ( size_type gauss = 0; gauss < GaussNodes::GaussNumber; ++ gauss )
	{
		// First of all constuct [B] matrix
		MATRIX(B,Element::VoigtNumber,size);
			for ( size_type i = 0; i < Element::VoigtNumber; ++ i )
				for ( size_type j = 0; j < size; ++ j )
					B(i,j) = b_matrix_proxy(el,gauss,i,j);
		MATRIX(Sigma,Element::VoigtNumber,Element::VoigtNumber);
		model_->stress_cauchy(graddefs_gauss_nodes_[el][gauss],Sigma);
		VECTOR(S,Element::VoigtNumber);
		for ( size_type l = 0; l < Element::VoigtNumber; ++ l )
		{
			S[l] = Sigma(g_voigt_mapping.find(l)->second.first,g_voigt_mapping.find(l)->second.second);
		}
		v += GaussNodes::weights[gauss]*prod(trans(B),S)
			* 2*PI*gauss_nodes_[el][gauss].dof[0];
	}
}


void axisymmetric6p3::export_msh(const char* filename,const ElementsArray& elements) const
{
	assert(initial_elements_.size() == elements.size());
	std::ofstream f(filename, std::ios_base::out | std::ios_base::trunc );
	if ( f )
	{
		for ( size_type i = 0; i < initial_elements_.size(); ++ i )
		{
			size_type j = 0;
			size_type k = 0;
			for ( ; j < 3; ++ j )
			{
				for ( k = 0; k < 2; ++ k )
				{
					f << initial_elements_[i].node(j).dof[k] << ",";
				}
			}
			for ( j = 0; j < 3; ++ j )
			{
				for ( k = 0; k < 2; ++ k )
				{
					f << elements[i].node(j).dof[k] << ",";
				}
			}
			f << "0,";
			Node nod = initial_elements_[i].center();
			MATRIX(F,3,3);
			deformation_gradient(i,elements,nod,F);
			MATRIX(S,3,3);
			model_->stress_cauchy(F,S);
			f << S(0,0) << "," << S(1,1) << "," << S(2,2) << "," << S(0,1) << std::endl;		
		}
		f.flush();
		f.close();
	}
}

void axisymmetric6p3::postprocessing(const char* filename,const ElementsArray& elements) const
{
//	char fname[50];
//	sprintf(fname,filename,parameter);
	export_msh(filename,elements);
}


