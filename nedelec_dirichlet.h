/**
 * Dirichlet boundary conditions for the rot-rot-discretizations based on the
 * Nedelec-type-1 elements (the Whitney-1 forms).
 *
 * Created on: 01.02.2013
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_DIRICHLET__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_DIRICHLET__

#include <map>
#include <vector>

// ug4 headers
#include "common/common.h"

// library-specific headers
#include "lib_grid/lg_base.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_grid/tools/subset_handler_interface.h"

#include "em_material.h"

namespace ug{
namespace Electromagnetism{

// Dirichlet BC for a rot-rot operator
/**
 * This class implements a Dirichlet boundary condition for the rot-rot operator.
 * The boundary condition has the form
 * \f{eqnarray*}{
 *  \mathbf{n} \times \mathbf{E} = \mathbf{n} \times \mathbf{E}_\mathrm{bc},
 * \f}
 * where \f$ \mathbf{E} \f$ is the unknown function, \f$ \mathbf{n} \f$ the
 * unit outer normal to the boundary segment (or patch), \f$ \mathbf{E}_\mathrm{bc} \f$
 * the given boundary condition. Note that \f$ \mathbf{E}_\mathrm{bc} \f$ specified
 * by the user is a vector (although the Nedelec DoFs of \f$ \mathbf{E} \f$ are
 * scalar).
 *
 * \remark For the time-harmonic problems (or other problems with several
 * components in every dof), the whole boundary segments are considered as
 * Dirichlet boundary for all the specified components if at least one of
 * the component is mentioned in the object of this class. For the components
 * that are not explicitly mentioned in the "add" functions, zero
 * Dirichlet values are imposed (so called "implicit specification").
 */
template <typename TDomain, typename TAlgebra>
class NedelecDirichletBC
	: public EMDirichlet<TDomain, TAlgebra>
{
public:
	///	Base type
		typedef IDomainConstraint<TDomain, TAlgebra> base_type;
		
	/// This type
		typedef NedelecDirichletBC<TDomain, TAlgebra> this_type;

	///	Type of domain
		typedef TDomain domain_type;

	///	World dimension
		static const int dim = domain_type::dim;

	///	Type of position coordinates (e.g. position_type)
		typedef typename domain_type::position_type position_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;
	
private:
	/// Iterator over edges
		typedef DoFDistribution::traits<Edge>::const_iterator t_edge_iterator;
	
public:
	///	Constructor
		NedelecDirichletBC
		(
			const char * funcNames ///< [in] names of the functions for the BC
		)
		:	m_vDirichletFunc (TokenizeString (funcNames))
		{};
		
	/// Adds a zero Dirichlet value on a subset
		void add_0 (const char * subsets);
	
	/// Adds a constant Dirichlet value on a subset
	///\{
		void add (MathVector<dim> & value, const char * function, const char * subsets);
		void add (std::vector<number> & vValue, const char * function, const char * subsets);
	///\}
	
	/// Adds a position-dependent Dirichlet BC based on UserData
	///\{
		void add (SmartPtr<UserData<MathVector<dim>, dim> > & func, const char * function, const char * subsets);
#ifdef UG_FOR_LUA
		void add (const char * name, const char * function, const char * subsets);
#endif
	///\}

	/// composes an array of subset indices of the Dirichlet boundary
		virtual void get_dirichlet_subsets
		(
			SubsetGroup & dirichlet_ssgrp ///< the group to update
		) const;
	
private:
// Data, auxiliary types and auxiliary functions:
	
	/// Coordinates of the grid vertices
	typename domain_type::position_accessor_type m_aaPos;
	
	/// Functions for which the Dirichlet conditions are imposed
	std::vector<std::string> m_vDirichletFunc;
	
	/** All the subsets with the Dirichlet BC.
	 * The map assignes the functions that have not been mentioned
	 * in the explicit specifications to the subsets with the
	 * Dirichlet BC.
	 */
	std::map<std::string, std::vector<std::string> > m_mDirichletSS;
	
	/// Structure for the Dirichlet BC that is constant over a patch
	struct TConstBC
	{
		MathVector<dim> bcValue; // the boundary condition
		std::string fctName; // grid function component name
		size_t fct; // grid function component
		std::string ssName; // subset name
		SubsetGroup ssGrp; // subset group
		
		TConstBC
		(
			MathVector<dim> & the_bcValue,
			std::string the_fctName,
			std::string the_ssName
		)
		: bcValue (the_bcValue), fctName (the_fctName), ssName (the_ssName)
		{}
		
		inline number operator ()
		(
			EdgeBase * edge, // the edge
			int si, // subset index
			typename domain_type::position_accessor_type & aaPos, // vertex coordinates
			number time // the time argument
		)
		{
			// Get the vector pointing along the edge from its beginning to its end:
			position_type edgeVector = aaPos [edge->vertex(1)];
			edgeVector -= aaPos [edge->vertex(0)];
			
			// Return the scalar product:
			return bcValue * edgeVector;
		}
	};
	
	/// Constant values as Dirichlet boundary conditions
	///\{
	std::vector<TConstBC> m_vConstBCData;
	std::map<int, std::vector<TConstBC *> > m_mConstBC;
	///\}
	
	/// Structure for the Dirichlet BC that are given by a function
	struct TUserDataBC
	{
		SmartPtr<UserData<MathVector<dim>, dim> > spBCFunc;
		std::string fctName; // grid function component name
		size_t fct; // grid function component
		std::string ssName; // subset name
		SubsetGroup ssGrp; // subset group
		
		TUserDataBC
		(
			SmartPtr<UserData<MathVector<dim>, dim> > the_spFunc,
			std::string the_fctName,
			std::string the_ssName
		)
		: spBCFunc (the_spFunc), fctName (the_fctName), ssName (the_ssName)
		{}
		
		inline number operator ()
		(
			EdgeBase * edge, // the edge
			int si, // subset index
			typename domain_type::position_accessor_type & aaPos, // vertex coordinates
			number time // the time argument
		)
		{
			number dofValue;
			MathVector<dim> func_value;
			
			// Get the vector pointing along the edge from its beginning to its end:
			position_type edgeVector = aaPos [edge->vertex(1)];
			edgeVector -= aaPos [edge->vertex(0)];
			
			// Evaluate the function
			(*spBCFunc) (func_value, aaPos [edge->vertex(0)], time, si);
			dofValue = func_value * edgeVector;
			
			(*spBCFunc) (func_value, aaPos [edge->vertex(1)], time, si);
			dofValue += func_value * edgeVector;
			
			return dofValue / 2;
		}
	};

	/// Dirichlet boundary conditions given as functions
	///\{
	std::vector<TUserDataBC> m_vUserDataBCData;
	std::map<int, std::vector<TUserDataBC *> > m_mUserDataBC;
	///\}
	
	/// List of all the Dirichlet subsets with implicitly specified Dirichlet BC: zero values
	std::map<int, FunctionGroup> m_mZeroBC;
	
	/// inserts a subset to the list of the Dirichlet subsets
	void add_subset
	(
		const char* f_name, ///< name of the function
		const char* ss_names ///< names of the subsets
	);
	
	///	sets the approximation space to work on
	void set_approximation_space
	(
		SmartPtr<ApproximationSpace<TDomain> > approxSpace
	)
	{
		base_type::set_approximation_space(approxSpace);
		m_aaPos = base_type::m_spApproxSpace->domain()->position_accessor();
	}

	/// Verifies whether the string input represents legal data
	void check_functions_and_subsets
	(
		FunctionGroup & functionGroup,
		SubsetGroup & subsetGroup
	);
	
	/// Composes the list of the implicit BC
	void extract_implicit ();
	
	/// Compiles the BC data into the map assigning the data to the subset id's
	///\{
	template <typename TUserData>
	void extract_data
	(
		std::map<int, std::vector<TUserData*> >& mvUserDataBndSegment,
		std::vector<TUserData>& vUserData
	);
	void extract_data ()
	{
		if (! base_type::m_spApproxSpace.valid ())
			UG_THROW ("NedelecDirichletBC: Approximation space not set.");
		extract_data (m_mConstBC, m_vConstBCData);
		extract_data (m_mUserDataBC, m_vUserDataBCData);
		extract_implicit ();
	}
	///\}
	
	/// Sets the specified Dirichlet values in the solution
	template <typename TUserData>
	void adjust_solution
	(
		const std::vector<TUserData *> & vUserData,
		int si,
		vector_type& u,
		ConstSmartPtr<DoFDistribution> dd,
		number time
	);
	
	/// Sets the zero Dirichlet values for the implicitly specified BC
	void adjust_solution_implicit
	(
		vector_type& u,
		ConstSmartPtr<DoFDistribution> dd,
		number time
	);

public:
// The interface:
	
	/// sets a unity row for all dirichlet indices
	void adjust_jacobian
	(
		matrix_type & J,
		const vector_type & u,
		ConstSmartPtr<DoFDistribution> dd,
		number time = 0.0,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
		const number s_a0 = 1.0
	);

	/// sets a zero value in the defect for all dirichlet indices
	void adjust_defect
	(
		vector_type & d,
		const vector_type & u,
		ConstSmartPtr<DoFDistribution> dd,
		number time = 0.0,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = NULL,
		const std::vector<number> * vScaleMass = NULL,
		const std::vector<number> * vScaleStiff = NULL
	);

	/// sets the dirichlet value in the solution for all dirichlet indices
	void adjust_solution
	(
		vector_type & u,
		ConstSmartPtr<DoFDistribution> dd,
		number time = 0.0
	);

	///	sets unity rows in A and dirichlet values in right-hand side b
	void adjust_linear
	(
		matrix_type & A,
		vector_type & b,
		ConstSmartPtr<DoFDistribution> dd,
		number time = 0.0
	)
	{
		this_type::adjust_jacobian (A, b, dd, time); // the 2nd arg. is dummy: A does not depend on u
		this_type::adjust_solution (b, dd, time);
	};

	///	sets the dirichlet value in the right-hand side
	void adjust_rhs
	(
		vector_type & b,
		const vector_type & u,
		ConstSmartPtr<DoFDistribution> dd,
		number time = 0.0
	)
	{
		this_type::adjust_solution (b, dd, time);
	};

	///	returns the type of the constraints
	int type () const {return CT_DIRICHLET;}
	
};

} // end namespace Electromagnetism
} // end namespace ug

// Implementation of the functions:
#include "nedelec_dirichlet_impl.h"

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_DIRICHLET__

/* End of File */
