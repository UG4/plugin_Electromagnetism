/*
 * Implementation of the functions for the Dirichlet boundary conditions
 * for the Nedelec-based discretizations of the rot-rot operator.
 *
 * Created on: 01.02.2013
 * Author: D. Logashenko
 */

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace Electromagnetism{

/**
 * Adds a subset to the list of the Dirichlet boundaries.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::add_subset
(
	const char* f_name, ///< name of the function
	const char* ss_names ///< names of the subsets
)
{
	typedef std::vector<std::string> t_f_vec;
	
	std::vector<std::string> subsets;
	TokenizeString (std::string (ss_names), subsets);
	
//	Check the data
	std::string s_f_name (f_name);
	size_t n_func;
	if ((n_func = TokenizeString(s_f_name).size ()) > 1)
		UG_THROW("NedelecDirichletBC:"
					" Only single functions function allowed in the specification of"
					" Dirichlet values, but the following functions given:"
					<< f_name);
	
//	Loop the subsets
	for (size_t i = 0; i < subsets.size (); i++)
	{
	//	If the map entry does not exist, create and initialize with all the functions
		if (m_mDirichletSS.find (subsets [i]) == m_mDirichletSS.end ())
			m_mDirichletSS [subsets [i]] = m_vDirichletFunc;
	//	Exclude the specified function from the list of the implicit specifications
		std::vector<std::string> & imp_funcs = m_mDirichletSS [subsets [i]];
		if (n_func == 0) continue;
		t_f_vec::iterator f_entry = std::find (imp_funcs.begin (), imp_funcs.end (), s_f_name);
		if (f_entry == imp_funcs.end ())
			UG_THROW ("NedelecDirichletBC: Wrong speficifation of the Dirichlet BC for subset "
				<< subsets [i] << ": func " << f_name <<
				" not present in the object or is specified twice for the subsets.");
		imp_funcs.erase (f_entry);
	}
}

/**
 * Adds a zero Dirichlet BC
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::add_0
(
	const char* subsets ///< the subset of the boundary patch
)
{
	add_subset ("", subsets);
}

/**
 * Adds a constant Dirichlet value on a given subset.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::add
(
	MathVector<dim> & value, ///< values of the BC (a dim-vector)
	const char* function, ///< the name of the function to impose the value at
	const char* subsets ///< the subset of the boundary patch
)
{
	add_subset (function, subsets);
	m_vConstBCData.push_back (TConstBC (value, function, subsets));
}

/**
 * Adds a constant Dirichlet value on a given subset.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::add
(
	std::vector<number> & vValue, ///< values of the BC (a dim-vector)
	const char* function, ///< the name of the function to impose the value at
	const char* subsets ///< the subset of the boundary patch
)
{
	if (vValue.size () != (size_t) dim)
		UG_THROW ("NedelecDirichletBC:"
					"Wrong dimensionality of a Dirichlet BC. Specify " << dim << " values.");
	MathVector<dim> value;
	for (size_t i = 0; i < (size_t) dim; i++) value[i] = vValue[i];
	add (value, function, subsets);
}

/**
 * Adds a constant Dirichlet value on a given subset.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::add
(
	SmartPtr<UserData<MathVector<dim>, dim> > & func, ///< values of the BC (as a function)
	const char* function, ///< the name of the function to impose the value at
	const char* subsets ///< the subset of the boundary patch
)
{
	add_subset (function, subsets);
	m_vUserDataBCData.push_back (TUserDataBC (func, function, subsets));
}

#ifdef UG_FOR_LUA
/**
 * Adds a constant Dirichlet value on a given subset (in a lua script).
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::add
(
	const char* name, ///< values of the BC (as a function specified by its name)
	const char* function, ///< the name of the function to impose the value at
	const char* subsets ///< the subset of the boundary patch
)
{
	if (! LuaUserData<MathVector<dim>, dim>::check_callback_returns (name))
		UG_THROW ("NedelecDirichletBC: Illegal BC specification."
					" A " << dim << "d vector-valued function should be specified");
	SmartPtr<UserData<MathVector<dim>, dim> > sp_userData =
			LuaUserDataFactory<MathVector<dim>, dim>::create (name);
	add (sp_userData, function, subsets);
}
#endif

/**
 * Composes an array of subset indices of the Dirichlet boundary
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::get_dirichlet_subsets
(
	SubsetGroup & dirichlet_ssgrp ///< the group to update
) const
{
	typedef std::map<std::string, std::vector<std::string> > t_ss_map;
	
// Loop the subset names
	t_ss_map::const_iterator iterEnd = m_mDirichletSS.end ();
	for (t_ss_map::const_iterator iter = m_mDirichletSS.begin (); iter != iterEnd; ++iter)
		dirichlet_ssgrp.add (iter->first);
}

/**
 * Verifies whether the string input represents proper data in the sence
 * of the approximation space.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::check_functions_and_subsets
(
	FunctionGroup& functionGroup, ///< the function group (should represent one scalar function)
	SubsetGroup& subsetGroup ///< subset group where the boundary condition is imposed
)
{
//	only number of functions allowed
	if (functionGroup.size () != 1)
		UG_THROW("NedelecDirichletBC:"
					" Only single functions function allowed in the specification of"
					" Dirichlet values, but the following functions given:"
					<< functionGroup);

//	get subsethandler
	ConstSmartPtr<ISubsetHandler> pSH = base_type::m_spApproxSpace->subset_handler ();

// 	loop subsets
	for(size_t si = 0; si < subsetGroup.size (); ++si)
	{
	//	get the subset and function indices
		const int subsetIndex = subsetGroup[si];
		const size_t fct = functionGroup[0];

	//	check that subsetIndex is valid
		if (subsetIndex < 0 || subsetIndex >= pSH->num_subsets ())
			UG_THROW ("NedelecDirichletBC:"
							" Invalid Subset Index " << subsetIndex << ". (Valid is"
							" 0, .. , " << pSH->num_subsets() <<").");

	// 	check if function exist
		if (fct >= base_type::m_spApproxSpace->function_pattern()->num_fct ())
			UG_THROW ("NedelecDirichletBC:"
						" Function " << fct << " does not exist in pattern.");

	// 	check that function is defined for segment
		if (!base_type::m_spApproxSpace->function_pattern()->is_def_in_subset (fct, subsetIndex))
			UG_THROW ("NedelecDirichletBC:"
							" Function " << fct << " not defined on subset " << subsetIndex);
	}
}

/**
 * Compiles the BC data into the map assigning the data to the subset id's
 * (i.e. fills the map) for a given BC data set.
 */
template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void NedelecDirichletBC<TDomain, TAlgebra>::extract_data
(
	std::map<int, std::vector<TUserData *> > & mvUserDataBndSegment, ///< the map to fill
	std::vector<TUserData> & vUserData /// the bc data set
)
{
//	Clear the extracted data
	mvUserDataBndSegment.clear ();

//	Fill the map
	for (size_t i = 0; i < vUserData.size (); ++i)
	{
		FunctionGroup fctGrp;
		
	//	create Function Group and Subset Group
		try
		{
			vUserData[i].ssGrp = base_type::m_spApproxSpace->subset_grp_by_name (vUserData[i].ssName.c_str());
		}
		UG_CATCH_THROW(" Subsets '" << vUserData[i].ssName << "' not all contained in ApproximationSpace.");

		try
		{
			fctGrp = base_type::m_spApproxSpace->fct_grp_by_name (vUserData[i].fctName.c_str());
		}
		UG_CATCH_THROW (" Functions '" << vUserData[i].fctName << "' not all contained in ApproximationSpace.");

	//	check functions and subsets
		check_functions_and_subsets (fctGrp, vUserData[i].ssGrp);

	//	set the function
		vUserData[i].fct = fctGrp[0];

	// 	loop subsets
		for (size_t si = 0; si < vUserData[i].ssGrp.size (); ++si)
			mvUserDataBndSegment[vUserData[i].ssGrp[si]].push_back (&vUserData[i]);
	}
}

/**
 * Composes the list of the implicit BC
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::extract_implicit ()
{
	typedef std::map<std::string, std::vector<std::string> > t_ss_map;
	
// Loop the subset names
	t_ss_map::const_iterator iterEnd = m_mDirichletSS.end ();
	for (t_ss_map::const_iterator iter = m_mDirichletSS.begin (); iter != iterEnd; ++iter)
	if (iter->second.size () != 0) // only if there are implicitly defined BCs
	{
	//	Get the subset group
		SubsetGroup ssGrp;
		try
		{
			ssGrp = base_type::m_spApproxSpace->subset_grp_by_name (iter->first.c_str ());
		}
		UG_CATCH_THROW(" Subsets '" << iter->first << "' not all contained in ApproximationSpace.");
		
	// 	Loop subsets
		for (size_t si = 0; si < ssGrp.size (); si++)
		{
			FunctionGroup & f_grp = m_mZeroBC [si];
			f_grp.set_function_pattern (base_type::m_spApproxSpace->function_pattern ());
			f_grp.add (iter->second);
		}
	}
}

/**
 * Sets the Dirichlet-rows of a given Jacobian to the identity.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::adjust_jacobian
(
	matrix_type & J,
	const vector_type & u,
	ConstSmartPtr<DoFDistribution> dd,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const number s_a0
)
{
	std::vector<DoFIndex>  multInd;
	SubsetGroup ssGrp (base_type::m_spApproxSpace->subset_handler ());
	FunctionGroup fctGrp (base_type::m_spApproxSpace->function_pattern (), m_vDirichletFunc);
	
//	Get the Dirichlet subsets
	NedelecDirichletBC::get_dirichlet_subsets (ssGrp);
	
//	Loop over the subsets
	for (size_t i = 0; i < ssGrp.size (); i++)
	{
		int si = ssGrp [i];
		
	//	Loop the edges
		t_edge_iterator iterEnd = dd->end<Edge> (si);
		for (t_edge_iterator iter = dd->begin<Edge> (si); iter != iterEnd; iter++)
		//	Loop functions
			for(size_t fct_i = 0; fct_i < fctGrp.size (); fct_i++)
			{
			//	Get the multiindices
				if (dd->inner_dof_indices (*iter, fctGrp [fct_i], multInd) != 1)
					UG_THROW ("NedelecDirichletBC: "
						"More than one DoF per edge. Not the Nedelec-type-1 element?");
	
			//	Set the Dirichlet row
				SetDirichletRow (J, multInd[0]);
			}
	}
}

/**
 * Sets the Dirichlet-entries of a given defect vector to 0.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::adjust_defect
(
	vector_type& d,
	const vector_type& u,
	ConstSmartPtr<DoFDistribution> dd,
	number time,
	ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
	const std::vector<number> * vScaleMass,
	const std::vector<number> * vScaleStiff
)
{
	std::vector<DoFIndex>  multInd;
	SubsetGroup ssGrp (base_type::m_spApproxSpace->subset_handler ());
	FunctionGroup fctGrp (base_type::m_spApproxSpace->function_pattern (), m_vDirichletFunc);
	
//	Get the Dirichlet subsets
	NedelecDirichletBC::get_dirichlet_subsets (ssGrp);
	
//	Loop over the subsets
	for (size_t i = 0; i < ssGrp.size (); i++)
	{
		int si = ssGrp [i];
		
	//	Loop the edges
		t_edge_iterator iterEnd = dd->end<Edge> (si);
		for (t_edge_iterator iter = dd->begin<Edge> (si); iter != iterEnd; iter++)
		//	Loop functions
			for(size_t fct_i = 0; fct_i < fctGrp.size (); fct_i++)
			{
			//	Get the multiindices
				if (dd->inner_dof_indices (*iter, fctGrp [fct_i], multInd) != 1)
					UG_THROW ("NedelecDirichletBC: "
						"More than one DoF per edge. Not the Nedelec-type-1 element?");
	
			//	Set the Dirichlet entry
				DoFRef (d, multInd[0]) = 0.0;
			}
	}
}

/**
 * Sets the Dirichlet-entries of a given solution vector to the boundary values.
 */
template <typename TDomain, typename TAlgebra>
template <typename TUserData>
void NedelecDirichletBC<TDomain, TAlgebra>::adjust_solution
(
	const std::vector<TUserData *> & vUserData, ///< the bc data set
	int si, ///< subset index
    vector_type& u, ///< solution vector to correct
    ConstSmartPtr<DoFDistribution> dd, ///< the DoF distribution
    number time ///< time argument
)
{
	std::vector<DoFIndex>  multInd;

//	Loop the edges
	t_edge_iterator iterEnd = dd->end<Edge> (si);
	for (t_edge_iterator iter = dd->begin<Edge> (si); iter != iterEnd; iter++)
	//	Loop dirichlet functions on this segment
		for(size_t i = 0; i < vUserData.size (); ++i)
		{
			EdgeBase * pEdge = *iter;
		//	Get the multiindices
			if (dd->inner_dof_indices (pEdge, vUserData[i]->fct, multInd) != 1)
				UG_THROW ("NedelecDirichletBC:"
					"More than one DoF per edge. Not the Nedelec-type-1 element?");

		//	Set the Dirichlet entry
			DoFRef (u, multInd[0]) = (* vUserData[i]) (pEdge, si, m_aaPos, time);
		}
}

/**
 * Sets the zero Dirichlet values for the DoFs that have not been mentioned
 * in the explicit specification.
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::adjust_solution_implicit
(
    vector_type & u, ///< solution vector to correct
    ConstSmartPtr<DoFDistribution> dd, ///< the DoF distribution
    number time ///< time argument
)
{
	typedef std::map<int, FunctionGroup> t_ss_map;
	
	std::vector<DoFIndex>  multInd;
	
// Loop the subset names
	t_ss_map::const_iterator ssIterEnd = m_mZeroBC.end ();
	for (t_ss_map::const_iterator ssIter = m_mZeroBC.begin (); ssIter != ssIterEnd; ++ssIter)
	{
		const FunctionGroup & fctGrp = ssIter->second;
		int si = ssIter->first;
		
	//	Loop the edges
		t_edge_iterator iterEnd = dd->end<Edge> (si);
		for (t_edge_iterator iter = dd->begin<Edge> (si); iter != iterEnd; iter++)
		{
			EdgeBase * pEdge = *iter;
			
		//	Loop dirichlet functions on this segment
			for(size_t i = 0; i < fctGrp.size (); ++i)
			{
			//	Get the multiindices
				if (dd->inner_dof_indices (pEdge, fctGrp[i], multInd) != 1)
					UG_THROW ("NedelecDirichletBC:"
						"More than one DoF per edge. Not the Nedelec-type-1 element?");
	
			//	Set the Dirichlet entry
				DoFRef (u, multInd[0]) = 0;
			}
		}
	}
}

/**
 * sets the dirichlet value in the solution for all dirichlet indices
 */
template <typename TDomain, typename TAlgebra>
void NedelecDirichletBC<TDomain, TAlgebra>::adjust_solution
(
	vector_type & u,
	ConstSmartPtr<DoFDistribution> dd, 
	number time
)
{
	extract_data ();
	
	// Loop boundary subsets (separately for every type of the BC specification)
	{
		typename std::map<int, std::vector<TConstBC*> >::const_iterator iter;
		for (iter = m_mConstBC.begin (); iter != m_mConstBC.end (); ++iter)
			adjust_solution<TConstBC> ((*iter).second, (*iter).first, u, dd, time);
	}
	{
		typename std::map<int, std::vector<TUserDataBC*> >::const_iterator iter;
		for (iter = m_mUserDataBC.begin (); iter != m_mUserDataBC.end (); ++iter)
			adjust_solution<TUserDataBC> ((*iter).second, (*iter).first, u, dd, time);
	}
	
	adjust_solution_implicit (u, dd, time);
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
