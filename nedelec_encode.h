/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/**
 * Computation of the vector of the Nedelec-element DoFs from a given
 * vector field.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_ENCODE__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_ENCODE__

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace Electromagnetism{

/**
 * Computes the value for the DoF associated with a given edge. Uses the
 * trapezoid quadrature rule.
 */
template <typename TDomain>
inline number ComputeNedelecDoF
(
	UserData<MathVector<TDomain::dim>, TDomain::dim> * function, ///< the function to encode
	Edge * edge, ///< the edge to compute the DoF for
	const typename TDomain::position_accessor_type& aaPos, ///< coordinates of the vertices
	const int si, ///< subset index
	number time ///< the time argument
)
{
	typedef typename TDomain::position_type position_type;
	
	number dofValue;
	MathVector<TDomain::dim> func_value;
	
//	Get the vector pointing along the edge from its beginning to its end:
	position_type edgeVector = aaPos [edge->vertex(1)];
	edgeVector -= aaPos [edge->vertex(0)];

//	Compute the quadrature

	(*function) (func_value, aaPos [edge->vertex(0)], time, si);
	dofValue = func_value * edgeVector;
	
	(*function) (func_value, aaPos [edge->vertex(1)], time, si);
	dofValue += func_value * edgeVector;
	
	return dofValue / 2;
}

/**
 * Computes the values of the DoFs of the Whitney-1 (Nedelec-type-1) basis
 * functions for a given vector function.
 */
template <typename TGridFunction>
void ComputeNedelecDoFs
(
	SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spFunction, ///< the function to encode
	SmartPtr<TGridFunction> spGridFct, ///< the grid function to store the result
	size_t fct, ///< index of the function component in the grid function
	const SubsetGroup& ssGrp, ///< subset group where to encode
	number time ///< the time argument (for the function to encode)
)
{
//	Domain type, position_type and iterators
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename TGridFunction::template traits<Edge>::const_iterator t_edge_iterator;

//	Get position accessor
	const typename domain_type::position_accessor_type& aaPos = spGridFct->domain()->position_accessor();

//	Multiindex to access the components
	std::vector<DoFIndex> ind;

	UG_ASSERT (spFunction.valid (), "Invalid smart pointer to the user data.");
	
//	Loop over subsets
	for (size_t i = 0; i < ssGrp.size (); i++)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if (! spGridFct->is_def_in_subset (fct, si))
		{
			UG_LOG ("ComputeNedelecDoFs warning: Function " << fct << " undefined in subdom " << si);
			continue;
		}
	
	//	Loop the edges
		t_edge_iterator iter = spGridFct->template begin<Edge> (si);
		t_edge_iterator iterEnd = spGridFct->template end<Edge> (si);
		for (; iter != iterEnd; iter++)
		{
			Edge * pEdge = *iter;
			
		//	Get the multiindices
			if (spGridFct->inner_dof_indices (pEdge, fct, ind) != 1)
				UG_THROW ("ComputeNedelecDoFs:"
					" More than one DoF per edge. Not the Nedelec-type-1 element?");
		
		//	Compute the value for the edge and write it into the grid function
			DoFRef ((* spGridFct), ind[0])
				= ComputeNedelecDoF<domain_type> (spFunction.get(), pEdge, aaPos, si, time);
		}
	}
}

/**
 * Computes the values of the DoFs of the Whitney-1 (Nedelec-type-1) basis
 * functions for a given vector function. This function parses the string
 * specifications of the function components and the subsets.
 */
///\{
template <typename TGridFunction>
void ComputeNedelecDoFs
(
	SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spFunction, ///< the function to encode
	SmartPtr<TGridFunction> spGridFct, ///< the grid function to store the result
	const char* cmp, ///< name of the component
    const char* subsets, ///< names of the subsets (may be NULL meaning "all")
    number time ///< the time argument
)
{
//	Get function id by name
	const size_t fct = spGridFct->fct_id_by_name (cmp);

//	Check that function exists
	if (fct > spGridFct->num_fct())
		UG_THROW("ComputeNedelecDoFs: Component with name '" << cmp << "' not found.");
	
//	Check the basis type
	if (spGridFct->local_finite_element_id(fct).type () != LFEID::NEDELEC)
		UG_THROW ("ComputeNedelecDoFs: Illegal grid function. It should be based on the Nedelec element.");
	
//	Create the subset group
	SubsetGroup ssGrp (spGridFct->domain()->subset_handler ());
	if (subsets != NULL)
		ssGrp.add (TokenizeString (subsets));
	else //	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
	
//	Compute the DoF values:
	ComputeNedelecDoFs (spFunction, spGridFct, fct, ssGrp, time);
}

template <typename TGridFunction>
void ComputeNedelecDoFs
(
	SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spFunction, ///< the function to encode
	SmartPtr<TGridFunction> spGridFct, ///< the grid function to store the result
	const char* cmp, ///< name of the component
    number time ///< the time argument
)
{
	ComputeNedelecDoFs (spFunction, spGridFct, cmp, NULL, time);
}

template <typename TGridFunction>
void ComputeNedelecDoFs
(
	SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spFunction, ///< the function to encode
	SmartPtr<TGridFunction> spGridFct, ///< the grid function to store the result
	const char* cmp, ///< name of the component
    const char* subsets ///< names of the subsets (may be NULL meaning "all")
)
{
	ComputeNedelecDoFs (spFunction, spGridFct, cmp, subsets, 0.0);
}

template <typename TGridFunction>
void ComputeNedelecDoFs
(
	SmartPtr<UserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > spFunction, ///< the function to encode
	SmartPtr<TGridFunction> spGridFct, ///< the grid function to store the result
	const char* cmp ///< name of the component
)
{
	ComputeNedelecDoFs (spFunction, spGridFct, cmp, NULL, 0.0);
}

///\}

#ifdef UG_FOR_LUA

/**
 * Computes the values of the DoFs of the Whitney-1 (Nedelec-type-1) basis
 * functions for a given vector function. This function can be called from
 * lua scripts.
 */
template <typename TGridFunction>
void ComputeNedelecDoFs
(
	const char* LuaFunction, ///< the simpolic lua function
    SmartPtr<TGridFunction> spGridFct, ///< the grid function to initialize
    const char* cmp, ///< name of the component of the grid function
	const char* subsets, ///< subsets where to compute the DoFs
	number time ///< time argument
)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<MathVector<dim>, dim> > spFunction
				= LuaUserDataFactory<MathVector<dim>, dim>::create (LuaFunction);
	ComputeNedelecDoFs (spFunction, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
void ComputeNedelecDoFs
(
	const char* LuaFunction, ///< the simpolic lua function
    SmartPtr<TGridFunction> spGridFct, ///< the grid function to initialize
    const char* cmp, ///< name of the component of the grid function
	number time ///< time argument
)
{
	ComputeNedelecDoFs (LuaFunction, spGridFct, cmp, NULL, time);
}
template <typename TGridFunction>
void ComputeNedelecDoFs
(
	const char* LuaFunction, ///< the simpolic lua function
    SmartPtr<TGridFunction> spGridFct, ///< the grid function to initialize
    const char* cmp, ///< name of the component of the grid function
	const char* subsets ///< subsets where to compute the DoFs
)
{
	ComputeNedelecDoFs (LuaFunction, spGridFct, cmp, subsets, 0);
}
template <typename TGridFunction>
void ComputeNedelecDoFs
(
	const char* LuaFunction, ///< the simpolic lua function
    SmartPtr<TGridFunction> spGridFct, ///< the grid function to initialize
    const char* cmp ///< name of the component of the grid function
)
{
	ComputeNedelecDoFs (LuaFunction, spGridFct, cmp, NULL, 0);
}

#endif // UG_FOR_LUA

/**
 * Computes the Whitney representation of the gradient of a vertex-centered
 * scalar grid function.
 */
template <typename TPotGridFunc, typename TGridFunc>
void NedelecGradPotential
(
	SmartPtr<TPotGridFunc> spPotGF, ///< the potential to compute the gradient for
	SmartPtr<TGridFunc> spGF, ///< to save the Whitney form of the gradient
	size_t fct ///< index of the function in spGF
)
{
	typedef typename TGridFunc::template traits<Edge>::const_iterator t_edge_iterator;
	
//	Check the basis type
	if (spGF->local_finite_element_id(fct).type () != LFEID::NEDELEC)
		UG_THROW ("NedelecGradPotential: Illegal grid function. It should be based on the Nedelec element.");
	
//	Arrays for the indices in the vectors:
	std::vector<size_t> vVertInd (1);
	std::vector<DoFIndex> vEdgeInd (1);
	
//	Loop the edges
	t_edge_iterator iterEnd = spGF->template end<Edge> ();
	for (t_edge_iterator iter = spGF->template begin<Edge> (); iter != iterEnd; iter++)
	{
		number corner_val [2];
		Edge * pEdge = *iter;
		
	//	Get the potential the ends of the edge:
		for (size_t i = 0; i < 2; i++)
		{
		//	Get the multiindices
			if (spPotGF->inner_algebra_indices (pEdge->vertex(i), vVertInd) != 1)
				UG_THROW ("NedelecGradPotential:"
					" More than one DoF per vertex for the potential.");

		//	Set the Dirichlet entry
			corner_val [i] = (* spPotGF) [vVertInd[0]];
		}
		
	//	Compute the gradient:
		if (spGF->inner_dof_indices (pEdge, fct, vEdgeInd) != 1)
			UG_THROW ("NedelecGradPotential:"
				" More than one DoF per edge. Not the Nedelec-Type-1 element?");
		DoFRef ((* spGF), vEdgeInd[0])
			-= corner_val [1] - corner_val [0];
	}
}

/**
 * Computes the Whitney representation of the gradient of a vertex-centered
 * scalar grid function. (This is the version with the symbolic function name.)
 */
template <typename TPotGridFunc, typename TGridFunc>
void NedelecGradPotential
(
	SmartPtr<TPotGridFunc> spPotGF, ///< the potential to compute the gradient for
	SmartPtr<TGridFunc> spGF, ///< to save the Whitney form of the gradient
	const char* cmp ///< name of the component
)
{
//	Get function id by name
	const size_t fct = spGF->fct_id_by_name (cmp);

//	Check that function exists
	if (fct > spGF->num_fct())
		UG_THROW("NedelecGradPotential: Component with name '" << cmp << "' not found.");
	
//	Distribute the potential
	NedelecGradPotential (spPotGF, spGF, fct);
}

/**
 * Sets a scalar potential field to a given value at vertices of a given
 * subset, and to 0 for all other subsets.
 */
template <typename TPotGridFunc>
void SetSubsetVertVal
(
	SmartPtr<TPotGridFunc> spPotGF, ///< the potential to compute the gradient for
	const SubsetGroup& ssGrp, ///< subset group where to set
	number value ///< value to set in the subset group
)
{
	typedef typename TPotGridFunc::template traits<Vertex>::const_iterator t_vert_iterator;
	
//	Array for the indices in the grid function:
	std::vector<size_t> vVertInd (1);
	
//	Set the function to 0
	spPotGF->set (0.0);
	
//	Loop over subsets
	for (size_t i = 0; i < ssGrp.size (); i++)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	Loop the vertices
		t_vert_iterator iter = spPotGF->template begin<Vertex> (si);
		t_vert_iterator iterEnd = spPotGF->template end<Vertex> (si);
		for (; iter != iterEnd; iter++)
		{
			Vertex * pVertex = *iter;
			
		//	Get the multiindices
			if (spPotGF->inner_algebra_indices (pVertex, vVertInd) != 1)
				UG_THROW ("SubsetPotential:"
					" More than one DoF per vertex for the potential.");
			
		//	Set the value
			(* spPotGF) [vVertInd[0]] = value;
		}
	}
}

/**
 * Sets a scalar potential field to a given value at vertices of a given
 * subset, and to 0 for all other subsets.
 */
template <typename TPotGridFunc>
void SetSubsetVertVal
(
	SmartPtr<TPotGridFunc> spPotGF, ///< the potential to compute the gradient for
    const char* subsets, ///< names of the subsets (may be NULL meaning "all")
	number value ///< value to set in the subset group
)
{
//	Create the subset group
	SubsetGroup ssGrp (spPotGF->domain()->subset_handler ());
	if (subsets != NULL)
		ssGrp.add (TokenizeString (subsets));
	else //	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
	
//	Set the potential
	SetSubsetVertVal (spPotGF, ssGrp, value);
}

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_ENCODE__

/* End of File */
