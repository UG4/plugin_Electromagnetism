/**
 * Computation of the vector of the Nedelec-element DoFs from a given
 * vector field.
 *
 * Created on: 01.02.2013
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_ENCODE__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_ENCODE__

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_disc/common/subset_group.h"
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
	typedef typename domain_type::position_type position_type;
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
					"More than one DoF per edge. Not the Nedelec-type-1 element?");
		
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
	if (spGridFct->local_finite_element_id (fct).type() != LFEID::NEDELEC)
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

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_ENCODE__

/* End of File */
