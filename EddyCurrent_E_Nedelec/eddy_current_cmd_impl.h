/*
 * Implementation of the functions and the commands for the eddy current simulations.
 *
 * Created on: 26.05.2014
 * Author: D. Logashenko
 */

// basic ug4 headers
#include "common/common.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

/* Discretization's headers: */
#include "../nedelec_local_ass.h"

namespace ug{
namespace Electromagnetism{

/// Computation of the magnetic flux through windings of a coil
/**
 * This function computes the magnetic flux through windings of a cylindric coil
 * and returns the electomotive force due to the varying magnetic flux.
 * The flux is a complex number (as well as the electric field that is used as
 * the base of the computation). The function returns the total area of all the
 * windings.
 */
template <typename TGridFunc>
number calc_EMF
(
	const TGridFunc * gfE, ///< [in] grid function with the electric field
	const size_t fct [2], ///< [in] function components (Re and Im)
	const SubsetGroup & ss_grp, ///< [in] subsets of the kerner of the coil
	const MathVector<TGridFunc::dim> & Normal, ///< [in] normal to the planes of the windings
	const typename TGridFunc::domain_type::position_type & base_pnt, ///< [in] point on the plane of the 1st winding
	const size_t n_pnt, ///< [in] number of the windings
	const typename TGridFunc::domain_type::position_type & d_pnt, ///< [in] thickness of the windings
	number emf [2] ///< [out] the computed electromotive force (Re and Im)
)
{
	typedef typename TGridFunc::domain_type domain_type;
	typedef typename domain_type::position_type position_type;
	static const int dim = TGridFunc::dim;
	typedef typename domain_traits<dim>::grid_base_object elem_type;
	typedef typename TGridFunc::template traits<elem_type>::const_iterator t_elem_iterator;
	static const size_t maxCorners = domain_traits<dim>::MaxNumVerticesOfElem;
	static const size_t maxEdges = (size_t) element_list_traits<typename domain_traits<dim>::DimElemList>::maxEdges;
	
	position_type corners [maxCorners];
	std::vector<DoFIndex> ind;
	number dofValues [2] [maxEdges];
	
//	check the function space of the grid function
	for (size_t part = 0; part < 2; part++)
		if (gfE->local_finite_element_id(fct[part]).type () != LFEID::NEDELEC)
			UG_THROW ("EddyCurrent: Illegal discretization space of the grid func. (must be NEDELEC)");
	
//	check whether the functions are defined in the specified subsets
	for (size_t ssgi = 0; ssgi < ss_grp.size (); ssgi++)
		for (size_t i = 0; i < 2; i++)
			if (! gfE->is_def_in_subset (i, ss_grp [ssgi]))
				UG_THROW ("EddyCurrent: Grid function is undefined in some of the subsets.");
	
//	normalize the normal
	number normal_len = VecLength (Normal);
	if (normal_len < 1e-16)
		UG_THROW ("EddyCurrent: Illegal normal vector specified (almost zero).");
	MathVector<dim> normal (Normal);
	normal *= 1 / normal_len;

//	get the domain
	const typename TGridFunc::domain_type * domain = gfE->domain().get ();
	
//	get position accessor
	const typename domain_type::position_accessor_type& aaPos = domain->position_accessor ();

//	compute the flux for every subset
	number total_area = 0;
	emf[0] = emf[1] = 0;
	for (size_t ssgi = 0; ssgi < ss_grp.size (); ssgi++)
	{
		int si = ss_grp [ssgi];
		
	//	Loop the (full-dim) elements
		t_elem_iterator iterEnd = gfE->template end<elem_type> (si);
		for (t_elem_iterator iter = gfE->template begin<elem_type> (si); iter != iterEnd; iter++)
		{
			elem_type * elem = *iter;
			
		//	get the corner coordinates
			for (size_t co = 0; co < elem->num_vertices (); co++)
				corners [co] = aaPos [elem->vertex (co)];
				
		//	get the DoFs
			for (size_t part = 0; part < 2; part++)
			{
				if (gfE->dof_indices (elem, fct[part], ind) > maxEdges)
					UG_THROW ("EddyCurrent: Illegal number of DoFs per element.");
				for (size_t sh = 0; sh < ind.size (); sh++)
					dofValues [part] [sh] = DoFRef (*gfE, ind [sh]);
			}
			
		//	loop over the windings
			position_type pnt = base_pnt;
			for (size_t i = 0; i < n_pnt; i++, pnt += d_pnt)
			{
				number elem_flux [2];
				number elem_area;
				if ((elem_area = NedelecInterpolation<domain_type, dim>::curl_flux
					(domain, (GridObject *) elem, corners, dofValues[0], normal, pnt, elem_flux[0])) != 0.0)
				{
					NedelecInterpolation<domain_type, dim>::curl_flux
						(domain, (GridObject *) elem, corners, dofValues[1], normal, pnt, elem_flux[1]);
					emf[0] += elem_flux[0];
					emf[1] += elem_flux[1];
					total_area += elem_area;
				}
			}
		}
	}
	
//	done
	return total_area;
};

/// Prints of the electromotive force due to the time variation of the magnetic flux through windings of a coil
template <typename TGridFunc>
void CalcEMF
(
	SmartPtr<TGridFunc> spGF, ///< [in] grid function with the electric field
    const char* cmps, ///< [in] names of the components of the grid function (for Re and Im)
	const char* subsets, ///< [in] subsets where to compute the flux
	const std::vector<number>& Normal, ///< [in] normal to the planes of the windings
	const std::vector<number>& base_pnt, ///< [in] point on the plane of the 1st winding
	const size_t n_pnt, ///< [in] number of the windings
	const std::vector<number>& d_pnt ///< [in] thickness of the windings
)
{
	typedef typename TGridFunc::domain_type domain_type;
	typedef typename domain_type::position_type position_type;
	static const int dim = TGridFunc::dim;
	
//	check the sizes of the arrays and copy the values
	if (Normal.size () != (size_t) dim || base_pnt.size () != (size_t) dim || d_pnt.size () != (size_t) dim)
		UG_THROW ("EddyCurrent: Illegal sizes of the vectors in the arguments.");
	MathVector<dim> normal_vec;
	position_type base_pnt_coord, d_pnt_coord;
	for (int i = 0; i < dim; i++)
	{
		normal_vec[i] = Normal[i];
		base_pnt_coord[i] = base_pnt[i];
		d_pnt_coord[i] = d_pnt[i];
	}
	
//	create the subset group
	if (subsets == NULL)
		UG_THROW ("EddyCurrent: No subsets specified.");
	SubsetGroup ssGrp (spGF->domain()->subset_handler ());
	{
		std::string ssString = std::string (subsets);
		std::vector<std::string> tokens;
		TokenizeString (ssString, tokens, ',');
		for (size_t i = 0; i < tokens.size (); i++)
			RemoveWhitespaceFromString (tokens [i]);
		ssGrp.add (tokens);
	}
	
//	get functions by names
	size_t fcts [2];
	{
		std::string fctString = std::string (cmps);
		std::vector<std::string> tokens;
		TokenizeString (fctString, tokens, ',');
		if (tokens.size () != 2)
			UG_THROW("EddyCurrent: Needed 2 components "
						"in symbolic function names (for Re and Im), "
							"but given: " << cmps);
		for (size_t i = 0; i < tokens.size (); i++)
			RemoveWhitespaceFromString (tokens [i]);
	
		for (int i = 0; i < 2; i++)
		try
		{
			fcts [i] = spGF->fct_id_by_name (tokens[i].c_str ());
		}
		UG_CATCH_THROW ("EddyCurrent: Cannot find symbolic function "
						"component for the name '" << tokens[i] << "'.");
	}
	
//	compute the flux and the area
	number emf [2], area;
	area = calc_EMF (spGF.get(), fcts, ssGrp, normal_vec,
		base_pnt_coord, n_pnt, d_pnt_coord, emf);
	
//	print the result
	UG_LOG ("--> Electromotive force in the coil with " << n_pnt << " windings (total area "
		<< area << "): " << emf[0] << " + " << emf[1] << " I.\n");
};

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
