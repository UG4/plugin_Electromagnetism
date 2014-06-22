/*
 * Implementation of the functions and the commands for the Whitney-1 (Nedelec) elements.
 *
 * Created on: 04.06.2014
 * Author: D. Logashenko
 */

// basic ug4 headers
#include "common/common.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

/* Discretization's headers: */
#include "nedelec_local_ass.h"

namespace ug{
namespace Electromagnetism{

/// Helper class for computation of the flux in elements. (Helper for ComputeFlux.)
template <typename TGridFunc, typename TElem>
class ComputeElemFluxHelper
{
private:
	/// Computation of the flux through sides of an element.
	static number compute_elem_flux
	(
		TGridFunc * pGF, ///< grid function of the Nedelec-DoFs of the vector field
		size_t fct, ///< index of the function in the grid function
		TElem * pElem, ///< the element
		SubsetGroup & faceSSG ///< the surface (the low-dim. subsets)
	)
	{
		typedef typename TGridFunc::domain_type domain_type;
		static const int dim = domain_type::dim;
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
		typedef typename domain_traits<dim>::side_type side_type;
		
		number flux = 0;
		
	//	loop the sides
		const ISubsetHandler * pIsh = pGF->subset_handler().get ();
		for (size_t i = 0; i < (size_t) ref_elem_type::numSides; i++)
		{
		//	check whether that is a side on the face
			side_type * pSide = pIsh->grid()->get_side (pElem, i);
			if (! faceSSG.contains (pIsh->get_subset_index (pSide)))
				continue;
			
		//---- There is typically only one side on the surface, and not for every
		//---- element. So we put the following section into the loop although it
		//---- does not depend on the side. This allows to skip getting coordinates
		//---- and dofs for elements that have no sides on the surface.
		
			const ref_elem_type & rRefElem = Provider<ref_elem_type>::get ();
			
		//	get position accessor and corner coordinates of the element
			const typename domain_type::position_accessor_type & aaPos
				= pGF->domain()->position_accessor();
			typename domain_type::position_type corners [ref_elem_type::numCorners];
			for (size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
				corners[co] = aaPos [pElem->vertex (co)];
		
		//	get the dof values of the function in the element
			std::vector<DoFIndex> ind;
			if (pGF->dof_indices (pElem, fct, ind) != (size_t) ref_elem_type::numEdges)
				UG_THROW ("Grid function type mismatch.");
			number dofValues [ref_elem_type::numEdges];
			for (size_t dof = 0; dof < (size_t) ref_elem_type::numEdges; dof++)
				dofValues[dof] = DoFRef (*pGF, ind[dof]);
		
		//----
			
		//	get the integration point
			MathVector<dim> loc_center;
			loc_center = 0.0;
			size_t co;
			for (co = 0; co < pSide->num_vertices (); co++)
			{
				int elem_co = rRefElem.id (dim-1, i, 0, co);
				UG_ASSERT (elem_co >= 0, "Index mismatch.");
				loc_center += rRefElem.corner (elem_co);
			}
			loc_center /= co;
			
		//	compute the value of the function at the ip
			MathVector<dim> val;
			NedelecInterpolation<domain_type, dim>::value
				(pGF->domain().get (), pElem, corners, dofValues, &loc_center, 1, &val);
			
		//	get the normal to the side (its length is the area of the side)
			MathVector<dim> normal;
			SideNormal<ref_elem_type, dim> (normal, i, corners);
		
		//	update the flux
			flux += VecDot (val, normal);
		}
		
		return flux;
	}
	
public:
	/// Computation of the flux through the sides of all elements of one type.
	static number compute_flux
	(
		TGridFunc * pGF, ///< grid function of the Nedelec-DoFs of the vector field
		size_t fct, ///< index of the function in the grid function
		SubsetGroup & volSSG, ///< full-dim. subsets (adjacent to the surface) to indicate the negative direction
		SubsetGroup & faceSSG ///< the surface (the low-dim. subsets)
	)
	{
		typedef typename TGridFunc::template traits<TElem>::const_iterator t_elem_iterator;
	
		number flux = 0;
		
	//	loop all volume subsets
		for (size_t i = 0; i < volSSG.size (); i++)
		{
			int si = volSSG[i];
		//	loop all the elements
			t_elem_iterator elem_iter = pGF->template begin<TElem> (si);
			t_elem_iterator end_iter = pGF->template end<TElem> (si);
			for (; elem_iter != end_iter; ++elem_iter)
				flux += compute_elem_flux (pGF, fct, *elem_iter, faceSSG);
		}
		
		return flux;
	}
};

/// Helper class for computation of the flux in elements. (Helper for ComputeFlux.)
template <typename TGridFunc>
class ComputeElemFluxHelper<TGridFunc, RegularEdge>
{
public:
	/// Computation of the flux through the sides of all elements of one type.
	static number compute_flux
	(
		TGridFunc * pGF, ///< grid function of the Nedelec-DoFs of the vector field
		size_t fct, ///< index of the function in the grid function
		SubsetGroup & volSSG, ///< full-dim. subsets (adjacent to the surface) to indicate the negative direction
		SubsetGroup & faceSSG ///< the surface (the low-dim. subsets)
	)
	{
		UG_THROW ("ComputeFlux: No flux computations in 1d.");
	}
};

/// Helper class for the loop over all the elements in the computation of the flux. (Helper for ComputeFlux.)
template <typename TGridFunc>
class ComputeFluxHelper
{
	TGridFunc * m_pGF;
	size_t m_fct;
	SubsetGroup & m_volSSG;
	SubsetGroup & m_faceSSG;
	number & m_flux;
	
public:

///	class constructor
	ComputeFluxHelper
	(
		TGridFunc * pGF, ///< grid function of the Nedelec-DoFs of the vector field
		size_t fct, ///< index of the function in the grid function
		SubsetGroup & volSSG, ///< full-dim. subsets (adjacent to the surface) to indicate the negative direction
		SubsetGroup & faceSSG, ///< the surface (the low-dim. subsets)
		number & flux ///< where to save the flux
	)
	: m_pGF (pGF), m_fct (fct), m_volSSG (volSSG), m_faceSSG (faceSSG), m_flux (flux)
	{
		flux = 0;
	}
	
	template <typename TElem> void operator() (TElem &)
	{
		m_flux += ComputeElemFluxHelper<TGridFunc, TElem>::compute_flux (m_pGF, m_fct, m_volSSG, m_faceSSG);
	}
};

///	returns the flux through a given surface
template <typename TGridFunc>
number ComputeFlux
(
	TGridFunc * pGF,
	size_t fct,
	SubsetGroup & volSSG,
	SubsetGroup & faceSSG
)
{
	typedef typename TGridFunc::domain_type domain_type;
	static const int dim = domain_type::dim;
	typedef typename domain_traits<dim>::DimElemList ElemList;
	
	if (pGF->local_finite_element_id (fct).type () != LFEID::NEDELEC)
		UG_THROW ("ComputeFlux: Not a Nedelec-element-based grid function specified.");
	
	number flux;
	boost::mpl::for_each<ElemList>
		(ComputeFluxHelper<TGridFunc> (pGF, fct, volSSG, faceSSG, flux));
	return flux;
}

///	prints the flux through a given surface in the shell
template <typename TGridFunc>
void ComputeFlux
(
	SmartPtr<TGridFunc> spGF,
	const char * fct_names,
	const char * vol_subsets,
	const char * face_subsets
)
{
	FunctionGroup fctGrp;
	try
	{
		fctGrp = spGF->fct_grp_by_name (fct_names);
	}
	UG_CATCH_THROW ("ComputeFlux: Functions '" << fct_names <<
		"' not all contained in the edge approximation space.");
	if (fctGrp.size () != 1)
		UG_THROW ("ComputeFlux: Only one function component per call supported");
	
	SubsetGroup volSSG, faceSSG;
	std::vector<std::string> vssNames;
	
	TokenizeString (vol_subsets, vssNames);
	for (size_t i = 0; i < vssNames.size(); i++)
		RemoveWhitespaceFromString (vssNames [i]);
	volSSG.set_subset_handler (spGF->subset_handler ());
	volSSG.add (vssNames);
	
	TokenizeString (face_subsets, vssNames);
	for (size_t i = 0; i < vssNames.size(); i++)
		RemoveWhitespaceFromString (vssNames [i]);
	faceSSG.set_subset_handler (spGF->subset_handler ());
	faceSSG.add (vssNames);
	
	number flux = ComputeFlux (spGF.get(), fctGrp[0], volSSG, faceSSG);
	UG_LOG ("--> Flux: " << flux << '\n');
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
