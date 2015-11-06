/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
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

/*
 * Implementation of the functions and the commands for the eddy current simulations.
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

//---- Computation of the power of the electromagnetic field ----//

/// Helper class for the integration of the power
/**
 * This class implements the function computes the integral
 * \f[
 *  \int_\Omega \overline{\mathbf{J}_G} \mathbf{E} \, dx
 * \f]
 * over the subdomain covered by grid elements of one type.
 *
 * \tparam	TGridFunc	grid function type
 * \tparam	TElem		grid element type
 */
template <typename TGridFunc, typename TElem>
class CalcVolPowerElemHelperClass
{
private:
	typedef typename TGridFunc::domain_type domain_type;
	typedef typename domain_type::position_accessor_type position_accessor_type;
	typedef typename domain_type::position_type position_type;
	
	static const size_t numCorners = NedelecT1_LDisc<domain_type, TElem>::numCorners;
	static const size_t numEdges = NedelecT1_LDisc<domain_type, TElem>::numEdges;
	static const size_t maxEdges = NedelecT1_LDisc<domain_type, TElem>::maxNumEdges;
		
///	computes the integral over one element
	inline static void calc_elem_power
	(
		TGridFunc * pJGGF, ///< grid function of the Nedelec-DoFs of the generator current \f$\mathbf{J}_G\f$
		size_t JG_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{J}_G\f$
		TGridFunc * pEGF, ///< grid function of the Nedelec-DoFs of the electric field \f$\mathbf{E}\f$
		size_t E_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{E}\f$
		position_accessor_type & aaPos, ///< position accessor
		TElem * pElem, ///< the element
		number pow[] ///< to add the integral (Re and Im parts)
	)
	{
	//	Get coordinates of the corners
		position_type corners [numCorners];
		for (size_t co = 0; co < numCorners; co++)
			corners [co] = aaPos [pElem->vertex (co)];
		
	//	Assemble the local mass matrix
		number M [maxEdges] [maxEdges];
		NedelecT1_LDisc<domain_type, TElem>::local_mass
			(pEGF->domain().get (), pElem, corners, M);
		
	//	Get the components of the grid functions
		number Re_JG_val [numEdges], Im_JG_val [numEdges];
		number Re_E_val [numEdges], Im_E_val [numEdges];
		std::vector<DoFIndex> ind;
		
		if (pJGGF->dof_indices (pElem, JG_fct[0], ind) != numEdges)
			UG_THROW ("J_G grid function type mismatch.");
		for (size_t dof = 0; dof < numEdges; dof++)
			Re_JG_val[dof] = DoFRef (*pJGGF, ind[dof]);
		
		if (pJGGF->dof_indices (pElem, JG_fct[1], ind) != numEdges)
			UG_THROW ("J_G grid function type mismatch.");
		for (size_t dof = 0; dof < numEdges; dof++)
			Im_JG_val[dof] = DoFRef (*pJGGF, ind[dof]);
		
		if (pEGF->dof_indices (pElem, E_fct[0], ind) != numEdges)
			UG_THROW ("E grid function type mismatch.");
		for (size_t dof = 0; dof < numEdges; dof++)
			Re_E_val[dof] = DoFRef (*pEGF, ind[dof]);
		
		if (pEGF->dof_indices (pElem, E_fct[1], ind) != numEdges)
			UG_THROW ("E grid function type mismatch.");
		for (size_t dof = 0; dof < numEdges; dof++)
			Im_E_val[dof] = DoFRef (*pEGF, ind[dof]);
		
	//	Compute the integral
		number Re_ME [numEdges], Im_ME [numEdges];
		
		memset (Re_ME, 0, numEdges * sizeof (number));
		memset (Im_ME, 0, numEdges * sizeof (number));
		for (size_t i = 0; i < numEdges; i++)
			for (size_t j = 0; j < numEdges; j++)
			{
				Re_ME [i] +=  M [i] [j] * Re_E_val [j];
				Im_ME [i] +=  M [i] [j] * Im_E_val [j];
			}
	
	//	a) the real part
		number Re_pow = 0;
		for (size_t i = 0; i < numEdges; i++)
			Re_pow += Re_JG_val[i] * Re_ME[i] + Im_JG_val[i] * Im_ME[i];
	//	b) the imaginary part
		number Im_pow = 0;
		for (size_t i = 0; i < numEdges; i++)
			Im_pow += Re_JG_val[i] * Im_ME[i] - Re_ME[i] * Im_JG_val[i];
	
	//	Done
		pow[0] += Re_pow; pow[1] += Im_pow;
	
	//-- For debugging only
	//	UG_LOG ("==> subset " << pEGF->domain()->subset_handler()->get_subset_index (pElem)
	//		<< ": Re = " << Re_pow << ", Im = " << Im_pow << "\n");
	//--
	}

public:

///	Computes the integral for all elements of one type
	static void calc_power
	(
		TGridFunc * pJGGF, ///< grid function of the Nedelec-DoFs of the generator current \f$\mathbf{J}_G\f$
		size_t JG_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{J}_G\f$
		SubsetGroup & JG_ssg, ///< (full-dim.) subsets where \f$\mathbf{J}_G\f$ is defined (non-zero and in the kernel)
		TGridFunc * pEGF, ///< grid function of the Nedelec-DoFs of the electric field \f$\mathbf{E}\f$
		size_t E_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{E}\f$
		number pow[] ///< to add the integral (Re and Im parts)
	)
	{
		typedef typename TGridFunc::template traits<TElem>::const_iterator t_elem_iterator;
		
		position_accessor_type & aaPos = pEGF->domain()->position_accessor();
		
	//	loop the subdomains
		for (size_t i = 0; i < JG_ssg.size (); i++)
		{
			int si = JG_ssg [i];
			
		//	loop all the elements of the given type in the subset
			t_elem_iterator elem_iter = pEGF->template begin<TElem> (si);
			t_elem_iterator end_iter = pEGF->template end<TElem> (si);
			for (; elem_iter != end_iter; ++elem_iter)
				CalcVolPowerElemHelperClass<TGridFunc, TElem>::calc_elem_power
					(pJGGF, JG_fct, pEGF, E_fct, aaPos, *elem_iter, pow);
		}
	}
};

/// Helper class for the computation of the power of the electromagnetic field
template <typename TGridFunc>
class CalcVolPowerHelperClass
{
	TGridFunc * m_pJGGF;
	size_t m_JG_fct[2];
	SubsetGroup & m_JG_ssg;
	TGridFunc * m_pEGF;
	size_t m_E_fct[2];
	number * m_pow;
	
public:
	
///	class constructor
	CalcVolPowerHelperClass
	(
		TGridFunc * pJGGF, ///< grid function of the Nedelec-DoFs of the generator current \f$\mathbf{J}_G\f$
		size_t JG_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{J}_G\f$
		SubsetGroup & JG_ssg, ///< (full-dim.) subsets where \f$\mathbf{J}_G\f$ is defined (non-zero and in the kernel)
		TGridFunc * pEGF, ///< grid function of the Nedelec-DoFs of the electric field \f$\mathbf{E}\f$
		size_t E_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{E}\f$
		number pow[] ///< to add the integral
	)
	: m_pJGGF (pJGGF), m_JG_ssg (JG_ssg), m_pEGF (pEGF), m_pow (pow)
	{
		m_JG_fct[0] = JG_fct[0]; m_JG_fct[1] = JG_fct[1];
		m_E_fct[0] = E_fct[0]; m_E_fct[1] = E_fct[1];
		m_pow[0] = m_pow[1] = 0;
	}
	
	template <typename TElem> void operator() (TElem &)
	{
		CalcVolPowerElemHelperClass<TGridFunc, TElem>::calc_power
			(m_pJGGF, m_JG_fct, m_JG_ssg, m_pEGF, m_E_fct, m_pow);
	}
};

/// Computes the power of the electromagnetic field (up to the contribution of the boundary)
template <typename TGridFunc>
void calc_power
(
	TGridFunc * pJGGF, ///< grid function of the Nedelec-DoFs of the generator current \f$\mathbf{J}_G\f$
	size_t JG_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{J}_G\f$
	SubsetGroup & JG_ssg, ///< (full-dim.) subsets where \f$\mathbf{J}_G\f$ is defined (non-zero and in the kernel)
	TGridFunc * pEGF, ///< grid function of the Nedelec-DoFs of the electric field \f$\mathbf{E}\f$
	size_t E_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{E}\f$
	number pow[] ///< to add the integral
)
{
	typedef typename TGridFunc::domain_type domain_type;
	static const int dim = domain_type::dim;
	typedef typename domain_traits<dim>::DimElemList ElemList;
	
	boost::mpl::for_each<ElemList>
		(CalcVolPowerHelperClass<TGridFunc> (pJGGF, JG_fct, JG_ssg, pEGF, E_fct, pow));
	
	pow [0] /= -2;
	pow [1] /= -2;
}

/// Prints the (complex-valued) power of the electromagnetic field
template <typename TGridFunc>
void CalcPower
(
	SmartPtr<TGridFunc> spJGGF, ///< [in] grid function with the generator current
    const char* JG_cmps, ///< [in] names of the components of the grid function (for Re and Im)
	const char* JG_ss, ///< (full-dim.) subsets where \f$\mathbf{J}_G\f$ is defined (non-zero and in the kernel)
	SmartPtr<TGridFunc> spEGF, ///< [in] grid function with the electric field
    const char* E_cmps ///< [in] names of the components of the grid function (for Re and Im)
)
{
//	Get functions by names
	size_t JG_fcts [2];
	{
		std::string fctString = std::string (JG_cmps);
		std::vector<std::string> tokens;
		TokenizeString (fctString, tokens, ',');
		if (tokens.size () != 2)
			UG_THROW("EddyCurrent: Needed 2 components "
						"in symbolic function names (for Re and Im) for JG, "
							"but given: " << JG_cmps);
		for (size_t i = 0; i < tokens.size (); i++)
			RemoveWhitespaceFromString (tokens [i]);
	
		for (int i = 0; i < 2; i++)
		try
		{
			JG_fcts [i] = spJGGF->fct_id_by_name (tokens[i].c_str ());
		}
		UG_CATCH_THROW ("EddyCurrent: Cannot find symbolic function "
						"component for the name '" << tokens[i] << "' (in JG).");
	}
	size_t E_fcts [2];
	{
		std::string fctString = std::string (E_cmps);
		std::vector<std::string> tokens;
		TokenizeString (fctString, tokens, ',');
		if (tokens.size () != 2)
			UG_THROW("EddyCurrent: Needed 2 components "
						"in symbolic function names (for Re and Im) for E, "
							"but given: " << E_cmps);
		for (size_t i = 0; i < tokens.size (); i++)
			RemoveWhitespaceFromString (tokens [i]);
	
		for (int i = 0; i < 2; i++)
		try
		{
			E_fcts [i] = spEGF->fct_id_by_name (tokens[i].c_str ());
		}
		UG_CATCH_THROW ("EddyCurrent: Cannot find symbolic function "
						"component for the name '" << tokens[i] << "' (in E).");
	}
	
	SubsetGroup JG_ssg;
	{
		std::string ssString = std::string (JG_ss);
		std::vector<std::string> tokens;
		TokenizeString (ssString, tokens, ',');
		for (size_t i = 0; i < tokens.size (); i++)
			RemoveWhitespaceFromString (tokens [i]);
		JG_ssg.set_subset_handler (spJGGF->subset_handler ());
		JG_ssg.add (tokens);
	}
	
//	Compute the power
	number power [2];
	calc_power (spJGGF.get (), JG_fcts, JG_ssg, spEGF.get (), E_fcts, power);
	
//	Print the result
	UG_LOG ("--> Power of the electromagnetic field: "
		<< power [0] << " + " << power [1] << " I.\n");
	UG_LOG ("Note: Contribution of the boundary is not taken into account.\n");
}

//---- Computation of magnetic flux through a cross-section ----//

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
