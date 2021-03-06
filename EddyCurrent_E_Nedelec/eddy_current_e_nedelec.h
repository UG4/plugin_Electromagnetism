/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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
 * FE-discretization of the time-harmonic E-based formulation of the eddy
 * current model.
 */

#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_E_NEDELEC__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_E_NEDELEC__

// basic ug4 headers
#include "common/common.h"

// library-specific headers
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "lib_disc/function_spaces/grid_function.h"

/* Discretization's headers: */
#include "eddy_current_traits.h"
#include "../em_material.h"
#include "../nedelec_local_ass.h"

namespace ug{
namespace Electromagnetism{

/// \ingroup lib_disc_elem_disc
/// @{

/// FE-discretization of the time-harmonic E-based formulation of the eddy current model.
/**
 * This class implements the local FE-discretization of the time-harmonic E-based
 * formulation of the eddy current model using the linear Nedelec-Type-1 (Whitney-1)
 * elements.
 * The model is formulated as follows:
 * \f{eqnarray*}{
 *  \mathbf{rot} \mu^{-1} \mathbf{rot} \mathbf{E} + i \omega \sigma \mathbf{E} = - i \omega \mathbf{J}_G,
 * \f}
 * where
 * <ul>
 * <li> \f$ \mathbf{E} \f$		(unknown, complex-valued scalar) electric field
 * <li> \f$ \mu \f$				(given) magnetic permeability (a real number)
 * <li> \f$ \sigma \f$			(given) conductivity (a real number)
 * <li> \f$ \omega \f$			(given) frequency (a real number)
 * <li> \f$ \mathbf{J}_G \f$	(given) generator current (non-zero only in the subdomain with \f$ \sigma = 0 \f$)
 * <li> \f$ i = \sqrt{-1} \f$
 * </ul>
 *
 * \remark
 * Note that every Nedelec shape function requires a scalar (complex) dof value.
 * These functions are thus vector-valued for the scalar dofs, each dof storing
 * two doubles (one for the real and one for the imaginary part of the value).
 *
 * References:
 * <ul>
 * <li> O. Sterz. Modellierung und Numerik zeitharmonischer Wirbelstromprobleme in 3D. PhD thesis, Heidelberg University, 2003.
 * <li> A. Bossavit. Computational Electromagnetism. Academic Press (Boston), 1998
 * </ul>
 *
 * \tparam	TDomain		Domain type
 * \tparam	TAlgebra	Algebra type
 */
template <typename TDomain, typename TAlgebra>
class EddyCurrent_E_Nedelec
	: public IElemDisc<TDomain>, public EddyCurrentTraits
{
private:
///	base class type
	typedef IElemDisc<TDomain> base_type;

///	own type
	typedef EddyCurrent_E_Nedelec<TDomain,TAlgebra> this_type;

/// type of grid functions (used for the sources)
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;

///	domain type
	typedef typename base_type::domain_type domain_type;

///	position type
	typedef typename base_type::position_type position_type;
	
///	world dimension
	static const int dim = base_type::dim;

/// max. number of the edges of the full-dimensional elements in the domain
	static const size_t maxNumEdges = element_list_traits<typename domain_traits<dim>::DimElemList>::maxEdges;

public:
///	class constructor
	EddyCurrent_E_Nedelec
	(
		const char * functions,
		ConstSmartPtr<EMaterial<domain_type> > spSubsetData,
		number frequency
	);

public:
/// adds a generator current item \f$ \mathbf{J}_{G,h} \f$ to the discretization
	void set_generator_current
	(
		SmartPtr<TGridFunction> spgfJG, ///< grid function with the data
		const char * cmp, ///< grid function components
		const char * ss_names = NULL ///< names of the subsets where the current is defined (NULL for "everywhere")
	);
/// adds a generator current item \f$ \mathbf{J}_{G,h} \f$ to the discretization
	void set_generator_current
	(
		SmartPtr<TGridFunction> spgfJG, ///< grid function with the data
		const char * cmp ///< grid function components
	)
	{
		set_generator_current (spgfJG, cmp, NULL);
	}

//---- Local discretization interface: ----
private:
	
///	check type of the grid and the trial space
	virtual void prepare_setting
	(
		const std::vector<LFEID> & vLfeID,
		bool bNonRegular
	);

//---- Assembling functions: ----
	
	template <typename TElem>
	void prepare_element_loop(ReferenceObjectID roid, int si);

	template <typename TElem>
	void prepare_element(const LocalVector& u, GridObject* elem, ReferenceObjectID roid, const position_type vCornerCoords[]);

	template <typename TElem>
	void finish_element_loop();

	template <typename TElem>
	void ass_JA_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_JM_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dA_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dM_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_rhs_elem(LocalVector& d, GridObject* elem, const position_type vCornerCoords[]);

//---- Registration of the template functions: ----
private:
	
	void register_all_loc_discr_funcs();

	struct RegisterLocalDiscr {
			RegisterLocalDiscr(this_type * pThis) : m_pThis(pThis){}
			this_type * m_pThis;
			template< typename TElem > void operator() (TElem &)
			{m_pThis->register_loc_discr_func<TElem> ();}
	};

	template <typename TElem>
	void register_loc_discr_func ();

//---- Auxiliary functions: ----
private:

/// composes the stiffness matrix of the stationary problem
	template<size_t numEdges>
	void ass_elem_stiffness
	(
		number perm, ///< the magnetic permeability
		number cond, ///< the electric conductivity
		number S [2][numEdges] [2][numEdges] ///< for the composed matrix
	);

//---- Physical parameters of the problem: ----
private:
/// frequency \f$\omega\f$ for the discretization
	number m_omega;

/// parameters of the materials in the domain
	ConstSmartPtr<EMaterial<domain_type> > m_spSubsetData;
	
///	class for a generator current (source) in a subdomain
	struct tGeneratorCurrent
	{
		SmartPtr<TGridFunction> m_spGf; ///< the grid function of the current
		size_t m_vFct[2]; ///< components of the grid function
		bool m_everywhere; ///< true iff the source is defined everywhere
		SubsetGroup m_ssGrp; ///< subsets where the source is defined (if ! m_everywhere)
		
	///	constructor
		tGeneratorCurrent
		(
			SmartPtr<TGridFunction> & spGf,
			size_t vFct_Re, size_t vFct_Im,
			SubsetGroup & ssGrp,
			bool ew = false
		)
		: m_spGf (spGf), m_everywhere (ew), m_ssGrp (ssGrp)
		{
			m_vFct[_Re_] = vFct_Re; m_vFct[_Im_] = vFct_Im;
		}
	};
	
///	array of all the sources (generator currents)
	std::vector<tGeneratorCurrent> m_vJG;

///	the source active in the current (assembled) subset
	tGeneratorCurrent * m_pSsJG;

//---- Temporary data used in the local discretization ----
private:
	
/// local stiffness matrix of the rot-rot operator
	number m_rot_rot_S [maxNumEdges][maxNumEdges];
/// local mass matrix of the rot-rot operator
	number m_rot_rot_M [maxNumEdges][maxNumEdges];

///	the magnetic permeability in the subdomain
	number m_permeability;

/// the electric conductivity in the subdomain
	number m_conductivity;

}; // end class EddyCurrent_E_Nedelec

///@}

} // end namespace Electromagnetism
} // end namespace ug

#include "eddy_current_e_nedelec_impl.h"

#endif /* __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_E_NEDELEC__ */

/* End of File */
