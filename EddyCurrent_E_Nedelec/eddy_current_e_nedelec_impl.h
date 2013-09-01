/*
 * Implementation of the class member functions templates of the complex-valued
 * FE discretization of the eddy current model.
 * Created on: 18.09.2012
 * Author: D. Logashenko
 */

/* UG4 headers: */
#include "common/util/provider.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace Electromagnetism{

////////////////////////////////////////////////////////////////////////////////
//	check the grid and the shape functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain, typename TAlgebra>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::prepare_setting
(
	const std::vector<LFEID> & vLfeID,
	bool bNonRegular
)
{
//	check the grid
	if (bNonRegular)
		UG_THROW ("ERROR in 'EddyCurrent_E_Nedelec::request_non_regular_grid':"
				" The discretization does not support hanging nodes.\n");

//	check number of the components
	if (vLfeID.size () != 2)
		UG_THROW ("EddyCurrent_E_Nedelec: The time-harmonic problem needs two components at every edge.");
		/* one component for the real part and one for the imaginary part */

//	check that these are the Nedelec elements
	if (vLfeID[0].type() != LFEID::NEDELEC || vLfeID[1].type() != LFEID::NEDELEC)
		UG_THROW ("EddyCurrent_E_Nedelec: This discretization works with the Nedelec-elements only.");
}

////////////////////////////////////////////////////////////////////////////////
// assembling
////////////////////////////////////////////////////////////////////////////////

/// prepares the loop over the elements: checks whether the parameters are set, ...
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::prepare_element_loop
(
	ReferenceObjectID roid, ///< only elements with this roid are looped over
	int si ///< and only in this subdomain
)
{
//	Get the data for the subset:
	if (m_spSubsetData->get_mu_sigma (si, m_permeability, m_conductivity))
		UG_THROW ("Cannot get material data for subset " << si << ".");
	
//	Get the source for this subset
	m_pSsJG = NULL;
	for (size_t i = 0; i < m_vJG.size (); i++)
	if (m_vJG[i].m_everywhere || m_vJG[i].m_ssGrp.contains (si))
	{
		if (m_pSsJG != NULL)
			UG_THROW ("More than one specification of the generator current in subset " << si << ".");
		m_pSsJG = & (m_vJG[i]);
	}
}

/// finalizes the loop over the elements: clear the source
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::finish_element_loop ()
{
	m_pSsJG = NULL;
}

/// prepares a given element for assembling: computes the discretization of the rot-rot operator
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::prepare_element
(
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GeometricObject * elem, ///< element to prepare
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
/* BEGIN code essential for the numerics */

//	compute the local matrices of the rot-rot operator:
	NedelecT1_LDisc<TDomain, TElem>::local_stiffness_and_mass
		(* EddyCurrent_E_Nedelec<TDomain,TAlgebra>::domain().get(),
			static_cast<TElem*> (elem), vCornerCoords, m_rot_rot_S, m_rot_rot_M);
			
/* END code essential for the numerics */
}

/// composes the local stiffness matrix
/**
 * This function composes the local stiffness matrix of the stationary problem
 * by combining the stiffness and the mass matrix of the discretization
 * of the rot-rot operator. Following O. Sterz, the system has the following
 * block form:
 * \f{eqnarray*}{
 * \left (
 *  \begin{array} {cc}
 *   - S_h (\mu^{-1}) & \omega M^{(1)}_h (\sigma) \\
 *   \omega M^{(1)}_h (\sigma) & S_h (\mu^{-1})
 *  \end{array}
 * \right )
 * \left (
 *  \begin{array} {c}
 *   \mathbf{Re} \mathbf{E}_h \\
 *   \mathbf{Im} \mathbf{E}_h
 *  \end{array}
 * \right )
 * =
 * \left (
 *  \begin{array} {c}
 *   - \omega M^{(1)}_h (1) \mathbf{Im} \mathbf{J}_{Gh} \\
 *   - \omega M^{(1)}_h (1) \mathbf{Re} \mathbf{J}_{Gh}
 *  \end{array}
 * \right ),
 * \f}
 * where
 * <ul>
 * <li> \f$S_h (\mu^{-1})\f$		the stiffness matrix of the \f$\mathbf{rot}\mu^{-1}\mathbf{rot}\f$,
 * <li> \f$ M^{(1)}_h (\sigma)\f$	the mass matrix of the \f$\sigma I\f$ operator,
 * <li> \f$\mathbf{E}_h\f$			the unknown electric field grid function (represented using the Nedelec element),
 * <li> \f$\mathbf{J}_{Gh}\f$		the generator current (represented using the Nedelec element).
 * </ul>
 *
 * \remark
 * For the solvability of the problem, the generator current must be numerically
 * weakly divergence free.
 */
template<typename TDomain, typename TAlgebra>
template<size_t numEdges>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_elem_stiffness
(
	number perm, ///< the magnetic permeability
	number cond, ///< the electric conductivity
	number S [2][numEdges] [2][numEdges] ///< for the composed matrix
)
{
/* BEGIN code essential for the numerics */
	
	cond *= m_omega;
	
// the lower triange (w.r.t. the edge dofs, not the Re/Im parts)
	for (size_t e_1 = 0; e_1 < numEdges; e_1++)
		for (size_t e_2 = 0; e_2 <= e_1; e_2++)
		{
			S [_Re_][e_1] [_Re_][e_2]
				= - (S [_Im_][e_1] [_Im_][e_2]
					= m_rot_rot_S [e_1] [e_2] / perm);
			S [_Re_][e_1] [_Im_][e_2]
				= S [_Im_][e_1] [_Re_][e_2]
					= m_rot_rot_M [e_1] [e_2] * cond;
		}
// the upper triangle
	for (size_t e_1 = 0; e_1 < numEdges - 1; e_1++)
		for (size_t e_2 = e_1+1; e_2 < numEdges; e_2++)
		{
			S [_Re_][e_1] [_Re_][e_2] = S [_Re_][e_2] [_Re_][e_1];
			S [_Re_][e_1] [_Im_][e_2] = S [_Re_][e_2] [_Im_][e_1];
			S [_Im_][e_1] [_Im_][e_2] = S [_Im_][e_2] [_Im_][e_1];
			S [_Im_][e_1] [_Re_][e_2] = S [_Im_][e_2] [_Re_][e_1];
		}
	
/* END code essential for the numerics */
}

/// transfers the precomputed local stiffness matrix to the global discretization
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_JA_elem
(
	LocalMatrix & J,
	const LocalVector & u,
	GeometricObject * elem,
	const position_type vCornerCoords []
)
{
	number S [2][NedelecT1_LDisc<TDomain, TElem>::numEdges]
	         [2][NedelecT1_LDisc<TDomain, TElem>::numEdges];
	
	ass_elem_stiffness<NedelecT1_LDisc<TDomain, TElem>::numEdges>
		(m_permeability, m_conductivity, S);

	for (size_t e_1 = 0; e_1 < NedelecT1_LDisc<TDomain, TElem>::numEdges; e_1++)
		for (size_t e_2 = 0; e_2 < NedelecT1_LDisc<TDomain, TElem>::numEdges; e_2++)
		{
			J(_Re_, e_1, _Re_, e_2) += S [_Re_][e_1] [_Re_][e_2];
			J(_Re_, e_1, _Im_, e_2) += S [_Re_][e_1] [_Im_][e_2];
			J(_Im_, e_1, _Im_, e_2) += S [_Im_][e_1] [_Im_][e_2];
			J(_Im_, e_1, _Re_, e_2) += S [_Im_][e_1] [_Re_][e_2];
		}
}

/// computes the right-hand side due to the generator currents
/**
 * The right-hand side has the form
 * \f{eqnarray*}{
 * \left (
 *  \begin{array} {c}
 *   - \omega M^{(1)}_h \mathbf{Im} \mathbf{J}_{Gh} \\
 *   - \omega M^{(1)}_h \mathbf{Re} \mathbf{J}_{Gh}
 *  \end{array}
 * \right ),
 * \f}
 * where
 * <ul>
 * <li> \f$ M^{(1)}_h\f$		the mass matrix of the \f$I\f$ (identity) operator,
 * <li> \f$\mathbf{J}_{Gh}\f$	the generator current (represented using the Nedelec element).
 * </ul>
 *
 * \remark
 * For the solvability of the problem, the generator current must be numerically
 * weakly divergence free.
 */
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_rhs_elem
(
	LocalVector& b,
	GeometricObject * elem,
	const position_type vCornerCoords []
)
{
	if (m_pSsJG == NULL) return; // no rhs
	
//	The generator current (source) at the dofs
	MathVector<2> vJG [NedelecT1_LDisc<TDomain, TElem>::numEdges];
	
//	Get the generator current from the grid functions
	TElem * pElem = static_cast<TElem*> (elem);
	for (size_t part = 0; part < 2; part++) // Re and Im
	{
		std::vector<MultiIndex<2> > ind;
		m_pSsJG->m_spGf->multi_indices (pElem, m_pSsJG->m_vFct[part], ind);
		for (size_t e = 0; e < NedelecT1_LDisc<TDomain, TElem>::numEdges; e++)
			vJG[e][part] = DoFRef (* (m_pSsJG->m_spGf), ind[e]);
	}
	
/* BEGIN code essential for the numerics */

//	Assemble the possibly non-zero rhs:
	for (size_t e_1 = 0; e_1 < NedelecT1_LDisc<TDomain, TElem>::numEdges; e_1++)
	{
		number Re_rhs = 0, Im_rhs = 0;
		for (size_t e_2 = 0; e_2 < NedelecT1_LDisc<TDomain, TElem>::numEdges; e_2++)
		{
			Re_rhs -= m_rot_rot_M[e_1][e_2] * vJG[e_2][_Im_];
			Im_rhs -= m_rot_rot_M[e_1][e_2] * vJG[e_2][_Re_];
		}
		b (_Re_, e_1) += Re_rhs; b (_Im_, e_1) += Im_rhs;
	}

/* END code essential for the numerics */
}

/// computes the local defect and transfer it to the global discretization
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_dA_elem
(
	LocalVector & d, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords []
)
{}

/// computes the mass matrix of a time-dependent problem
/**
 * This is a stationary problem, so there is no mass matrix in this sence.
 * For this, the following function does nothing.
 */
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_JM_elem
(
	LocalMatrix & J, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords []
)
{}

/// computes the mass part of the defect of a time-dependent problem
/**
 * This is a stationary problem, so there is no mass matrix in this sence.
 * For this, the following function does nothing.
 */
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_dM_elem
(
	LocalVector & d, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords []
)
{}

////////////////////////////////////////////////////////////////////////////////
//	register assembling functions
////////////////////////////////////////////////////////////////////////////////

/// registers the local assembler functions for all the elements and dimensions
template<typename TDomain, typename TAlgebra>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::register_all_loc_discr_funcs ()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::AllElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList> (RegisterLocalDiscr (this));
}

/// registers the local assembler functions for a given element
template<typename TDomain, typename TAlgebra>
template<typename TElem> // the element to register for
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::register_loc_discr_func ()
{
	static const ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	this->enable_fast_add_elem (true);
	this->set_prep_elem_loop_fct(id, & this_type::template prepare_element_loop<TElem>);
	this->set_prep_elem_fct		(id, & this_type::template prepare_element<TElem>);
	this->set_fsh_elem_loop_fct	(id, & this_type::template finish_element_loop<TElem>);
	this->set_add_jac_A_elem_fct(id, & this_type::template ass_JA_elem<TElem>);
	this->set_add_jac_M_elem_fct(id, & this_type::template ass_JM_elem<TElem>);
	this->set_add_def_A_elem_fct(id, & this_type::template ass_dA_elem<TElem>);
	this->set_add_def_M_elem_fct(id, & this_type::template ass_dM_elem<TElem>);
	this->set_add_rhs_elem_fct	(id, & this_type::template ass_rhs_elem<TElem>);
}

////////////////////////////////////////////////////////////////////////////////
//	obtaining the parameters
////////////////////////////////////////////////////////////////////////////////

/// adds a generator current item \f$ \mathbf{J}_{G,h} \f$ to the discretization
template<typename TDomain, typename TAlgebra>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::set_generator_current
(
	SmartPtr<TGridFunction> spgfJG, ///< pointer to the grid function
	const char * cmp, ///< names of the components
	const char * ss_names ///< names of the subsets (or NULL if defined everywhere)
)
{
	std::vector<std::string> tokens;
	
	if (spgfJG.invalid ())
		UG_THROW("EddyCurrent_E_Nedelec: Invalid grid function for the generator source");
	
//	Data to get:
	size_t vfctJG [2];
	bool ew;
	SubsetGroup ssGrp (spgfJG->approx_space()->subset_handler ());
	
//	Get function names:
	std::string fctString (cmp);

//	Tokenize the strings and select functions
	TokenizeString (fctString, tokens, ',');

	if ((int) tokens.size () != 2)
		UG_THROW("EddyCurrent_E_Nedelec: Needed 2 components "
				 "in symbolic function names for the generator current (for Re and Im), "
				 "but given: " << cmp);

	for (size_t i = 0; i < 2; i++)
		RemoveWhitespaceFromString (tokens [i]);

//	Get function id's by names:
	for (int i = 0; i < 2; i++)
	try
	{
		vfctJG [i] = spgfJG->fct_id_by_name (tokens[i].c_str ());
	}
	UG_CATCH_THROW ("EddyCurrent_E_Nedelec: Cannot find symbolic function "
					"component for the name '" << tokens[i] << "'.");
	
//	Check the function space of the grid function:
	for (int i = 0; i < 2; i++)
		if (spgfJG->local_finite_element_id(vfctJG[i]) != LFEID(LFEID::NEDELEC, dim, 1))
			UG_THROW ("EddyCurrent_E_Nedelec: The function space of component "
							<< tokens[i] << " of the grid function does not correspond "
							"to the Nedelec element.");
	
//	Get the subsets:
	if (ss_names == NULL)
		ew = true;
	else
	{
		ew = false;
		std::string ssString (ss_names);
		TokenizeString (ssString, tokens, ',');
		for (size_t i = 0; i < tokens.size (); i++)
			RemoveWhitespaceFromString (tokens [i]);
		try
		{
			ssGrp.add (tokens);
		}
		UG_CATCH_THROW
			("EddyCurrent_E_Nedelec: Failed to add subsets to the source");
	}
		
//	Create the source:
	m_vJG.push_back (tGeneratorCurrent (spgfJG, vfctJG[_Re_], vfctJG[_Im_], ssGrp, ew));
}

////////////////////////////////////////////////////////////////////////////////
///	class constructor
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain, typename TAlgebra>
EddyCurrent_E_Nedelec<TDomain,TAlgebra>::
EddyCurrent_E_Nedelec
(
	const char * functions, ///< grid function names
	ConstSmartPtr<EMaterial<domain_type> > spSubsetData, ///< material data object
	number frequency ///< frequency \f$\omega\f$
)
:	IElemDisc<TDomain> (functions, spSubsetData->subset_names ()),
	m_omega (frequency),
	m_spSubsetData (spSubsetData)
{
//	check number of functions
	if (this->num_fct () != 2)
		UG_THROW ("Wrong number of functions: The ElemDisc 'EddyCurrent_E_Nedelec'"
					" needs exactly 2 symbolic function"
					" (one for the real part and one for the imaginary one).");

//	register assemble functions
	register_all_loc_discr_funcs ();
}
	
} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
