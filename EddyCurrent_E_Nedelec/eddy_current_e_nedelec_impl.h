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
	if (vLfeID[0] != LFEID (LFEID::NEDELEC, 1)
		|| vLfeID[1] != LFEID (LFEID::NEDELEC, 1))
		UG_THROW ("EddyCurrent_E_Nedelec: This discretization works with the Nedelec-elements only.");
}

////////////////////////////////////////////////////////////////////////////////
// assembling
////////////////////////////////////////////////////////////////////////////////

// prepares the loop over the elements: checks whether the imports are set, ...
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::prepare_element_loop
(
	ReferenceObjectID roid, // only elements with this roid are looped over
	int si // and only in this subdomain
)
{
//	get the data for the subdomain
	if (m_spSubsetData->get_mu_sigma (si, m_permeability, m_conductivity))
		UG_THROW ("Cannot get material data for subset " << si << ".");
}

// finishes the loop over the elements: nothing to do
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::finish_element_loop ()
{}

// prepares a given element for assembling: computes the discretization of the rot-rot operator
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::prepare_element
(
	const LocalVector & u, // local solution at the dofs associated with elem
	GeometricObject * elem, // element to prepare
	const position_type vCornerCoords [] // coordinates of the corners of the element
)
{
/* BEGIN code essential for the numerics */

//	compute the local matrices of the rot-rot operator:
	NedelecT1_LDisc<TDomain, TElem>::local_stiffness_and_mass
		(* EddyCurrent_E_Nedelec<TDomain,TAlgebra>::domain().get(),
			static_cast<TElem*> (elem), vCornerCoords, m_rot_rot_S, m_rot_rot_M);
			
/* END code essential for the numerics */
}

// compose the local stiffness matrix
/**
 * This function composes the local stiffness matrix of the stationary problem
 * by combining the stiffness and the mass matrix of the discretization
 * of the rot-rot operator. Followin O. Sterz, the system has the following
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
 *   - \omega \mathbf{Im} \mathbf{J}_{Gh} \\
 *   - \omega \mathbf{Re} \mathbf{J}_{Gh}
 *  \end{array}
 * \right ),
 * \f}
 * where \f$S_h (\mu^{-1})\f$ is the stiffness matrix of the \f$\mathbf{rot}\mu^{-1}\mathbf{rot}\f$,
 * \f$ M^{(1)}_h (\sigma)\f$ the mass matrix of the \f$\sigma I\f$ operator,
 * \f$\mathbf{E}_h\f$ the unknown electric field grid function,
 * \f$\mathbf{J}_{Gh}\f$ the generator current.
 */
template<typename TDomain, typename TAlgebra>
template<size_t numEdges>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_elem_stiffness
(
	number perm, // the magnetic permeability
	number cond, // the electric conductivity
	number S [2][numEdges] [2][numEdges] // for the composed matrix
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

// transfer the local stiffness matrix to the global discretization
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

// compute the right-hand side due to the generator currents
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_rhs_elem
(
	LocalVector& b,
	GeometricObject * elem,
	const position_type vCornerCoords []
)
{
	if (m_spgfJG.invalid ()) return; // no rhs
	
//	the generator current at the dofs
	MathVector<2> vJG [NedelecT1_LDisc<TDomain, TElem>::numEdges];
	
//	get the generator currend from the grid function
	TElem * pElem = static_cast<TElem*> (elem);
	for (size_t i = 0; i < 2; i++) // Re and Im
	{
		std::vector<MultiIndex<2> > ind;
		m_spgfJG->multi_indices (pElem, m_vfctJG[i], ind);
		for (size_t e = 0; e < NedelecT1_LDisc<TDomain, TElem>::numEdges; e++)
			vJG[e][i] = DoFRef (*m_spgfJG, ind[e]);
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

// compute the local defect and transfer it to the global discretization
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_dA_elem
(
	LocalVector & d, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords []
)
{}

// This is a stationary problem, so there is no mass matrix in this sence.
// For this, the following two functions do nothing.
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::ass_JM_elem
(
	LocalMatrix & J, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords []
)
{}
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

// register the local assembler functions for all the elements and dimensions
template<typename TDomain, typename TAlgebra>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::register_all_loc_discr_funcs ()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<dim>::AllElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList> (RegisterLocalDiscr (this));
}

// register the local assembler functions for a given element
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

// the generator current \f$ \mathbf{J}_{G,h} \f$
template<typename TDomain, typename TAlgebra>
void EddyCurrent_E_Nedelec<TDomain,TAlgebra>::set_generator_current
(
	SmartPtr<TGridFunction> spgfJG, ///< pointer to the grid function
	const char* cmp ///< names of the components
)
{
//	get strings
	std::string fctString = std::string(cmp);

//	tokenize strings and select functions
	std::vector<std::string> tokens;
	TokenizeString(fctString, tokens, ',');

	if((int)tokens.size() != 2)
		UG_THROW("EddyCurrent_E_Nedelec: Needed 2 components "
				 "in symbolic function names for the generator current (for Re and Im), "
				 "but given: " << cmp);

	for (size_t i = 0; i < tokens.size(); i++)
		RemoveWhitespaceFromString(tokens[i]);

//	get function id of name
	for (int i = 0; i < 2; i++)
		try
		{
			m_vfctJG[i] = spgfJG->fct_id_by_name(tokens[i].c_str());
		}
		UG_CATCH_THROW ("EddyCurrent_E_Nedelec: Cannot find symbolic function "
						"component for the name '" << tokens[i] << "'.");
	
//	check the function space of the grid function
	for (int i = 0; i < 2; i++)
		if (spgfJG->local_finite_element_id(m_vfctJG[i]) != LFEID(LFEID::NEDELEC, 1))
			UG_THROW ("EddyCurrent_E_Nedelec: The function space of component "
							<< tokens[i] << " of the grid function does not correspond "
							"to the Nedelec element.");
	
//	set the inner pointer to the grid function
	m_spgfJG = spgfJG;
}

////////////////////////////////////////////////////////////////////////////////
//	constructor
////////////////////////////////////////////////////////////////////////////////
template<typename TDomain, typename TAlgebra>
EddyCurrent_E_Nedelec<TDomain,TAlgebra>::
EddyCurrent_E_Nedelec
(
	const char* functions, ///< grid function names
	ConstSmartPtr<EMaterial<domain_type> > spSubsetData, ///< material data object
	number frequency ///< frequency \f$\omega\f$
)
:	IElemDisc<TDomain> (functions, spSubsetData->subset_names ()),
	m_omega (frequency),
	m_spSubsetData (spSubsetData)
{
//	check number of functions
	if(this->num_fct() != 2)
		UG_THROW("Wrong number of functions: The ElemDisc 'EddyCurrent_E_Nedelec'"
					" needs exactly 2 symbolic function"
					" (one for the real part and one for the imaginary one).");

//	register assemble functions
	register_all_loc_discr_funcs ();
}
	
} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
