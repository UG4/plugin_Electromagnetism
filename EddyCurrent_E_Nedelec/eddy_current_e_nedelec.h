/**
 * FE-discretization of the time-harmonic E-based formulation of the eddy
 * current model.
 *
 * Created on: 17.09.2012
 * Author: D. Logashenko
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
 * Note that every Nedelec shape function requires a scalar (complex) dof value.
 * These functions are thus vector-valued for the scalar dofs, each dof storing
 * two doubles (one for the real and one for the imaginary part of the value).
 *
 * References:
 * <ul>
 * <li> O. Sterz. Modellierung und Numerik zeitharmonischer Wirbelstromprobleme in 3D. PhD thesis, 2003.
 * <li> A. Bossavit. Computational Electromagnetism. Academic Press (Boston), 1998 (available in the Internet)
 * </ul>
 *
 * \tparam	TDomain		Domain type
 * \tparam	TAlgebra	Algebra type
 */
template <typename TDomain, typename TAlgebra>
class EddyCurrent_E_Nedelec
	: public IElemDisc<TDomain>
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

///	world dimension
	static const int dim = base_type::dim;

///	position type
	typedef typename base_type::position_type position_type;

/// indices of the real and the imaginary parts in the grid functions
	static const size_t _Re_ = 0;
	static const size_t _Im_ = 1;

public:
///	constructor
	EddyCurrent_E_Nedelec
	(
		const char* functions,
		ConstSmartPtr<EMaterial<domain_type> > spSubsetData,
		number frequency
	);

private:
/// frequency \f$\omega\f$ for the discretization
	number m_omega;

/// parameters of the materials in the domain
	ConstSmartPtr<EMaterial<domain_type> > m_spSubsetData;
	
/// the generator current \f$ \mathbf{J}_{G,h} \f$
	SmartPtr<TGridFunction> m_spgfJG; ///< the grid function
	size_t m_vfctJG[2]; ///< components of the grid function

public:
/// sets the generator current \f$ \mathbf{J}_{G,h} \f$
	void set_generator_current(SmartPtr<TGridFunction> spgfJG, const char* cmp);

//---- Local discretization interface: ----
private:
///	check type of the grid and the trial space
	virtual void prepare_setting
	(
		const std::vector<LFEID> & vLfeID,
		bool bNonRegular
	);

/// assembling functions
/// \{
	template <typename TElem>
	void prepare_element_loop(ReferenceObjectID roid, int si);

	template <typename TElem>
	void prepare_element(const LocalVector& u, GeometricObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void finish_element_loop();

	template <typename TElem>
	void ass_JA_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_JM_elem(LocalMatrix& J, const LocalVector& u, GeometricObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dA_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_dM_elem(LocalVector& d, const LocalVector& u, GeometricObject* elem, const position_type vCornerCoords[]);

	template <typename TElem>
	void ass_rhs_elem(LocalVector& d, GeometricObject* elem, const position_type vCornerCoords[]);
/// \}

	private:
//---- Registration of the template functions: ----
	void register_all_loc_discr_funcs();

	struct RegisterLocalDiscr {
			RegisterLocalDiscr(this_type* pThis) : m_pThis(pThis){}
			this_type* m_pThis;
			template< typename TElem > void operator()(TElem&)
			{m_pThis->register_loc_discr_func<TElem>();}
	};

	template <typename TElem>
	void register_loc_discr_func();

	private:
//---- Auxiliary functions: ----

/// composes the stiffness matrix of the stationary problem
	template<size_t numEdges>
	void ass_elem_stiffness
	(
		number perm, ///< the magnetic permeability
		number cond, ///< the electric conductivity
		number S [2][numEdges] [2][numEdges] ///< for the composed matrix
	);

//---- Temporary data used in the local discretization ----
	
/// local stiffness matrix of the rot-rot operator
	number m_rot_rot_S [ROT_ROT_MAX_EDGES][ROT_ROT_MAX_EDGES];
/// local mass matrix of the rot-rot operator
	number m_rot_rot_M [ROT_ROT_MAX_EDGES][ROT_ROT_MAX_EDGES];

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
