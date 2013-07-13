/**
 * UserData class for computations with the time-harmonic Nedelec-based
 * representation of the electric field.
 *
 * Created on: 23.06.2013
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_GF_USER_DATA__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_GF_USER_DATA__

#include <map>
#include <vector>

// ug4 headers
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"

// Nedelec-type-1 headers:
#include "../nedelec_local_ass.h"

#include "../em_material.h"

namespace ug{
namespace Electromagnetism{

/**
 * UserData based class that computes the field of the heat sources
 * \f$Q = \frac{1}{2} \sigma \mathbf{E} \overline{\mathbf{E}}\f$, where
 * \f$\overline{\mathbf{E}}\f$ is the complex conjugate to \f$\mathbf{E}\f$.
 * The electric field \f$\mathbf{E}\f$ in the grid function should be
 * represented by the Nedelec-type-1 elements.
 */
template <typename TGridFunc>
class EddyCurrentHeat
	: public StdDependentUserData<EddyCurrentHeat<TGridFunc>, number, TGridFunc::dim>
{
/// index of the real part in the grid functions
	static const size_t _Re_ = 0;
/// index of the imaginary part in the grid functions
	static const size_t _Im_ = 1;

public:
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

///	Type of position coordinates (e.g. position_type)
	typedef typename domain_type::position_type position_type;

private:
///	grid function for \f$ \mathbf{E} \f$
	SmartPtr<TGridFunc> m_spGF;

///	components (Re and Im) of the grid function
	size_t m_fct [2];
	
///	subdomain-dependent data (propertiels of the materials)
	SmartPtr<EMaterial<domain_type> > m_spEMaterial;

public:
/// constructor
	EddyCurrentHeat
	(
		SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
		const char* cmp, ///< the components (Re and Im) of the grid function keeping the scalar DoFs
		SmartPtr<EMaterial<domain_type> > emMatherial ///< properties of the materials
	)
	: m_spGF (spGridFct), m_spEMaterial (emMatherial)
	{
	//	get strings, tokenize them and select functions
		std::string fctString = std::string (cmp);
		std::vector<std::string> tokens;
		TokenizeString (fctString, tokens, ',');
		if ((int) tokens.size () != 2)
			UG_THROW("EddyCurrentHeat: Needed 2 components "
						"in symbolic function names (for Re and Im), "
							"but given: " << cmp);
		for (size_t i = 0; i < tokens.size (); i++)
			RemoveWhitespaceFromString (tokens [i]);
	
	//	get function id's by names
		for (int i = 0; i < 2; i++)
		try
		{
			m_fct [i] = m_spGF->fct_id_by_name (tokens[i].c_str ());
		}
		UG_CATCH_THROW ("EddyCurrentHeat: Cannot find symbolic function "
						"component for the name '" << tokens[i] << "'.");
		
	//	check the function space of the grid function
		for (int i = 0; i < 2; i++)
			if (m_spGF->local_finite_element_id(m_fct[i]).type () != LFEID::NEDELEC)
				UG_THROW ("EddyCurrentHeat: The function space of component "
							<< tokens[i] << " of the grid function does not correspond "
								"to the Nedelec element.");
	};
	
///	The vector field retrieved from the Nedelec-type 1 (Whitney-1) dofs are not continuous
	virtual bool continuous () const {return false;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return true;}

/// Performs the main computations:
	template <int refDim>
	void eval_and_deriv
	(
		number vValue[], ///< for the computed value
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GeometricObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		bool bDeriv,
		int s,
		std::vector<std::vector<number> > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
	//	Derivatives are not implemented
		if (bDeriv)
			UG_THROW ("EddyCurrentHeat: Derivatives are not implemented.");
		
	//	Get the material data
		number mu, sigma;
		int elem_si = m_spGF->approx_space()->subset_handler()->get_subset_index (elem);
		if (m_spEMaterial->get_mu_sigma (elem_si, mu, sigma))
			UG_THROW ("EddyCurrentHeat: No material data set to subset" << si
				<< " (or this is a low-dim. domain).");
		if (sigma == 0) // no currents in insulators
		{
			for (size_t i = 0; i < nip; i++) vValue [i] = 0.0;
			return;
		}
		
	//	Compute the values
		std::vector<MathVector<dim> > E [2]; // Re and Im parts of E
		std::vector<MultiIndex<2> > ind;
		std::vector<number> dofValues;
	
		for (size_t part = 0; part < 2; part++) // part: Re or Im
		{
			E[part].resize (nip);
			
		//	Get multiindices of element
			m_spGF->multi_indices (elem, m_fct [part], ind);
		
		//	The DoF values of the grid function
			dofValues.resize (ind.size ());
			for (size_t sh = 0; sh < dofValues.size (); ++sh)
				dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
			
		//	Compute the values of the grid function
			NedelecInterpolation<domain_type, refDim>::value
				(*(m_spGF->domain().get()), elem, vCornerCoords, &(dofValues[0]),
					vLocIP, nip, & (E[part][0]));
		}
		
	//	Compute \f$ \frac{1}{2} \sigma \mathbf{E} \overline{\mathbf{E}} \f$
		for (size_t i = 0; i < nip; i++)
			vValue [i] = sigma * (VecTwoNormSq (E[_Re_][i]) + VecTwoNormSq (E[_Im_][i])) / 2;
	};
};

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_GF_USER_DATA__

/* End of File */
