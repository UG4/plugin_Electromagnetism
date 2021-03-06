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
 * UserData class for computations with the time-harmonic Nedelec-based
 * representation of the electric field.
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

// Further local headers
#include "eddy_current_traits.h"
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
:	public StdDependentUserData<EddyCurrentHeat<TGridFunc>, number, TGridFunc::dim>,
	public EddyCurrentTraits
{
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
		const char * cmp, ///< the components (Re and Im) of the grid function keeping the scalar DoFs
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
		GridObject* elem,
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
		std::vector<DoFIndex> ind;
		std::vector<number> dofValues;
	
		for (size_t part = 0; part < 2; part++) // part: Re or Im
		{
			E[part].resize (nip);
			
		//	Get multiindices of element
			m_spGF->dof_indices (elem, m_fct [part], ind);
		
		//	The DoF values of the grid function
			dofValues.resize (ind.size ());
			for (size_t sh = 0; sh < dofValues.size (); ++sh)
				dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
			
		//	Compute the values of the grid function
			NedelecInterpolation<domain_type, refDim>::value
				(m_spGF->domain().get(), elem, vCornerCoords, &(dofValues[0]),
					vLocIP, nip, & (E[part][0]));
		}
		
	//	Compute \f$ \frac{1}{2} \sigma \mathbf{E} \overline{\mathbf{E}} \f$
		for (size_t i = 0; i < nip; i++)
			vValue [i] = sigma * (VecTwoNormSq (E[_Re_][i]) + VecTwoNormSq (E[_Im_][i])) / 2;
	};
};

/**
 * Template of UserData based classes for computation of vector fields that
 * depend on curls of components of E. The TImpl class should implement function
 * 'void get_value (const MathVector<dim> & curlE, MathVector<dim> & value) const'
 * that transforms \f$ \mathbf{rot} \, \mathbf {E} \f$ to the value.
 *
 * \tparam TImpl	implementation of the class (should be derived from EddyCurrentCurlEDependentCmpUserData)
 * \tparam ReIm		whether the UserData depends on \f$ \mathrm{Re} \, \mathbf{rot} \, \mathbf {E} \f$ or \f$ \mathrm{Im} \, \mathbf{rot} \, \mathbf {E} \f$
 * \tparam TGFunc	grid function type
 */
template <typename TImpl, size_t ReIm, typename TGFunc>
class EddyCurrentCurlEDependentCmpUserData
:	public StdDependentUserData <TImpl, MathVector<TGFunc::dim>, TGFunc::dim>,
	public EddyCurrentTraits
{
public:
///	Type of domain
	typedef typename TGFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

///	Type of position coordinates (e.g. position_type)
	typedef typename domain_type::position_type position_type;

protected:
///	grid function
	SmartPtr<TGFunc> m_spGF;
	
///	components (Re and Im) of the grid function
	size_t m_fct [2];

private:
///	const 'this' pointer of the implementation class
	const TImpl * this_impl () const {return static_cast<const TImpl*> (this);}

public:
/// constructor
	EddyCurrentCurlEDependentCmpUserData
	(
		SmartPtr<TGFunc> spGridFct, ///< grid function with the DoFs
		const char * cmp ///< the components (Re and Im) of the grid function keeping the scalar DoFs
	)
	: m_spGF (spGridFct)
	{
	//	get strings, tokenize them and select functions
		std::string fctString = std::string (cmp);
		std::vector<std::string> tokens;
		TokenizeString (fctString, tokens, ',');
		if ((int) tokens.size () != 2)
			UG_THROW("EddyCurrentReB: Needed 2 components "
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
		UG_CATCH_THROW ("EddyCurrent: Cannot find symbolic function "
						"component for the name '" << tokens[i] << "'.");
		
	//	check the function space of the grid function
		for (int i = 0; i < 2; i++)
			if (m_spGF->local_finite_element_id(m_fct[i]).type () != LFEID::NEDELEC)
				UG_THROW ("EddyCurrent: The function space of component "
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
		MathVector<dim> vValue[], ///< for the computed value
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject * elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector * u,
		bool bDeriv,
		int s,
		std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		const MathMatrix<refDim, dim> * vJT = NULL
	) const
	{
	//	Derivatives are not implemented
		if (bDeriv)
			UG_THROW ("EddyCurrentReB: Derivatives are not implemented.");
		
	//	Get multiindices of element
		std::vector<DoFIndex> ind;
		m_spGF->dof_indices (elem, m_fct [ReIm], ind);
	
	//	The DoF values of the grid function
		std::vector<number> dofValues (ind.size());
		for (size_t sh = 0; sh < dofValues.size (); ++sh)
			dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
		
	//	Compute the curl
		MathVector<dim> curl;
		NedelecInterpolation<domain_type, refDim>::curl
			(m_spGF->domain().get(), elem, vCornerCoords, &(dofValues[0]), curl);
		for (size_t ip = 0; ip < nip; ip++)
			this_impl()->get_value (curl, vValue [ip]);
	};
};

/**
 * UserData based class for the computation of the real part of the magnetic field
 * \f$ \mathbf{B} = \frac{i}{\omega} \mathbf{rot} \, \mathbf{E} =
 * - \frac{1}{\omega} \mathbf{rot} \, \mathrm{Re} \, \mathbf{E}
 * + i \frac{1}{\omega} \mathbf{rot} \, \mathrm{Im} \, \mathbf{E} \f$,
 * i.e. the vector field \f$ \mathrm{Re} \, \mathbf{B} = - \frac{1}{\omega}
 * \mathbf{rot} \, \mathrm{Im} \, \mathbf{E} \f$.
 */
template <typename TGridFunc>
class EddyCurrentReBofEUserData
	: public EddyCurrentCurlEDependentCmpUserData
		<EddyCurrentReBofEUserData<TGridFunc>, EddyCurrentTraits::_Im_, TGridFunc>
{
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

private:
///	Frequency \f$ \omega \f$
	number m_omega;

public:
/// constructor
	EddyCurrentReBofEUserData
	(
		SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
		const char * cmp, ///< the components (Re and Im) of the grid function keeping the scalar DoFs
		number omega ///< the frequency of the field
	)
	:	EddyCurrentCurlEDependentCmpUserData
		 <EddyCurrentReBofEUserData<TGridFunc>, EddyCurrentTraits::_Im_, TGridFunc> (spGridFct, cmp),
		m_omega (omega)
	{};

///	\f$ \mathrm{Re} \, \mathbf{B} = - \frac{1}{\omega} * \mathbf{rot} \, \mathrm{Im} \, \mathbf{E} \f$
	inline void get_value (const MathVector<dim> & ImCurlE, MathVector<dim> & ReB) const
	{
		ReB = ImCurlE;
		ReB /= - m_omega;
	}
};

/**
 * UserData based class for the computation of the imaginary part of the magnetic
 * field \f$ \mathbf{B} = \frac{i}{\omega} \mathbf{rot} \, \mathbf{E} =
 * - \frac{1}{\omega} \mathbf{rot} \, \mathrm{Re} \, \mathbf{E}
 * + i \frac{1}{\omega} \mathbf{rot} \, \mathrm{Im} \, \mathbf{E} \f$,
 * i.e. the vector field \f$ \mathrm{Im} \, \mathbf{B} = \frac{1}{\omega}
 * \mathbf{rot} \, \mathrm{Re} \, \mathbf{E} \f$.
 */
template <typename TGridFunc>
class EddyCurrentImBofEUserData
	: public EddyCurrentCurlEDependentCmpUserData
		<EddyCurrentImBofEUserData<TGridFunc>, EddyCurrentTraits::_Re_, TGridFunc>
{
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

private:
///	Frequency \f$ \omega \f$
	number m_omega;

public:
/// constructor
	EddyCurrentImBofEUserData
	(
		SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
		const char * cmp, ///< the components (Re and Im) of the grid function keeping the scalar DoFs
		number omega ///< the frequency of the field
	)
	:	EddyCurrentCurlEDependentCmpUserData
		 <EddyCurrentImBofEUserData<TGridFunc>, EddyCurrentTraits::_Re_, TGridFunc> (spGridFct, cmp),
		m_omega (omega)
	{};

///	\f$ \mathrm{Im} \, \mathbf{B} = \frac{1}{\omega} * \mathbf{rot} \, \mathrm{Re} \, \mathbf{E} \f$
	inline void get_value (const MathVector<dim> & ReCurlE, MathVector<dim> & ImB) const
	{
		ImB = ReCurlE;
		ImB /= m_omega;
	}
};

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_GF_USER_DATA__

/* End of File */
