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

/**
 * UserData class for convertion of the Nedelec-based data into the vector field.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_GF_USER_DATA__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_GF_USER_DATA__

#include <map>
#include <vector>

// ug4 headers
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/grid_function_user_data.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"

// Nedelec-type-1 headers:
#include "nedelec_local_ass.h"

#include "em_material.h"

namespace ug{
namespace Electromagnetism{

/**
 * UserData based class that computes the vector values of grid functions
 * keeping a Nedelec-element (Whitney-1) based representation of a vector field.
 */
template <typename TGridFunc>
class NedelecGridFunctionData
	: public StdDependentUserData
				<NedelecGridFunctionData<TGridFunc>, MathVector<TGridFunc::dim>, TGridFunc::dim>
{
public:
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

///	Type of position coordinates (e.g. position_type)
	typedef typename domain_type::position_type position_type;

private:
///	grid function
	SmartPtr<TGridFunc> m_spGF;

///	component of function
	size_t m_fct;

public:
/// constructor
	NedelecGridFunctionData
	(
		SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
		const char * cmp ///< the component of the grid function keeping the scalar DoFs
	)
	: m_spGF (spGridFct)
	{
	//	Get function id by name:
		try
		{
			m_fct = m_spGF->fct_id_by_name (cmp);
		}
		UG_CATCH_THROW ("NedelecGridFunctionData: Function space does not contain"
					" a function with name " << cmp << ".");

	//	Get local finite element id and check whether this is the Nedelec element:
		if (m_spGF->local_finite_element_id(m_fct).type () != LFEID::NEDELEC)
			UG_THROW ("NedelecGridFunctionData: Function " << cmp
					<< "is not based on the Nedelec element.");
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
		static const size_t maxEdges = NedelecInterpolation<domain_type, refDim>::maxNumEdges;
		
	//	Derivatives are not implemented
		if (bDeriv)
			UG_THROW ("NedelecGridFunctionData: Derivatives are not implemented.");
		
	//	Get multiindices of element
		std::vector<DoFIndex> ind;
		m_spGF->dof_indices (elem, m_fct, ind);
		if (ind.size() > maxEdges)
			UG_THROW ("NedelecGridFunctionData: Illegal number of DoFs!");
	
	//	The DoF values of the grid function
		number dofValues [maxEdges];
		for (size_t sh = 0; sh < ind.size(); ++sh)
			dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
		
	//	Compute the values
		NedelecInterpolation<domain_type, refDim>::value
			(m_spGF->domain().get(), elem, vCornerCoords, dofValues, vLocIP, nip, vValue);
	};
};

/**
 * UserData based class that computes the curl vector of grid functions
 * keeping a Nedelec-element (Whitney-1) based representation of a vector field.
 */
template <typename TGridFunc>
class NedelecCurlData
	: public StdDependentUserData
				<NedelecCurlData<TGridFunc>, MathVector<TGridFunc::dim>, TGridFunc::dim>
{
public:
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

///	Type of position coordinates (e.g. position_type)
	typedef typename domain_type::position_type position_type;

private:
///	grid function
	SmartPtr<TGridFunc> m_spGF;

///	component of function
	size_t m_fct;

public:
/// constructor
	NedelecCurlData
	(
		SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
		const char * cmp ///< the component of the grid function keeping the scalar DoFs
	)
	: m_spGF (spGridFct)
	{
	//	Get function id by name:
		try
		{
			m_fct = m_spGF->fct_id_by_name (cmp);
		}
		UG_CATCH_THROW ("NedelecCurlData: Function space does not contain"
					" a function with name " << cmp << ".");

	//	Get local finite element id and check whether this is the Nedelec element:
		if (m_spGF->local_finite_element_id(m_fct).type () != LFEID::NEDELEC)
			UG_THROW ("NedelecCurlData: Function " << cmp
					<< "is not based on the Nedelec element.");
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
		static const size_t maxEdges = NedelecInterpolation<domain_type, refDim>::maxNumEdges;
		
	//	Derivatives are not implemented
		if (bDeriv)
			UG_THROW ("NedelecCurlData: Derivatives are not implemented.");
		
	//	Get multiindices of element
		std::vector<DoFIndex> ind;
		m_spGF->dof_indices (elem, m_fct, ind);
		if (ind.size() > maxEdges)
			UG_THROW ("NedelecGridFunctionData: Illegal number of DoFs!");
	
	//	The DoF values of the grid function
		number dofValues [maxEdges];
		for (size_t sh = 0; sh < ind.size(); ++sh)
			dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
		
	//	Compute the curl
		MathVector<dim> curl;
		NedelecInterpolation<domain_type, refDim>::curl
			(m_spGF->domain().get(), elem, vCornerCoords, dofValues, curl);
		for (size_t ip = 0; ip < nip; ip++)
			vValue[ip] = curl;
	};
};

/**
 * UserData based class that computes vector values of grid functions multiplied
 * by the conductivity \f$ \sigma \f$. The grid functions should keep a
 * Nedelec-element (Whitney-1) based representation of a vector field. If
 * the grid functions represent the electric field E, then this UserData
 * computes the current (not taking into account the generator current).
 */
template <typename TGridFunc>
class NedelecSigmaEData
	: public StdDependentUserData
				<NedelecSigmaEData<TGridFunc>, MathVector<TGridFunc::dim>, TGridFunc::dim>
{
public:
///	Type of domain
	typedef typename TGridFunc::domain_type domain_type;

///	World dimension
	static const int dim = domain_type::dim;

///	Type of position coordinates (e.g. position_type)
	typedef typename domain_type::position_type position_type;

private:
///	grid function
	SmartPtr<TGridFunc> m_spGF;

///	component of function
	size_t m_fct;
	
///	subdomain-dependent data (propertiels of the materials)
	SmartPtr<EMaterial<domain_type> > m_spEMaterial;

public:
/// constructor
	NedelecSigmaEData
	(
		SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
		const char * cmp, ///< the component of the grid function keeping the scalar DoFs
		SmartPtr<EMaterial<domain_type> > emMatherial ///< properties of the materials
	)
	: m_spGF (spGridFct), m_spEMaterial (emMatherial)
	{
	//	Get function id by name:
		try
		{
			m_fct = m_spGF->fct_id_by_name (cmp);
		}
		UG_CATCH_THROW ("NedelecGridFunctionData: Function space does not contain"
						" a function with name " << cmp << ".");

	//	Get local finite element id and check whether this is the Nedelec element:
		if (m_spGF->local_finite_element_id(m_fct).type () != LFEID::NEDELEC)
			UG_THROW ("NedelecGridFunctionData: Function " << cmp
					<< "is not based on the Nedelec element.");
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
		static const size_t maxEdges = NedelecInterpolation<domain_type, refDim>::maxNumEdges;
		
	//	Derivatives are not implemented
		if (bDeriv)
			UG_THROW ("NedelecSigmaEData: Derivatives are not implemented.");
		
	//	Get the material data
		number mu, sigma;
		if (m_spEMaterial->get_mu_sigma (si, mu, sigma))
			UG_THROW ("NedelecSigmaEData: No material data set to subset" << si
				<< " (or this is a low-dim. domain).");
		if (sigma == 0) // no currents in insulators
		{
			for (size_t i = 0; i < nip; i++) vValue [i] = 0.0;
			return;
		}
		
	//	Get multiindices of element
		std::vector<DoFIndex> ind;
		m_spGF->dof_indices (elem, m_fct, ind);
		if (ind.size() > maxEdges)
			UG_THROW ("NedelecGridFunctionData: Illegal number of DoFs!");
	
	//	The DoF values of the grid function
		number dofValues [maxEdges];
		for (size_t sh = 0; sh < ind.size(); ++sh)
			dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
		
	//	Compute the values of the grid function
		NedelecInterpolation<domain_type, refDim>::value
			(m_spGF->domain().get(), elem, vCornerCoords, dofValues, vLocIP, nip, vValue);
		
	//	Multiply the values by \f$ \sigma \f$
		for (size_t i = 0; i < nip; i++)
			vValue [i] *= sigma;
	};
};

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_GF_USER_DATA__

/* End of File */
