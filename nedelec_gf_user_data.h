/**
 * UserData class for convertion of the Nedelec-based data into the vector field.
 *
 * Created on: 20.02.2013
 * Author: D. Logashenko
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
		GeometricObject * elem,
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
			UG_THROW ("NedelecGridFunctionData: Derivatives are not implemented.");
		
	//	Get multiindices of element
		std::vector<DoFIndex> ind;
		m_spGF->dof_indices (elem, m_fct, ind);
	
	//	The DoF values of the grid function
		std::vector<number> dofValues (ind.size());
		for (size_t sh = 0; sh < dofValues.size (); ++sh)
			dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
		
	//	Compute the values
		NedelecInterpolation<domain_type, refDim>::value
			(*(m_spGF->domain().get()), elem, vCornerCoords, &(dofValues[0]), vLocIP, nip, vValue);
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
		GeometricObject * elem,
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
			UG_THROW ("NedelecCurlData: Derivatives are not implemented.");
		
	//	Get multiindices of element
		std::vector<DoFIndex> ind;
		m_spGF->dof_indices (elem, m_fct, ind);
	
	//	The DoF values of the grid function
		std::vector<number> dofValues (ind.size());
		for (size_t sh = 0; sh < dofValues.size (); ++sh)
			dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
		
	//	Compute the curl
		MathVector<dim> curl;
		NedelecInterpolation<domain_type, refDim>::curl
			(*(m_spGF->domain().get()), elem, vCornerCoords, &(dofValues[0]), curl);
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
		GeometricObject * elem,
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
	
	//	The DoF values of the grid function
		std::vector<number> dofValues (ind.size());
		for (size_t sh = 0; sh < dofValues.size (); ++sh)
			dofValues[sh] = DoFRef (*m_spGF, ind[sh]);
		
	//	Compute the values of the grid function
		NedelecInterpolation<domain_type, refDim>::value
			(*(m_spGF->domain().get()), elem, vCornerCoords, &(dofValues[0]), vLocIP, nip, vValue);
		
	//	Multiply the values by \f$ \sigma \f$
		for (size_t i = 0; i < nip; i++)
			vValue [i] *= sigma;
	};
};

/**
 * User data based class for the indicator function.
 * The subset names specified in the constructor define a subdomain \f$ S \f$
 * of the whole domain \f$ \Omega \f$. This user data class computes the
 * indicator function that is \f$ 1 \f$ in \f$ S \f$ and zero in
 * \f$ \Omega \setminus S \f$.
 */
template <typename TDomain>
class SubsetIndicatorUserData
	: public StdUserData<SubsetIndicatorUserData<TDomain>, number, TDomain::dim, void, UserData<number, TDomain::dim, void> >
{
public:
///	Type of domain
	typedef TDomain domain_type;
	
///	World dimension
	static const int dim = domain_type::dim;
	
/// subset handler type
	typedef typename domain_type::subset_handler_type subset_handler_type;

private:
	/// subset group representing the specified subdomain
	SubsetGroup m_ssGrp;
	
public:

///	Constructor
	SubsetIndicatorUserData
	(
		ConstSmartPtr<domain_type> domain, ///< domain of the problem
		const char * ss_names ///< subset names
	)
	: m_ssGrp (domain->subset_handler ())
	{
	// Parse the subset names:
		std::vector<std::string> vssNames;
		try
		{
			TokenizeString (ss_names, vssNames);
			for (size_t k = 0; k < vssNames.size (); k++)
				RemoveWhitespaceFromString (vssNames [k]);
			m_ssGrp.clear ();
			m_ssGrp.add (vssNames);
		} UG_CATCH_THROW ("SubsetIndicatorUserData: Failed to parse subset names.");
	}

///	Indicator functions are discontinuous
	virtual bool continuous () const {return false;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return true;}

///	Evaluator
	template <int refDim>
	inline void evaluate
	(
		number vValue [],
		const MathVector<dim> vGlobIP [],
		number time,
		int si,
		GeometricObject * elem,
		const MathVector<dim> vCornerCoords [],
		const MathVector<refDim> vLocIP [],
		const size_t nip,
		LocalVector * u,
		const MathMatrix<refDim, dim> * vJT = NULL
	) const
	{
	//	Get the subset index of the element
		int elem_si = m_ssGrp.subset_handler()->get_subset_index (elem);
	//	Check if the element is in one of the specified subsets:
		number indicator = (m_ssGrp.contains (elem_si))? 1 : 0;
	//	Return the indicator:
		for (size_t i = 0; i < nip; i++)
			vValue [i] = indicator;
	};
	
///	This function should not be used
	void operator()
	(
		number & vValue,
		const MathVector<dim> & globIP,
		number time,
		int si
	)
	const
	{
		UG_THROW("SubsetIndicatorUserData: Element required for evaluation, but not passed. Cannot evaluate.");
	}

///	This function should not be used
	void operator()
	(
		number vValue [],
		const MathVector<dim> vGlobIP [],
		number time,
		int si,
		const size_t nip
	) const
	{
		UG_THROW("SubsetIndicatorUserData: Element required for evaluation, but not passed. Cannot evaluate.");
	}
};

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_GF_USER_DATA__

/* End of File */
