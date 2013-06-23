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

// Nedelec-type-1 headers:
#include "nedelec_local_ass.h"

namespace ug{
namespace Electromagnetism{

/**
 * UserData base class that computes the vector values of grid functions
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
	// grid function
		SmartPtr<TGridFunc> m_spGF;

	//	component of function
		size_t m_fct;

public:
	/// constructor
		NedelecGridFunctionData
		(
			SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
			const char* cmp ///< the component of the grid function keeping the scalar DoFs
		)
		: m_spGF (spGridFct)
		{
		//	Get function id of name and check if the function exists:
			if((m_fct = m_spGF->fct_id_by_name (cmp)) >= m_spGF->num_fct())
				UG_THROW ("NedelecGridFunctionData: Function space does not contain"
						" a function with name " << cmp << ".");

		//	Get local finite element id and check whether this is the Nedelec element:
			if (m_spGF->local_finite_element_id (m_fct).type() != LFEID::NEDELEC)
				UG_THROW ("NedelecGridFunctionData: Function " << cmp
						<< "is not based on the Nedelec element.");
		};
		
	///	The vector field retrieved from the Nedelec-type 1 (Whitney-1) dofs are not continuous
		virtual bool continuous () const {return false;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct () const {return true;}

	/// Performs the main computations:
		template <int refDim>
		void eval_and_deriv
		(
			MathVector<dim> vValue[], ///< for the computed value
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
			std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		    const MathMatrix<refDim, dim>* vJT = NULL
		) const
		{
		//	Derivatives are not implemented
			if (bDeriv)
				UG_THROW ("NedelecGridFunctionData: Derivatives are not implemented.");
			
		//	Get multiindices of element
			std::vector<MultiIndex<2> > ind;
			m_spGF->multi_indices (elem, m_fct, ind);
		
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
 * UserData base class that computes the curl vector of grid functions
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
	// grid function
		SmartPtr<TGridFunc> m_spGF;

	//	component of function
		size_t m_fct;

public:
	/// constructor
		NedelecCurlData
		(
			SmartPtr<TGridFunc> spGridFct, ///< grid function with the DoFs
			const char* cmp ///< the component of the grid function keeping the scalar DoFs
		)
		: m_spGF (spGridFct)
		{
		//	Get function id of name and check if the function exists:
			if((m_fct = m_spGF->fct_id_by_name (cmp)) >= m_spGF->num_fct())
				UG_THROW ("NedelecCurlData: Function space does not contain"
						" a function with name " << cmp << ".");

		//	Get local finite element id and check whether this is the Nedelec element:
			if (m_spGF->local_finite_element_id (m_fct).type() != LFEID::NEDELEC)
				UG_THROW ("NedelecCurlData: Function " << cmp
						<< "is not based on the Nedelec element.");
		};
		
	///	The vector field retrieved from the Nedelec-type 1 (Whitney-1) dofs are not continuous
		virtual bool continuous () const {return false;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct () const {return true;}

	/// Performs the main computations:
		template <int refDim>
		void eval_and_deriv
		(
			MathVector<dim> vValue[], ///< for the computed value
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
			std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		    const MathMatrix<refDim, dim>* vJT = NULL
		) const
		{
		//	Derivatives are not implemented
			if (bDeriv)
				UG_THROW ("NedelecCurlData: Derivatives are not implemented.");
			
		//	Get multiindices of element
			std::vector<MultiIndex<2> > ind;
			m_spGF->multi_indices (elem, m_fct, ind);
		
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

} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_GF_USER_DATA__

/* End of File */
