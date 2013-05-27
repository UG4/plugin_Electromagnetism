/**
 * Base classes for problem settings for the Maxwell equations.
 *
 * Created on: 13.03.2013
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__EM_DATA__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__EM_DATA__

#include <vector>
#include <map>

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/subset_util.h"

#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

namespace ug{
namespace Electromagnetism{

/// Class for subdomain-dependent data for the E-based formulated problems.
/**
 * This class stores the magnetic permeabilities and electric conductivities
 * for the subdomains. It is used as a parameter set for the discretizations.
 * Furthermore, this class finds out the connectivity of the conductors.
 * This information is used (in particular) for the projection of the
 * solution to the divergence-free space.
 */
template <typename TDomain>
class EMaterial
{
private:
/// own type
	typedef EMaterial<TDomain> this_type;
	
///	domain type
	typedef TDomain domain_type;

/// subset handler type
	typedef typename domain_type::subset_handler_type subset_handler_type;
	
///	world dimension
	static const int dim = domain_type::dim;

public:
/// Constructor
	EMaterial
	(
		ConstSmartPtr<domain_type> domain ///< [in] domain of the problem
	);
	
/// adds a generic subset data item
	void add
	(
		const char * subsets, ///< [in] names of the subsets
		number mu, ///< [in] (nonzero) magnetic permeability
		number sigma ///< [in] electric conductivity
	);

/// adds a insulator
	void add
	(
		const char * subsets, ///< [in] names of the subsets
		number mu ///< [in] (nonzero) magnetic permeability
	)
	{
		add (subsets, mu, (number) 0);
	}

/// finalizes the object
	void close ();
	
/// returns the string of the subset names
	const char * subset_names () const
	{
		if (! m_bClosed)
			UG_THROW ("EMaterial: The object has not been closed.");
		return m_sSsNames.c_str ();
	}

/// returns true iff closed
	bool finalized () const {return m_bClosed;}
	
/// returns pointer to the domain
	ConstSmartPtr<domain_type> domain () const {return m_spDomain;}
	
///	reads the data for a subdomain from the data item
	bool get_mu_sigma
	(
		int si, ///< [in] subset index
		number& mu, ///< [out] magnetic permeability
		number& sigma ///< [out] electric conductivity
	) const;

/** 
 * The connectivity of the conductors:
 * This array is indexed by the subset indices, so that
 * every entry of a conductor is initialized by the
 * smallest conductor subset id connected to it
 * (and by -1 for insulators as well as <= -2 the other subdomains)
 */
	const std::vector<int> & min_conductor_index () const
	{
		return m_minCondInd;
	}
	
/**
 * Base conductor indices: The list of the smallest indices
 * of the subsets occupied by conductors. The length of this
 * list is the 2nd Betti number: the number of the insulator
 * separated subdomains occupied by the conductors. Note that
 * the Dirichlet BC is not taken into account here.
 */
	const std::vector<int> & base_conductors () const
	{
		return m_baseConductors;
	}
	
	/// constant access to the subset handler
	ConstSmartPtr<subset_handler_type> subset_handler () const {return m_spDomain->subset_handler ();};

private:
	
	/// computes the connectivity of the conductions
	void connectivity
	(
		std::vector<int> & minCondInd ///< [out] min. subset ind. of the conductor connected to a given one (or -1)
	);
	
	/// a helper function for 'connectivity'
	template <typename TElem>
	void get_elem_connectivity
	(
		std::vector<int> & minCondInd ///< [out] min. subset ind. of the conductor connected to a given one (or -1)
	);
	
	/// a helper class for 'connectivity'
	struct GetElemConnectivity
	{
		this_type * m_pThis;
		std::vector<int> & m_minCondInd;
		GetElemConnectivity (this_type* pThis, std::vector<int> & minCondInd)
		 : m_pThis (pThis), m_minCondInd (minCondInd) {}
		template <typename TElem> void operator() (TElem &)
			{m_pThis->get_elem_connectivity<TElem> (m_minCondInd);}
	};
	
	/// analyzes the conductor topology of the domain
	void analyze_topology ();

private:
	/// domain
	ConstSmartPtr<domain_type> m_spDomain;
	
	///	data item type
	struct TSubdomData
	{
		std::string ssNames; ///< subset names
		SubsetGroup ssGrp; ///< subset group
		
		number mu; ///< magnetic permeability
		number sigma; ///< electric conductivity
		
		/// Constructor:
		TSubdomData
		(
			ConstSmartPtr<subset_handler_type> pSH,
			const char * names, // the names
			number the_mu, // magnetic permeability
			number the_sigma // electric conductivity
		)
		:	ssNames (names), ssGrp (pSH),
			mu (the_mu), sigma (the_sigma)
		{};
	};

	/// Subdomain data items
	std::vector<TSubdomData> m_vSdD;
	
	///	Data map type
	typedef std::map<int, TSubdomData *> t_data_map;

	/// Map assigning subdomain indices to the subdomain data items
	t_data_map m_mUserDataBC;
	
	/// Flag that indicates that the description of the domain has been completed
	bool m_bClosed;
	
	/// String of all the subset names mentioned in the data items
	std::string m_sSsNames;
	
	/** 
	 * The connectivity of the conductors:
	 * This array is indexed by the subset indices, so that
	 * every entry of a conductor is initialized by the
	 * smallest conductor subset id connected to it
	 * (and by -1 for insulators as well as <= -2 for the other subdomains)
	 */
	std::vector<int> m_minCondInd;
	
	/**
	 * Base conductor indices: The list of the smalles indices
	 * of the subsets occupied by conductors. The length of this
	 * list is the 2nd Betti number: the number of the insulator
	 * separated subdomains occupied by the conductors. Note that
	 * the Dirichlet BC is not taken into account here.
	 */
	std::vector<int> m_baseConductors;
};

/// Common interface to get the Dirichlet boundary conditions
/**
 * The class providing a common interface to get low-dimensional
 * subsets with Dirichlet boundary conditions.
 */
template <typename TDomain, typename TAlgebra>
class EMDirichlet
	: public IDomainConstraint<TDomain, TAlgebra>
{
public:
	
	/// should extend the given subset group with the Dirichlet subsets
	virtual void get_dirichlet_subsets
	(
		SubsetGroup & dirichlet_ssgrp ///< the group to update
	) const = 0;
};

} // end namespace Electromagnetism
} // end namespace ug

#include "em_material_impl.h"

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__EM_DATA__

/* End of File */
