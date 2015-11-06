/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
 * Base classes for problem settings for the Maxwell equations.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__EM_DATA__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__EM_DATA__

#include <vector>
#include <map>

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_grid/tools/subset_group.h"
#include "lib_grid/algorithms/subset_dim_util.h"

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
	
///	constant access to the subset handler
	ConstSmartPtr<subset_handler_type> subset_handler () const
	{
		return m_spDomain->subset_handler ();
	}

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
 * This function returns the array indexed by the subset indices, so that
 * to every subset index of a conductor, it assignes the smallest conductor
 * subset index connected to it (and -1 to insulators as well as <= -2 to
 * the other subdomains)
 */
	const std::vector<int> & min_conductor_ssi () const
	{
		return m_minCondSsI;
	}
	
/** 
 * For every subset index of a conductor, this function returns the smallest
 * conductor subset index connected to the given conductor (and -1 for
 * insulators as well as <= -2 for the other subdomains)
 */
	int min_conductor_ssi
	(
		int si ///< the subset index of the conductor
	) const
	{
		return m_minCondSsI [si];
	}
	
/**
 * This function returns the array of the base conductor subset indices:
 * The list of the smallest indices of the subsets representing conductors.
 * The length of this array is the 2nd Betti number: the number of the
 * insulator-separated subdomains occupied by the conductors. Note that
 * the Dirichlet BC is not taken into account here and some of the conductors
 * may be grounded.
 */
	const std::vector<int> & base_conductors () const
	{
		return m_baseConductors;
	}
	
/**
 * This function returns the array that assignes to every subset index of
 * a conductor the index of its base conductor in the list of based
 * conductors (as returned by base_conductors). Thus, for a conductor
 * in subset si, min_conductor_ssi (si) == base_conductors ()
 * [base_conductor_index () [si]]. If the subset represents an
 * insulator, the array assignes -1 to it. To other subsets, the
 * array assignes <= -2.
 */
	const std::vector<int> & base_conductor_index () const
	{
		return m_baseCondInd;
	}
	
/**
 * Returns index of the base conductor in base_conductors () for a conductor
 * identified by the subset index of its subset. If the subset represents an
 * insulator, the function returns -1. For other subsets, <= -2 is returned.
 */
	int base_conductor_index
	(
		int si ///< the subset index of the conductor
	) const
	{
		return m_baseCondInd [si];
	}
	
private:
	
	/// computes the connectivity of the conductions
	void connectivity
	(
		std::vector<int> & minCondInd ///< [out] min. subset ind. of the conductor connected to a given one (or -1)
	);
	
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
	 * This array is indexed by the subset indices and assignes to
	 * every subset of a conductor the smallest conductor subset index
	 * connected to this conductor (and -1 to insulators as well
	 * as <= -2 to the other subdomains)
	 */
	std::vector<int> m_minCondSsI;
	
	/**
	 * Base conductor subset indices: The list of the smallest indices
	 * of the subsets representing conductors. The length of this
	 * list is the 2nd Betti number: the number of the insulator-
	 * separated subdomains occupied by the conductors. Note that
	 * the Dirichlet BC is not taken into account here so that some
	 * of the base conductors may be grounded.
	 */
	std::vector<int> m_baseConductors;
	
	/**
	 * Correspondence of subset indices to the indices of the base
	 * conductors in m_baseConductors: For every conductor with subset
	 * index si, m_minCondSsI [si] == m_baseConductors [m_baseCondInd [si]].
	 * For other subsets, this array keeps either -1 (for insulators) or
	 * <= -2 (for low-dimensional subsets etc).
	 */
	std::vector<int> m_baseCondInd;
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
