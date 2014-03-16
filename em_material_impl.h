/**
 * Implementation of base classes for problem settings for the Maxwell equations.
 *
 * Created on: 14.03.2013
 * Author: D. Logashenko
 */
#include "lib_disc/domain_traits.h"

namespace ug{
namespace Electromagnetism{

/**
 * Constructor: initializes the object and checks the validity of its parts
 */
template <typename TDomain>
EMaterial<TDomain>::EMaterial
(
	ConstSmartPtr<domain_type> domain ///< [in] domain of the problem
)
:	m_spDomain (domain),
	m_bClosed (false)
{
	if (domain.invalid () || subset_handler().invalid ())
		UG_THROW ("EMaterial: Invalid domain or subset handler.");
}

/**
 * Gets the data from a subset data item. This function returns true if the
 * specified subset index is not listed in the object, or the required data
 * are not supplied for this subset. (For example, it should return true for
 * low-dimensional subsets.) Otherwise the function returns false.
 * If the object is not finalized, the function throws and exception.
 */
template <typename TDomain>
bool EMaterial<TDomain>::get_mu_sigma
(
	int si, ///< [in] subset index
	number& mu, ///< [out] magnetic permeability
	number& sigma ///< [out] electric conductivity
) const
{
	TSubdomData * pSdD;
	
	if (! m_bClosed)
		UG_THROW ("EMaterial: Attempt to get data from a unfinalized object.")
	
	typename t_data_map::const_iterator iter = m_mUserDataBC.find (si);
	if (iter == m_mUserDataBC.end () || (pSdD = iter->second) == NULL)
		return true;
	
	mu = pSdD->mu;
	sigma = pSdD->sigma;
	return false;
}
	
/**
 * Adds a data item to the object.
 */
template <typename TDomain>
void EMaterial<TDomain>::add
(
	const char * subsets, ///< [in] names of the subsets
	number mu, ///< [in] (nonzero) magnetic permeability
	number sigma ///< [in] electric conductivity
)
{
//	Check if already finalized
	if (m_bClosed)
		UG_THROW ("EMaterial::add:"
			" Attempt to add a subset data item to a finalized domain description.");
	
//	Add a new element
	m_vSdD.push_back (TSubdomData (subset_handler (), subsets, mu, sigma));
}

/**
 * Performs finalization steps for the object.
 */
template <typename TDomain>
void EMaterial<TDomain>::close ()
{
	std::vector<std::string> vssNames;
	
// Clear the map and remove all the subset indices:
	m_bClosed = false;
	m_mUserDataBC.clear ();
	m_sSsNames = "";
	
	if (m_vSdD.size () == 0)
		UG_THROW ("No data items specified.");
	
// Parse the subset names:
	try
	{
		for (size_t i = 0; i < m_vSdD.size (); i++)
		{
			TSubdomData & sdD = m_vSdD [i];
			TokenizeString (sdD.ssNames, vssNames);
			for (size_t k = 0; k < vssNames.size (); k++)
				RemoveWhitespaceFromString (vssNames [k]);
			sdD.ssGrp.clear ();
			sdD.ssGrp.add (vssNames);
		}
	} UG_CATCH_THROW ("EMaterial::close: Failed to parse subset names.");
	
// Get the subset handler:
	const subset_handler_type * ss_handler = subset_handler().get ();

// Fill the map with the appropriate subsets:
	for (int si = 0; si < ss_handler->num_subsets (); si++)
		if (DimensionOfSubset (*ss_handler, si) == dim) // do not consider low-dimensional subsets
			m_mUserDataBC [si] = 0;
	
// Fill the map with the data items:
	for (size_t i = 0; i < m_vSdD.size (); i++)
	{
		TSubdomData * psdD = & (m_vSdD [i]);
		const SubsetGroup & ssg = psdD->ssGrp;
		for (size_t k = 0; k < ssg.size (); k++)
		{
			typename t_data_map::iterator iter = m_mUserDataBC.find (ssg [k]);
			if (iter == m_mUserDataBC.end ())
				UG_THROW ("EMaterial::close: Refered subset " << ssg.name (k) << "has an illegal dimension.");
			if (iter->second != NULL)
				UG_THROW ("EMaterial::close: Two data items refer to subset " << ssg.name (k) << ".");
			iter->second = psdD;
		}
	}
	
// Check if there are full-dim. subsets that have not been mentioned:
	for (typename t_data_map::iterator iter = m_mUserDataBC.begin ();
		 iter != m_mUserDataBC.end (); ++iter)
		if (iter->second == NULL)
			UG_THROW ("EMaterial::close: Subset "
				<< ss_handler->get_subset_name (iter->first)
					<< " not mentioned in the description.");
	
// Analyze the topology
	analyze_topology ();
	
// Compose the string of all the names
	m_sSsNames = m_vSdD[0].ssNames;
	for (size_t i = 1; i < m_vSdD.size (); i++)
	{
		m_sSsNames += ',';
		m_sSsNames += m_vSdD[i].ssNames;
	}
	
// Mark the object as finalized:
	m_bClosed = true;
	
// Print the 'greeting':
	UG_LOG ("Materials specified. " << m_baseConductors.size () << " conductive parts found.\n");
}

/**
 * Computes the connectivity of the conductors and therefore finds out the
 * topology of the domain.
 *
 * \param[out] minCondInd	an array indexed by the subset indices, so that
 *							every entry of a conductor is initialized with the
 *							smallest conductor subset id connected to it
 *							(and with < 0 for the other subdomains including insulators)
 */
template <typename TDomain>
void EMaterial<TDomain>::connectivity
(
	std::vector<int> & minCondInd
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<dim>::grid_base_object t_base_object;
	
//	Initialize the marks of the conductor subsets:
	std::vector<bool> isConductor (subset_handler()->num_subsets ());
	
	for (size_t si = 0; si < isConductor.size (); si++)
	{
		typename t_data_map::iterator iter = m_mUserDataBC.find (si);
		if (iter == m_mUserDataBC.end () || iter->second == NULL) // if no data
			isConductor [si] = false; // skip it
		else if (iter->second->sigma == (number) 0) // if insulator
			isConductor [si] = false; // skip it, too
		else
			isConductor [si] = true; // this one should be considered
	}

//	Find out the subset connectivity:
	FindSubsetGroups<t_base_object> (minCondInd, isConductor, * subset_handler().get(),
		NHT_VERTEX_NEIGHBORS);
	
//	Note that the low-dimensional subsets have been "not-marked" and got '-1'
//	in minCondInd. As '-1' should denote the insulators only, we revise the 
//	array. Furthermore, there should be no marks '-2' at all because none of
//	the low-dimensional subsets have been marked. We check it to be on the
//	safe side:
	for (size_t si = 0; si < minCondInd.size (); si++)
	if (minCondInd [si] == -1) // check whether this is an insulator
	{
		if (m_mUserDataBC.find (si) == m_mUserDataBC.end ())
			minCondInd [si] = -2; // this is no insulator; otherwise, there would be data for it
	}
	else if (minCondInd [si] < -1) // check for marked low-dimensional subdomains
		UG_THROW ("EMaterial::connectivity: Low-dimensional conductor found (subset index " << si << ").");
};

/**
 * Analyses the conductor topology of the domain.
 */
template <typename TDomain>
void EMaterial<TDomain>::analyze_topology ()
{
//	Get the connectivity:
	connectivity (m_minCondSsI);

//	Compose the list of the minimum conductor subset indices:
	std::vector<bool> is_min_index (m_minCondSsI.size ());
	for (size_t i = 0; i < m_minCondSsI.size (); i++)
	{
		is_min_index [i] = false; // note that (m_minCondSsI[i] <= i)
		if (m_minCondSsI[i] >= 0)
			is_min_index [m_minCondSsI[i]] = true;
	}
	m_baseConductors.clear ();
	for (size_t i = 0; i < is_min_index.size (); i++)
		if (is_min_index [i])
			m_baseConductors.push_back (i);
	
//	Compose the references to the base conductors
	m_baseCondInd.resize (m_minCondSsI.size ());
	for (size_t i = 0; i < m_minCondSsI.size (); i++)
	{
		int min_cond_ssi = m_minCondSsI [i];
		if (min_cond_ssi <= -2)
			m_baseCondInd [i] = min_cond_ssi;
		else
			m_baseCondInd [i] = -1; // to be on the safe side
	}
	for (size_t base_cond = 0; base_cond < m_baseConductors.size (); base_cond++)
	{
		int base_cond_si = m_baseConductors [base_cond];
		for (size_t i = 0; i < m_minCondSsI.size (); i++)
			if (m_minCondSsI [i] == base_cond_si)
				m_baseCondInd [i] = base_cond;
	}
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
