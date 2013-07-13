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

/**
 * Computes the connectivity of the conductors and therefore finds out the
 * topology of the domain.
 *
 * \param[out] minCondInd	an array indexed by the subset indices, so that
 *							every entry of a conductor is initialized by the
 *							smallest conductor subset id connected to it
 *							(and by -1 for the other subdomains including insulators)
 */
template <typename TDomain>
void EMaterial<TDomain>::connectivity
(
	std::vector<int> & minCondInd
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<dim>::DimElemList ElemList;
	
//	Initialize the array:
	minCondInd.resize (subset_handler()->num_subsets ());
	
	for (size_t si = 0; si < minCondInd.size (); si++)
	{
		typename t_data_map::iterator iter = m_mUserDataBC.find (si);
		if (iter == m_mUserDataBC.end () || iter->second == NULL) // if no data
			minCondInd [si] = -2; // skip it
		else if (iter->second->sigma == (number) 0) // if insulator
			minCondInd [si] = -1; // skip it, too
		else
			minCondInd [si] = si; // this is the current minimum
	}

//	Call all the instances of 'get_elem_connectivity':
	boost::mpl::for_each<ElemList> (GetElemConnectivity (this, minCondInd));
};

/// An auxiliary class for the grid that "describes" the geometry
/**
 * This class gets the access to the elements that describe the basic
 * topology of the domain, for ex., to the coarsest grid in the grid
 * hierarcy. The general template definition provides only the error
 * messages. Cf. the specializations below.
 */
///\{
template <typename TElem, typename TSubsetHandler>
class TopologyDesc
{
public:
	/// Constructor
	TopologyDesc
	(
		ConstSmartPtr<TSubsetHandler> pSH ///< the subset handler
	)
	{
		UG_THROW ("Initialization of a TopologyDesc for an unknown SubsetHandler type.");
	}
	
	/// returns a pointer to the subset handler
	const GridSubsetHandler * subset_handler () const
	{
		UG_THROW ("Attempt to use a TopologyDesc for an unknown SubsetHandler type.");
	}
	
	///	returns a pointer to the grid on which the subset-handler works
	Grid * grid () const
	{
		UG_THROW ("Attempt to use a TopologyDesc for an unknown SubsetHandler type.");
	}

	///	returns the begin-iterator for the elements of type TElem in the given subset.
	typename geometry_traits<TElem>::const_iterator begin (int subsetIndex) const
	{
		UG_THROW ("Attempt to use a TopologyDesc for an unknown SubsetHandler type.");
	}
	
	///	returns the end-iterator for the elements of type TElem in the given subset.
	typename geometry_traits<TElem>::const_iterator end (int subsetIndex) const
	{
		UG_THROW ("Attempt to use a TopologyDesc for an unknown SubsetHandler type.");
	}
};

template <typename TElem>
class TopologyDesc <TElem, GridSubsetHandler>
{
public:
	/// Constructor
	TopologyDesc
	(
		ConstSmartPtr<GridSubsetHandler> pSH ///< the subset handler
	)
	:	m_pSH (pSH)
	{
		if (m_pSH.invalid ())
			UG_THROW ("TopologyDesc: specify a valid subset handler.");
	}
	
	/// returns a pointer to the subset handler
	const GridSubsetHandler * subset_handler () const {return m_pSH.get ();}
	
	///	returns a pointer to the grid on which the subset-handler works
	Grid * grid () const
	{
		return m_pSH->grid ();
	}

	///	returns the begin-iterator for the elements of type TElem in the given subset.
	typename geometry_traits<TElem>::const_iterator begin (int subsetIndex) const
	{
		return m_pSH->template begin<TElem> (subsetIndex);
	}
	
	///	returns the end-iterator for the elements of type TElem in the given subset.
	typename geometry_traits<TElem>::const_iterator end (int subsetIndex) const
	{
		return m_pSH->template end<TElem> (subsetIndex);
	}
	
private:
	ConstSmartPtr<GridSubsetHandler> m_pSH; ///< the subset handler
};

template <typename TElem>
class TopologyDesc <TElem, MultiGridSubsetHandler>
{
public:
	/// Constructor
	TopologyDesc
	(
		ConstSmartPtr<MultiGridSubsetHandler> pSH ///< the subset handler
	)
	:	m_pSH (pSH)
	{
		if (m_pSH.invalid ())
			UG_THROW ("TopologyDesc: specify a valid subset handler.");
		if (! ((MultiGridSubsetHandler *) m_pSH.get ())->subset_inheritance_enabled ())
			UG_THROW ("TopologyDesc: Subset inheritance should be enabled in the subset handler.");
	}
	
	/// returns a pointer to the subset handler
	const MultiGridSubsetHandler * subset_handler () const {return m_pSH.get ();}
	
	///	returns a pointer to the grid on which the subset-handler works
	Grid * grid () const
	{
		return m_pSH->grid ();
	}

	///	returns the begin-iterator for the elements of type TElem in the given subset.
	typename geometry_traits<TElem>::const_iterator begin (int subsetIndex) const
	{
		return m_pSH->template begin<TElem> (subsetIndex, 0);
	}
	
	///	returns the end-iterator for the elements of type TElem in the given subset.
	typename geometry_traits<TElem>::const_iterator end (int subsetIndex) const
	{
		return m_pSH->template end<TElem> (subsetIndex, 0);
	}
	
private:
	ConstSmartPtr<MultiGridSubsetHandler> m_pSH; ///< the subset handler
};
///\}

/**
 * Helper for 'connectivity'. The sematics of the argument is the same as for
 * 'connectivity'. Note that the proper size of minCondInd should be set before
 * the call. Furthermore, this array is updated, not reset.
 */
template <typename TDomain>
template <typename TElem>
void EMaterial<TDomain>::get_elem_connectivity
(
	std::vector<int> & minCondInd
)
{
	typedef typename geometry_traits<TElem>::geometric_base_object base_object;
	typedef typename geometry_traits<TElem>::const_iterator elem_iterator;
	
	UG_ASSERT (((int) minCondInd.size ()) == subset_handler()->num_subsets (), "get_elem_connectivity: array size mismatch");
	
	std::vector<base_object*> neighbours;
	
//	Access to the elements:
	TopologyDesc<TElem, subset_handler_type> base_grid (subset_handler ());
	
//	Loop over the subsets:
	for (size_t si = 0; si < minCondInd.size (); si++)
	{
		int min_si;
		
	//	Conductor?
		if ((min_si = minCondInd [si]) < 0)
			continue; // no, we do not treat this subset
		
	//	Yes, loop over the elements in the subdomain:
		elem_iterator e_end = base_grid.end (si);
		for (elem_iterator e_iter = base_grid.begin (si); e_iter != e_end; ++e_iter)
		{
		//	Loop over the neighbours:
			CollectNeighbors (neighbours, *e_iter, *base_grid.grid(), NHT_VERTEX_NEIGHBORS);
			for (size_t k = 0; k < neighbours.size (); k++)
			{
				int min_nbr_si;
				int nbr_si = base_grid.subset_handler()->get_subset_index (neighbours [k]);
				
				if (nbr_si < 0 || nbr_si >= (int) minCondInd.size ())
					UG_THROW ("get_elem_connectivity: Illegal neighbour subset index.");
				if ((min_nbr_si = minCondInd [nbr_si]) < 0)
					continue; // we do not treat this subset
				
				if (min_nbr_si < min_si)
					for (size_t l = 0; l < minCondInd.size (); l++)
						if (minCondInd [l] == min_si)
							minCondInd [l] = min_nbr_si;
			}
		}
	}
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
