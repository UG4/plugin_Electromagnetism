/**
 * nedelec_source_impl.h - implementation of class members for computation of
 * divergence-free sources in form of Nedelec-type-1 element DoFs.
 *
 * Created on: Jun. 30, 2013
 * Author: D. Logashenko
 */
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

#include "nedelec_local_ass.h"

namespace ug{
namespace Electromagnetism{

/**
 * Class constructor:
 */
template <typename TDomain, typename TAlgebra>
NedelecLoopCurrent<TDomain, TAlgebra>::NedelecLoopCurrent
(
	const char * ssNames, ///< names of the subsets of the source (up to the subset of the pos. dir.)
	const char * posSsNames, ///< names of the subsets of the positive direction
	const char * cutSsNames, ///< names of the surfaces on the cut of the loop
	SmartPtr<ApproximationSpace<TDomain> > vertApproxSpace, ///< vertex-centered approx. space
	SmartPtr<ILinearOperatorInverse<pot_vector_type> > potSolver ///< linear solver for the potential
)
:	m_allSsNames ((std::string (ssNames) + ',') + posSsNames), m_posSsNames (posSsNames), m_cutSsNames (cutSsNames),
	m_spVertApproxSpace (vertApproxSpace),
	m_auxLocLaplace (new AuxLaplaceLocAss (*this)),
	m_outOfSource (new OutOfSource (*this)),
	m_auxLaplaceAss (new DomainDiscretization<TDomain, TPotAlgebra> (vertApproxSpace)),
	m_auxLaplaceOp (new AssembledLinearOperator<TPotAlgebra> (SmartPtr<IAssemble<TPotAlgebra> >(m_auxLaplaceAss))),
	m_potSolver (potSolver)
{
//	Check the parameters:
	if (m_spVertApproxSpace.invalid ())
		UG_THROW ("NedelecLoopCurrent: Illegal vert.-centered approx. space.");
	if (m_spVertApproxSpace->num_fct () != 1)
		UG_THROW ("NedelecLoopCurrent: Exactly one function should be defined in the vert.-centered approx. space.");
	if (! m_spVertApproxSpace->is_def_everywhere (0))
		UG_THROW ("NedelecLoopCurrent: The function in the vert.-centered approx. space must be defined everywhere.");
	if (m_potSolver.invalid ())
		UG_THROW ("NedelecLoopCurrent: Illegal solver for the auxiliary problems.");
	
//	Fill the subset groups:
	std::vector<std::string> vssNames;
	ConstSmartPtr<subset_handler_type> spIsh = vertApproxSpace->subset_handler ();
	
	TokenizeString (m_allSsNames, vssNames);
	for (size_t i = 0; i < vssNames.size(); i++)
		RemoveWhitespaceFromString (vssNames [i]);
	m_allSsGrp.set_subset_handler (spIsh); m_allSsGrp.add (vssNames);
	
	TokenizeString (m_posSsNames, vssNames);
	for (size_t i = 0; i < vssNames.size(); i++)
		RemoveWhitespaceFromString (vssNames [i]);
	m_posSsGrp.set_subset_handler (spIsh); m_posSsGrp.add (vssNames);
	
	TokenizeString (m_cutSsNames, vssNames);
	for (size_t i = 0; i < vssNames.size(); i++)
		RemoveWhitespaceFromString (vssNames [i]);
	m_cutSsGrp.set_subset_handler (spIsh); m_cutSsGrp.add (vssNames);
	
//	Compose the global discretization of the auxiliary equations:
	m_auxLaplaceAss->add (SmartPtr<IElemDisc<TDomain> >(m_auxLocLaplace));
	m_auxLaplaceAss->add
		(SmartPtr<IDomainConstraint<TDomain, TPotAlgebra> >(m_outOfSource));
}

/**
 * Setting the electric current value 
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::set
(
	const char * fctNames, ///< names of the components
	number I ///< electric current for the functions
)
{
	m_vSrcData.push_back (TSrcData (fctNames, I));
}

/**
 * Computation of the source by updating a given grid function
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::compute
(
	SmartPtr<GridFunction<TDomain, TAlgebra> > sp_u ///< the grid function for the source
)
{
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	
//	Do we have data to compute?
	if (m_vSrcData.size () == 0)
		UG_THROW ("NedelecLoopCurrent: No electric currents specified.");
	
//	Check the grid function:
	if (sp_u.invalid ())
		UG_THROW ("NedelecLoopCurrent: Illegal grid function specification.");
	
//	Get the DoF distributions:
	GridFunction<TDomain, TAlgebra> & u = * sp_u.get ();
	SmartPtr<DoFDistribution> edgeDD = u.dof_distribution ();
	const GridLevel g_lev (u.grid_level ());
	SmartPtr<DoFDistribution> vertDD = m_spVertApproxSpace->dof_distribution (g_lev);
	
//	Compute the potential of the source:
	pot_vector_type pot_u;
	compute_potential (* vertDD.get (), g_lev, pot_u);
	
//	Compute the normalization factor of the potential (to scale the current to 1)
	boost::mpl::for_each<ElemList>
		(GetFluxOfPotential (this, * m_spVertApproxSpace->domain().get (), pot_u,
			* vertDD.get (), m_pot_scaling));
	m_pot_scaling = - m_pot_scaling;
	
//	Loop over the source data
	for (size_t i_data = 0; i_data < m_vSrcData.size (); i_data++)
	{
	//	Get the component indices:
		FunctionGroup fctGrp;
		try
		{
			fctGrp = u.fct_grp_by_name (m_vSrcData[i_data].fctNames.c_str ());
		}
		UG_CATCH_THROW ("NedelecLoopCurrent: Functions '" << m_vSrcData[i_data].fctNames <<
			"' not all contained in the edge approximation space.");
		
	//	Check the functions
		for (size_t i_fct = 0; i_fct < fctGrp.size (); i_fct++)
			if (u.local_finite_element_id(fctGrp[i_fct]).type () != LFEID::NEDELEC)
				UG_THROW ("NedelecLoopCurrent: Not a Nedelec-element-based grid function specified for the source.");
		
	//	Compute the gradients of the potential
		number value = m_vSrcData[i_data].I;
		for (size_t i_fct = 0; i_fct < fctGrp.size (); i_fct++)
			distribute_source_potential (* vertDD.get (), pot_u, * edgeDD.get (),
				fctGrp[i_fct], value, u);
	}
}

/**
 * Marks edges that belong to the loop source domain (for one type of elements).
 * For edges of elements belonging to the source domain, the entries of
 * 'in_source' are set to non-zero. (Other entries are not set.) For those
 * entries that correspond to edges belonging to elements in the 'positive
 * direction' subdomain and having one of the ends at the cut surface, the
 * 1st and the 2nd bits of the entry is set to 1 (depending whether the
 * beginning or the end lyies at the cut). For all other edges in the source,
 * these bits are set to zero (so that only the 0th bit is 1).
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::mark_source_edges
(
	const DoFDistribution & edgeDD, ///< [in] the edge DD
	VariableArray1<char> & in_source ///< [out] the arrays of flags to update
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
	const ISubsetHandler * pIsh = edgeDD.subset_handler().get ();
	std::vector<size_t> vEdgeInd (1);
	
//	Loop over the source subsets:
	for (size_t i = 0; i < m_allSsGrp.size (); i++)
	{
		int si = m_allSsGrp [i];
		bool in_pos = m_posSsGrp.contains (si);
	//	Loop over all the elements of the given type in the subset
		iterator e_end = edgeDD.template end<TElem> (si);
		for (iterator elem_iter = edgeDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
			TElem * pElem = *elem_iter;
			Grid::edge_traits::secure_container edge_list;
			pIsh->grid()->associated_elements_sorted (edge_list, pElem);
			UG_ASSERT ((edge_list.size () == (size_t) ref_elem_type::numEdges),
				"Mismatch of numbers of corners and vertices of an element");
			for (size_t e = 0; e < (size_t) ref_elem_type::numEdges; e++)
			{
				EdgeBase * pEdge = edge_list[e];
				char flag = 1;
				if (edgeDD.inner_algebra_indices (pEdge, vEdgeInd) != 1)
					UG_THROW ("NedelecLoopCurrent: Edge DoF distribution mismatch. Not the Nedelec-type-1 element?");
				if (in_pos)
				{
					if (m_cutSsGrp.contains (pIsh->get_subset_index (pEdge->vertex (0))))
						flag |= 2;
					if (m_cutSsGrp.contains (pIsh->get_subset_index (pEdge->vertex (1))))
						flag |= 4;
				}
				in_source [vEdgeInd [0]] = flag;
			}
		}
	}
}

/**
 * Computes the potential of the source
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::compute_potential
(
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	const GridLevel & grid_lev, ///< [in] grid level of the vector to correct
	pot_vector_type & pot_u ///< [out] a vector for the potential
)
{
	pot_vector_type pot_rhs;
	
//	Resize the vectors for the auxiliary system:
	size_t aux_num_ind = vertDD.num_indices ();
	pot_rhs.resize (aux_num_ind);
	pot_u.resize (aux_num_ind);
	
//	Assemble the matrix of the auxiliary problem:
	pot_u.set (0.0);
	m_auxLaplaceOp->set_level (grid_lev);
	m_auxLaplaceOp->init_op_and_rhs (pot_rhs);
	
//	Initizlize the solver:
	m_potSolver->init (m_auxLaplaceOp);
	
//	Assemble and solve the Laplace equation:
	m_potSolver->apply (pot_u, pot_rhs);
}

/**
 * Computes the gradient of the potential for one function
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::distribute_source_potential
(
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	pot_vector_type & src_pot, ///< [in] potential to distribute
	const DoFDistribution & edgeDD, ///< [in] the edge DD
	size_t func, ///< [in] index of the function
	number value, ///< [in] value of the source
	vector_type & src_field ///< [out] the computed source field
)
{
	typedef DoFDistribution::traits<Edge>::const_iterator t_edge_iter;
	
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	
//	Mark the edges in the source:
	VariableArray1<char> in_source (edgeDD.num_indices ());
	memset (& (in_source.at (0)), 0, in_source.size () * sizeof (char));
	boost::mpl::for_each<ElemList> (MarkSourceEdges (this, edgeDD, in_source));
	
//	Scaling of the potential:
	value /= m_pot_scaling;
	
//	Compute the gradient:
	std::vector<size_t> vVertInd (1);
	std::vector<DoFIndex> vEdgeInd (1);
	t_edge_iter edgeIterEnd = edgeDD.end<Edge> ();
	for (t_edge_iter edgeIter = edgeDD.begin<Edge> (); edgeIter != edgeIterEnd; ++edgeIter)
	{
		EdgeBase * pEdge = *edgeIter;
		char edge_flag;
		number nd_pot [2];
		
		// Get the edge DoF and check whether the edge is in the source:
		if (edgeDD.inner_multi_indices (pEdge, func, vEdgeInd) != 1)
			UG_THROW ("NedelecLoopCurrent: Edge DoF distribution mismatch. Not the Nedelec-Type-1 element?");
		if ((edge_flag = in_source [vEdgeInd [0] [0]]) == 0)
		{ // we are not in the source
			DoFRef (src_field, vEdgeInd [0]) = 0;
			continue;
		}
		
		// Get the values at the ends:
		for (size_t i = 0; i < 2; i++)
		{
			VertexBase * pVertex = pEdge->vertex (i);
			if (vertDD.inner_algebra_indices (pVertex, vVertInd) != 1)
				UG_THROW ("NedelecLoopCurrent: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			nd_pot [i] = src_pot [vVertInd [0]];
			if ((edge_flag & (2 << i)) != 0)
				nd_pot [i] += 1.0; // the jump of the potential at the cut
		}
		
		// Compute the gradient:
		DoFRef (src_field, vEdgeInd [0]) = value * (nd_pot [1] - nd_pot [0]);
	}
}

/*----- Assembling the auxiliary Poisson problems 1: class AuxLaplaceLocAss -----*/

/**
 * Computes the local discretization of the Laplace operator
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::LocLaplaceA<TElem>::stiffness
(
	GeometricObject * elem, ///< [in] element to prepare
	const position_type vCornerCoords [], ///< [in] coordinates of the corners of the element
	number loc_A [numCorners] [numCorners] ///< [out] the local stiffness matrix
)
{
	typedef FEGeometry<TElem, WDim, LagrangeLSFS<ref_elem_type, 1>, GaussQuadrature<ref_elem_type, 1> > TFEGeom;
	
//	request the finite element geometry
	TFEGeom & geo = GeomProvider<TFEGeom>::get ();
	// we assume that this is the simplest discretization:
	UG_ASSERT (geo.num_ip () == 1, "Only the simplest quadrature is supported here.");
	
//	initialize the fe geometry
	geo.update (elem, vCornerCoords);
	
//	compute the upper triangle of the local stiffness matrix
	for (size_t from_co = 0; from_co < (size_t) numCorners; from_co++)
		for (size_t to_co = from_co; to_co < (size_t) numCorners; to_co++)
			loc_A [from_co] [to_co]
				= VecDot (geo.global_grad (0, from_co), geo.global_grad (0, to_co))
					* geo.weight (0);
//	use the symmetry
	for (size_t from_co = 1; from_co < (size_t) numCorners; from_co++)
		for (size_t to_co = 0; to_co < from_co; to_co++)
			loc_A [from_co] [to_co] = loc_A [to_co] [from_co];
}

/**
 * Checks whether corners are on boundary (returns true if yes)
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
bool NedelecLoopCurrent<TDomain, TAlgebra>::LocLaplaceA<TElem>::bnd_corners
(
	GeometricObject * elem, ///< [in] element to compute for
	const SubsetGroup & bndGrp, ///< [in] boundary subsets
	bool bnd_flag [numCorners] ///< [out] true for bnd nodes, false for others
)
{
	TElem * pElem = (TElem *) elem;
	const ISubsetHandler * pIsh = bndGrp.subset_handler().get ();
	bool bnd = false;
	
	for (size_t i = 0; i < (size_t) numCorners; i++)
		if ((bnd_flag [i] =  bndGrp.contains (pIsh->get_subset_index (pElem->vertex (i)))))
			bnd = true;
	return bnd;
}

/**
 * Class constructor that gets the master NedelecLoopCurrent object and extracts
 * the data necessary for the base class:
 */
template <typename TDomain, typename TAlgebra>
NedelecLoopCurrent<TDomain, TAlgebra>::AuxLaplaceLocAss::AuxLaplaceLocAss
(
	NedelecLoopCurrent & master
)
:	IElemDisc<TDomain> (master.m_spVertApproxSpace->name(0).c_str(), master.m_allSsNames.c_str()),
	m_master (master),
	m_posSsGrp (master.m_spVertApproxSpace->subset_handler (), master.m_posSsNames),
	m_cutSsGrp (master.m_spVertApproxSpace->subset_handler (), master.m_cutSsNames),
	m_in_pos_subset (false)
{
//	check the data
	if (m_posSsGrp.empty ())
		UG_THROW ("NedelecLoopCurrent: No positive direction subsets specified");
	if (m_cutSsGrp.empty ())
		UG_THROW ("NedelecLoopCurrent: No cut specified");
//	register assemble functions
	register_all_loc_discr_funcs ();
}

// check the basis and the grid
template<typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::AuxLaplaceLocAss::prepare_setting
(
	const std::vector<LFEID> & vLfeID,
	bool bNonRegular
)
{
	if (bNonRegular)
		UG_THROW ("NedelecLoopCurrent:"
				" The discretization of the auxiliary systems does not support hanging nodes.\n");

	if (vLfeID.size () != 1)
		UG_THROW ("NedelecLoopCurrent:"
			" Only scalar grid functions are supported for the potential");

	if (vLfeID[0].type() != LFEID::LAGRANGE || vLfeID[0].order() !=  1)
		UG_THROW ("NedelecLoopCurrent:"
			" Only Largange-1 functions are supported for the potential");
}

// register the local assembler functions for all the elements and dimensions
template<typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::AuxLaplaceLocAss::register_all_loc_discr_funcs ()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<WDim>::DimElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList> (RegisterLocalDiscr (this));
}

// register the local assembler functions for a given element
template<typename TDomain, typename TAlgebra>
template<typename TElem> // the element to register for
void NedelecLoopCurrent<TDomain, TAlgebra>::AuxLaplaceLocAss::register_loc_discr_func ()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;

	this->enable_fast_add_elem (true);
	this->set_prep_elem_loop_fct(id, & AuxLaplaceLocAss::template prepare_element_loop<TElem>);
	this->set_prep_elem_fct		(id, & AuxLaplaceLocAss::template prepare_element<TElem>);
	this->set_fsh_elem_loop_fct	(id, & AuxLaplaceLocAss::template finish_element_loop<TElem>);
	this->set_add_jac_A_elem_fct(id, & AuxLaplaceLocAss::template ass_JA_elem<TElem>);
	this->set_add_jac_M_elem_fct(id, & AuxLaplaceLocAss::template ass_JM_elem<TElem>);
	this->set_add_def_A_elem_fct(id, & AuxLaplaceLocAss::template ass_dA_elem<TElem>);
	this->set_add_def_M_elem_fct(id, & AuxLaplaceLocAss::template ass_dM_elem<TElem>);
	this->set_add_rhs_elem_fct	(id, & AuxLaplaceLocAss::template ass_rhs_elem<TElem>);
}

// prepares the loop over the elements: check whether we are in the positive subset
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::AuxLaplaceLocAss::prepare_element_loop
(
	ReferenceObjectID roid, // [in] only elements with this roid are looped over
	int si // [in] and only in this subdomain
)
{
	m_in_pos_subset = m_posSsGrp.contains (si);
}

/// transfer the local stiffness matrix to the global discretization
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::AuxLaplaceLocAss::ass_JA_elem
(
	LocalMatrix & J, ///< [in] the matrix to update
	const LocalVector & u, ///< [in] the local solution (not used here)
	GeometricObject * elem, ///< [in] element to prepare
	const position_type vCornerCoords [] ///< [in] coordinates of the corners of the element
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	assemble the local matrix	
	number loc_A [ref_elem_type::numCorners] [ref_elem_type::numCorners];
	NedelecLoopCurrent<TDomain, TAlgebra>::template LocLaplaceA<TElem>::stiffness (elem, vCornerCoords, loc_A);
	
//	add the local matrix to the global one
	for (size_t from_co = 0; from_co < (size_t) ref_elem_type::numCorners; from_co++)
		for (size_t to_co = 0; to_co < (size_t) ref_elem_type::numCorners; to_co++)
			J (_C_, from_co, _C_, to_co) += loc_A [from_co] [to_co];
}

// computes the local right-hand side
template<typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::AuxLaplaceLocAss::ass_rhs_elem
(
	LocalVector & d, ///< [in] the right-hand side to assemble
	GeometricObject * elem, ///< [in] element to prepare
	const position_type vCornerCoords [] ///< [in] coordinates of the corners of the element
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	
	if (! m_in_pos_subset) return; // rhs is nonzero only in one subset

//	assemble the local matrix
	bool at_cut [ref_elem_type::numCorners];
	if (! NedelecLoopCurrent<TDomain, TAlgebra>::template LocLaplaceA<TElem>::bnd_corners
		(elem, m_master.m_cutSsGrp, at_cut))
		return; // no corners on the cut

//	assemble the local matrix and compose the right-hand side of it
	number loc_A [ref_elem_type::numCorners] [ref_elem_type::numCorners];
	NedelecLoopCurrent<TDomain, TAlgebra>::template LocLaplaceA<TElem>::stiffness
		(elem, vCornerCoords, loc_A);
	for (size_t i = 0; i < (size_t) ref_elem_type::numCorners; i++)
		for (size_t j = 0; j < (size_t) ref_elem_type::numCorners; j++)
			if (at_cut [j])
				d (_C_, i) -= loc_A [i] [j]; // to simulate the jump of the potential
}

/*----- Assembling the auxiliary Poisson problems 2: class OutOfSource -----*/

/**
 * Class constructor:
 */
template <typename TDomain, typename TAlgebra>
NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::OutOfSource
(
	NedelecLoopCurrent & master
)
:	m_master (master)
{}

/**
 * Marks the vertices belonging to elements in the source (for one element type)
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::mark_source_vertices_elem_type
(
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	VariableArray1<bool> & in_source ///< [out] the arrays of flags to update
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
	std::vector<size_t> vVertInd (ref_elem_type::numCorners);
	
//	Loop over the source subsets:
	for (size_t i = 0; i < m_master.m_allSsGrp.size (); i++)
	{
		int si = m_master.m_allSsGrp [i];
	//	Loop over all the elements of the given type in the subset
		iterator e_end = vertDD.template end<TElem> (si);
		for (iterator elem_iter = vertDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
			TElem * pElem = *elem_iter;
			if (vertDD.algebra_indices (pElem, vVertInd) != (size_t) ref_elem_type::numCorners)
				UG_THROW ("NedelecLoopCurrent: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			for (size_t i = 0; i < (size_t) ref_elem_type::numCorners; i++)
				in_source [vVertInd [i]] = true;
		}
	}
}

/**
 * Marks the vertices belonging to elements in the source (for all element types)
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::mark_source_vertices
(
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	VariableArray1<bool> & in_source ///< [out] the arrays of flags to update
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	
	UG_ASSERT (vertDD.num_indices () == in_source.size (),
		"NedelecLoopCurrent: Size mismatch of the vertex DoF distribution.");
	memset (& (in_source.at (0)), 0, in_source.size () * sizeof (bool));
	boost::mpl::for_each<ElemList> (MarkSourceVertices (this, vertDD, in_source));
}

/**
 * Sets to identity all the matrix rows that do not belong to the closure
 * of the source subdomain
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::adjust_matrix
(
	VariableArray1<bool> & in_source, ///< the arrays of flags
	pot_matrix_type & A ///< the matrix to adjust
)
{
	UG_ASSERT (A.num_rows () == in_source.size (), "Matrix size mismatch.")
	for (size_t i = 0; i < in_source.size (); i++)
	if (! in_source[i])
		SetDirichletRow (A, i);
}

/**
 * Sets to 0 all the entries of a vector that do not belong to the closure
 * of the source subdomain
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::adjust_vector
(
	VariableArray1<bool> & in_source, ///< the arrays of flags
	pot_vector_type & u ///< the vector to adjust
)
{
	UG_ASSERT (u.size () == in_source.size (), "Vector size mismatch.")
	for (size_t i = 0; i < in_source.size (); i++)
	if (! in_source[i])
		u [i] = 0;
}

/*----- Computation of the flux of the potential over the cut -----*/

/**
 * Computes the flux of the gradient of the potential over the cut
 * (for one type of elements)
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::get_flux_of_pot
(
	const domain_type & domain, ///< [in] the domain
	const pot_vector_type & pot, ///< [in] the potential field
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	number & flux ///< [out] the flux to update
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
//	Get the positions of the grid points:
	const typename TDomain::position_accessor_type & aaPos = domain.position_accessor ();

//	Get the 'negative' subset indices
	SubsetGroup negSsGrp (m_allSsGrp);
	negSsGrp.remove (m_posSsGrp);
	
//	Loop the elements in the subset group:
	for (size_t i = 0; i < negSsGrp.size (); i++)
	{
		int si = negSsGrp [i];
		iterator e_end = vertDD.template end<TElem> (si);
		for (iterator elem_iter = vertDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
			TElem * pElem = *elem_iter;
			bool bnd_flag [ref_elem_type::numCorners];
			
		//	Check whether we are at the cut
			if (! LocLaplaceA<TElem>::bnd_corners (pElem, m_cutSsGrp, bnd_flag))
				continue; // this is not the case
			
		//	Get the corner positions:
			position_type aCorners [ref_elem_type::numCorners];
			for (size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
				aCorners [co] = aaPos [pElem->vertex (co)];
			
		//	Assemble the local Laplacian:
			number loc_A [ref_elem_type::numCorners] [ref_elem_type::numCorners];
			LocLaplaceA<TElem>::stiffness (pElem, aCorners, loc_A);
		
		//	Compute the local contributions to the flux:
			std::vector<size_t> vVertInd (1);
			for (size_t i = 0; i < (size_t) ref_elem_type::numCorners; i++)
			if (bnd_flag [i]) // consider only the corners at the cut
				for (size_t j = 0; j < (size_t) ref_elem_type::numCorners; j++)
				{
					if (vertDD.inner_algebra_indices (pElem->vertex (j), vVertInd) != 1)
						UG_THROW ("NedelecLoopCurrent: Illegal vertex-centered DoF distribution");
					flux += loc_A [i] [j] * pot [vVertInd [0]];
				}
		}
	}
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
