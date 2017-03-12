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
 * nedelec_source_impl.h - implementation of class members for computation of
 * divergence-free sources in form of Nedelec-type-1 element DoFs.
 */
#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"

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
	pot_gf_type pot_u (m_spVertApproxSpace, g_lev);
#	ifdef UG_PARALLEL
	pot_u.set_storage_type (PST_CONSISTENT);
#	endif
	compute_potential (pot_u);
	
//	Compute the normalization factor of the potential (to scale the current to 1)
	number pot_scaling;
	boost::mpl::for_each<ElemList>
		(GetFluxOfPotential (this, * m_spVertApproxSpace->domain().get (), pot_u,
			* vertDD.get (), pot_scaling));
#	ifdef UG_PARALLEL
	{
		pcl::ProcessCommunicator proc_comm;
		pot_scaling = proc_comm.allreduce (pot_scaling, PCL_RO_SUM);
	}
#	endif
	pot_scaling = - pot_scaling;
	
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
		number value = m_vSrcData[i_data].I / pot_scaling;
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
	aa_edge_flag_type & in_source ///< [out] the array of flags to update
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
	const ISubsetHandler * pIsh = edgeDD.subset_handler().get ();
	Grid::edge_traits::secure_container edge_list;
	
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
			pIsh->grid()->associated_elements (edge_list, pElem);
			UG_ASSERT ((edge_list.size () == (size_t) ref_elem_type::numEdges),
				"Mismatch of numbers of corners and vertices of an element");
			for (size_t e = 0; e < (size_t) ref_elem_type::numEdges; e++)
			{
				Edge * pEdge = edge_list[e];
				char flag = 1;
				if (in_pos)
				{
					if (m_cutSsGrp.contains (pIsh->get_subset_index (pEdge->vertex (0))))
						flag |= 2;
					if (m_cutSsGrp.contains (pIsh->get_subset_index (pEdge->vertex (1))))
						flag |= 4;
				}
				in_source [pEdge] = flag;
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
	pot_gf_type & pot_u ///< a grid function for the potential
)
{
	pot_gf_type pot_rhs (pot_u.approx_space (), pot_u.dof_distribution ());
	
//	Prepare the attachment for the flags:
	MultiGrid * mg = pot_u.dd()->multi_grid().get ();
	a_vert_flag_type a_in_source;
	mg->attach_to_vertices (a_in_source);
	m_outOfSource->init (pot_u.domain().get (), a_in_source);
	
//	Assemble the matrix of the auxiliary problem:
	pot_u.set (0.0);
	m_auxLaplaceOp->set_level (pot_u.grid_level ());
	m_auxLaplaceOp->init_op_and_rhs (pot_rhs);
	
//	Initizlize the solver:
	m_potSolver->init (m_auxLaplaceOp);
	
//	Assemble and solve the Laplace equation:
	m_potSolver->apply (pot_u, pot_rhs);
	
//	Release the attachment:
	mg->detach_from_vertices (a_in_source);
}

/**
 * Computes the gradient of the potential for one function.
 *
 * Remark: The potential is considered to be properly scaled. Otherwise,
 * the scaling should be included into the value of the source
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::distribute_source_potential
(
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	pot_vector_type & src_pot, ///< [in] potential to distribute
	DoFDistribution & edgeDD, ///< [in] the edge DD
	size_t func, ///< [in] index of the function
	number value, ///< [in] value of the source (probably divided by the scaling of the potential)
	vector_type & src_field ///< [out] the computed source field
)
{
	typedef DoFDistribution::traits<Edge>::const_iterator t_edge_iter;
	
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	
//	Multigrid and iterators:
	SmartPtr<MultiGrid> sp_mg = edgeDD.multi_grid ();
	t_edge_iter edgeIterBeg = edgeDD.begin<Edge> ();
	t_edge_iter edgeIterEnd = edgeDD.end<Edge> ();
	
//	Mark the edges in the source:
	a_edge_flag_type a_in_source;
	sp_mg->attach_to_edges (a_in_source);
	aa_edge_flag_type aa_in_source (*sp_mg, a_in_source);
	SetAttachmentValues (aa_in_source, edgeIterBeg, edgeIterEnd, 0);
	boost::mpl::for_each<ElemList> (MarkSourceEdges (this, edgeDD, aa_in_source));
#	ifdef UG_PARALLEL
	AttachmentAllReduce<Edge> (*sp_mg, a_in_source, PCL_RO_BOR);
#	endif
	
//	Compute the gradient:
	std::vector<size_t> vVertInd (1);
	std::vector<DoFIndex> vEdgeInd (1);
	for (t_edge_iter edgeIter = edgeIterBeg; edgeIter != edgeIterEnd; ++edgeIter)
	{
		Edge * pEdge = *edgeIter;
		char edge_flag;
		number nd_pot [2];
		
	// Get the edge DoF and check whether the edge is in the source:
		if (edgeDD.inner_dof_indices (pEdge, func, vEdgeInd) != 1)
			UG_THROW ("NedelecLoopCurrent: Edge DoF distribution mismatch. Not the Nedelec-Type-1 element?");
		if ((edge_flag = aa_in_source [pEdge]) == 0)
		{ // we are not in the source
			DoFRef (src_field, vEdgeInd [0]) = 0;
			continue;
		}
		
	// Get the values at the ends:
		for (size_t i = 0; i < 2; i++)
		{
			Vertex * pVertex = pEdge->vertex (i);
			if (vertDD.inner_algebra_indices (pVertex, vVertInd) != 1)
				UG_THROW ("NedelecLoopCurrent: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			nd_pot [i] = src_pot [vVertInd [0]];
			if ((edge_flag & (2 << i)) != 0)
				nd_pot [i] += 1.0; // the jump of the potential at the cut
		}
		
	// Compute the gradient:
		DoFRef (src_field, vEdgeInd [0]) = value * (nd_pot [1] - nd_pot [0]);
	}
	
//	Release the attachment:
	sp_mg->detach_from_edges (a_in_source);
}

/*----- Assembling the auxiliary Poisson problems 1: class AuxLaplaceLocAss -----*/

/**
 * Computes the local discretization of the Laplace operator
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecLoopCurrent<TDomain, TAlgebra>::LocLaplaceA<TElem>::stiffness
(
	GridObject * elem, ///< [in] element to prepare
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
	GridObject * elem, ///< [in] element to compute for
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
	
	this->clear_add_fct(id);
	
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
	GridObject * elem, ///< [in] element to prepare
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
	GridObject * elem, ///< [in] element to prepare
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
	const TDomain * dom ///< [in] the domain
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename geometry_traits<TElem>::const_iterator iterator;
	
//	Get the multigrid and the subset handler:
	const MultiGrid * mg = dom->grid().get ();
	const MGSubsetHandler * ssh = dom->subset_handler().get ();
	size_t n_levels = mg->num_levels ();
	
//	Loop over the source subsets:
	for (size_t i = 0; i < m_master.m_allSsGrp.size (); i++)
	{
		int si = m_master.m_allSsGrp [i];
		for (size_t lev = 0; lev < n_levels; lev++)
		{
		//	Loop over all the elements of the given type in the subset
			iterator e_end = ssh->template end<TElem> (si, lev);
			for (iterator elem_iter = ssh->template begin<TElem> (si, lev);
					elem_iter != e_end; ++elem_iter)
			{
				TElem * pElem = *elem_iter;
				for (size_t i = 0; i < (size_t) ref_elem_type::numCorners; i++)
					m_in_source [pElem->vertex (i)] = true;
			}
		}
	}
}

/**
 * Marks the vertices belonging to elements in the source (for all element types)
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::mark_source_vertices
(
	const TDomain * dom ///< [in] the domain
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	
//	Reset the flags:
	const MultiGrid * mg = dom->grid().get ();
	SetAttachmentValues (m_in_source, mg->begin<Vertex> (), mg->end<Vertex> (), false);
		
//	Mark the vertices:
	boost::mpl::for_each<ElemList> (MarkSourceVertices (this, dom));
}

/**
 * Sets to identity all the matrix rows that do not belong to the closure
 * of the source subdomain
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::adjust_matrix
(
	const DoFDistribution & vertDD, ///< the vertex DD
	pot_matrix_type & A ///< the matrix to adjust
)
{
	typedef DoFDistribution::traits<Vertex>::const_iterator iterator;
	std::vector<size_t> vVertInd (1);
	
//	Loop over all the vertices out of the source
	iterator vert_end = vertDD.end<Vertex> ();
	for (iterator vert_iter = vertDD.begin<Vertex> (); vert_iter != vert_end; ++vert_iter)
	{
		Vertex * pVertex = *vert_iter;
		if (m_in_source [pVertex])
			continue; // the vertex is in the source
		
		if (vertDD.inner_algebra_indices (pVertex, vVertInd) != 1)
			UG_THROW ("NedelecLoopCurrent: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
		
		SetDirichletRow (A, vVertInd[0]);
	}
}

/**
 * Sets to 0 all the entries of a vector that do not belong to the closure
 * of the source subdomain
 */
template <typename TDomain, typename TAlgebra>
void NedelecLoopCurrent<TDomain, TAlgebra>::OutOfSource::adjust_vector
(
	const DoFDistribution & vertDD, ///< the vertex DD
	pot_vector_type & u ///< the vector to adjust
)
{
	typedef DoFDistribution::traits<Vertex>::const_iterator iterator;
	std::vector<size_t> vVertInd (1);
	
//	Loop over all the vertices out of the source
	iterator vert_end = vertDD.end<Vertex> ();
	for (iterator vert_iter = vertDD.begin<Vertex> (); vert_iter != vert_end; ++vert_iter)
	{
		Vertex * pVertex = *vert_iter;
		if (m_in_source [pVertex])
			continue; // the vertex is in the source
		
		if (vertDD.inner_algebra_indices (pVertex, vVertInd) != 1)
			UG_THROW ("NedelecLoopCurrent: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
		
		u [vVertInd[0]] = 0;
	}
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
