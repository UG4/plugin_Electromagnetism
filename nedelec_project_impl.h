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
 * nedelec_project_impl.h - implementation of class members of the classes for
 * the projection of the functions based on the Nedelec element to the space of
 * the divergence-free functions.
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
NedelecProject<TDomain, TAlgebra>::NedelecProject
(
	SmartPtr<EMaterial<TDomain> > emMatherial, ///< properties of the materials
	SmartPtr<ApproximationSpace<TDomain> > vertApproxSpace, ///< vertex-centered approx. space
	SmartPtr<ILinearOperatorInverse<pot_vector_type> > potSolver ///< linear solver for the potential
)
:	m_spEmMaterial (emMatherial),
	m_spVertApproxSpace (vertApproxSpace),
	m_bDampDVFs (true),
	m_auxLocLaplace (new AuxLaplaceLocAss (*this)),
	m_auxLaplaceRHS (new AuxLaplaceRHS (*this)),
	m_auxLaplaceAss (new DomainDiscretization<TDomain, TPotAlgebra> (vertApproxSpace)),
	m_auxLaplaceOp (new AssembledLinearOperator<TPotAlgebra> (SmartPtr<IAssemble<TPotAlgebra> >(m_auxLaplaceAss))),
	m_potSolver (potSolver)
{
//	Check the parameters:
	if (m_spEmMaterial.invalid ())
		UG_THROW ("NedelecProject: Object of the matherial properties not specified.");
	if (m_spVertApproxSpace.invalid ())
		UG_THROW ("NedelecProject: Illegal vert.-centered approx. space.");
	if (m_spVertApproxSpace->num_fct () != 1)
		UG_THROW ("NedelecProject: Exactly one function should be defined in the vert.-centered approx. space.");
	if (! m_spVertApproxSpace->is_def_everywhere (0))
		UG_THROW ("NedelecProject: The function in the vert.-centered approx. space must be defined everywhere.");
	if (m_potSolver.invalid ())
		UG_THROW ("NedelecProject: Illegal solver for the auxiliary problems.");
	
//	Compose the global discretization of the auxiliary equations:
	m_auxLaplaceAss->add (SmartPtr<IElemDisc<TDomain> >(m_auxLocLaplace));
	m_auxLaplaceAss->add
		(SmartPtr<IDomainConstraint<TDomain, TPotAlgebra> >(m_auxLaplaceRHS));
}

/**
 * Performs the projection
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::apply
(
	SmartPtr<GridFunction<TDomain, TAlgebra> > sp_u, ///< the grid function to project
	const char * fct_names ///< the function name
)
{
//	Check the data:
	if (sp_u.invalid ())
		UG_THROW ("NedelecProject: Illegal grid function specification.");
	GridFunction<TDomain, TAlgebra> & u = * sp_u.get ();
	
	if (! m_spEmMaterial->finalized ())
		UG_THROW ("NedelecProject: The material data structure has not been finalized.");
	
//	Get the domain:
	SmartPtr<domain_type> domain = u.domain ();
	if (domain.get () != m_spVertApproxSpace->domain().get ())
		UG_THROW ("NedelecProject: The approximation spaces are based on different domains.");
	
//	Get the function indices:
	FunctionGroup fctGrp;
	try
	{
		fctGrp = u.fct_grp_by_name (fct_names);
	}
	UG_CATCH_THROW ("NedelecProject: Functions '" << fct_names << "' not all contained in the edge approximation space.");
	
	for (size_t i_fct = 0; i_fct < fctGrp.size (); i_fct++)
		if (u.local_finite_element_id(fctGrp[i_fct]).type () != LFEID::NEDELEC)
			UG_THROW ("NedelecProject: Not a Nedelec-element-based grid function specified for the projection.");
	
//	Get the DoF distributions:
	SmartPtr<DoFDistribution> edgeDD = u.dof_distribution ();
	SmartPtr<DoFDistribution> vertDD = m_spVertApproxSpace->dof_distribution (u.grid_level ());
	
//	Create temporary grid functions for the auxiliary problem
	pot_gf_type aux_rhs (m_spVertApproxSpace, vertDD);
	pot_gf_type aux_cor (m_spVertApproxSpace, vertDD);
	
//	Assemble the matrix of the auxiliary problem:
	aux_cor.set (0.0);
	m_auxLaplaceOp->set_level (u.grid_level ());
	m_auxLaplaceOp->init (aux_cor);

//	Initizlize the solver:
	m_potSolver->init (m_auxLaplaceOp);

//	Compute the Dirichlet vector fields:
	if (m_bDampDVFs)
	{
		alloc_DVFs (domain, aux_rhs);
		compute_DVFs (aux_rhs);
		compute_DVF_potential_coeffs (domain, vertDD);
	}
	
//	Project every function:
	for (size_t i_fct = 0; i_fct < fctGrp.size (); i_fct++)
		project_func (domain, edgeDD, u, fctGrp[i_fct], vertDD, aux_rhs, aux_cor);
	
//	Release the Dirichlet vector fields:
	if (m_bDampDVFs)
	{
		for (size_t i = 0; i < m_DVF_phi.size (); i++)
			delete m_DVF_phi [i];
		m_DVF_phi.resize (0);
	}
}

/**
 * Computes weak divergence in insulators and saves it in a vertex-centered grid function
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::compute_div
(
	SmartPtr<GridFunction<TDomain, TAlgebra> > sp_u, ///< [in] the vector field grid function
	const char * u_fct_name, ///< [in] the function name of the Nedelec DoFs
	SmartPtr<GridFunction<TDomain, TPotAlgebra> > sp_div ///< [out] the grid function for the divergence
)
{
//	Check the data:
	if (sp_u.invalid ())
		UG_THROW ("NedelecProject: Illegal input grid function specification.");
	GridFunction<TDomain, TAlgebra> & u = * sp_u.get ();
	
	if (sp_div.invalid ())
		UG_THROW ("NedelecProject: Illegal output grid function specification.");
	pot_vector_type & div = * sp_div.get ();
	
	if (! m_spEmMaterial->finalized ())
		UG_THROW ("The material data structure has not been finalized.");
	
//	Get the domain:
	SmartPtr<domain_type> domain = u.domain ();
	if (domain.get () != m_spVertApproxSpace->domain().get ())
		UG_THROW ("NedelecProject: The approximation spaces are based on different domains.");
	
//	Get the function index of the vector field:
	size_t u_fct;
	try
	{
		u_fct = u.fct_id_by_name (u_fct_name);
	}
	UG_CATCH_THROW (" Function '" << u_fct_name << "' not contained in the edge approximation space.");
	
	if (u.local_finite_element_id(u_fct).type () != LFEID::NEDELEC)
		UG_THROW ("NedelecProject: Not a Nedelec-element-based input grid function.");
	
//	Get the DoF distributions:
	SmartPtr<DoFDistribution> edgeDD = u.dof_distribution ();
	SmartPtr<DoFDistribution> vertDD = sp_div->dof_distribution ();
	if (vertDD.get () != m_spVertApproxSpace->dof_distribution(u.grid_level ()).get ())
		UG_THROW ("NedelecProject: DoFDistribution mismatch, illegal output grid function.");
	
//	Compute the weak divergence:
	DenseVector<VariableArray1<number> > charge;
	div.set (0.0);
	assemble_div (* domain.get(), * edgeDD.get(), u, u_fct, * vertDD.get(), div, charge);
}

/**
 * Allocates memory for the DVFs associated with the conductors that
 * do not touch the Dirichlet boundary
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::alloc_DVFs
(
	const SmartPtr<TDomain> & domain, ///< [in] the domain
	pot_gf_type & aux_rhs ///< [in] the grid function of the auxiliary rhs
)
{
	const DoFDistribution * vertDD = aux_rhs.dof_distribution().get ();
	const std::vector<int> & cond_index = m_spEmMaterial->base_conductor_index ();
	const std::vector<int> & base_cond = m_spEmMaterial->base_conductors ();
	size_t n_cond = base_cond.size ();
	std::vector<bool> grounded (n_cond);
	
	for (size_t i = 0; i < n_cond; i++)
		grounded [i] = false;
	
//	Exclude grounded conductors
	if (m_spDirichlet.valid ())
	{
		typedef DoFDistribution::traits<Edge>::const_iterator t_edge_iterator;
		
		const grid_type * grid = domain->grid().get ();
		const subset_handler_type * ss_handler = domain->subset_handler().get ();
		Grid::volume_traits::secure_container volume_list;
		
		SubsetGroup dirichlet_ssgrp (domain->subset_handler());
		m_spDirichlet->get_dirichlet_subsets (dirichlet_ssgrp);
		
	//	Loop the Dirichlet subsets
		for (size_t j = 0; j < dirichlet_ssgrp.size (); j++)
		{
			int si = dirichlet_ssgrp [j];
			
		//	Loop the Dirichlet edges in the subset
			t_edge_iterator iterEnd = vertDD->end<Edge> (si);
			for (t_edge_iterator iter = vertDD->begin<Edge> (si); iter != iterEnd; iter++)
			{
			//	Loop the adjacent volumes
				((grid_type*) grid)->associated_elements (volume_list, *iter);
				for (size_t v = 0; v < volume_list.size (); v++)
				{
					int v_b_cnd;
					if ((v_b_cnd = cond_index [ss_handler->get_subset_index (volume_list [v])]) < 0)
						continue;
					grounded [v_b_cnd] = true;
				}
			}
		}
	}
	
//	Allocate memory for the conductors
	m_DVF_phi.resize (n_cond);
	for (size_t i = 0; i < n_cond; i++)
	if (grounded [i])
		m_DVF_phi [i] = NULL; // the conductor is grounded, skip it
	else
		m_DVF_phi [i] = new pot_gf_type (aux_rhs.approx_space (), aux_rhs.dof_distribution ());
}

/**
 * Computes the Dirichlet vector fields
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::compute_DVFs
(
	pot_gf_type & aux_rhs //< grid function for the auxiliary rhs
)
{
	for (size_t i = 0; i < m_DVF_phi.size (); i++)
	{
		pot_gf_type * phi = m_DVF_phi [i];
		if (phi == NULL) continue; // a grounded conductor
		
	// 1. Compose the right-hand side:
		m_auxLaplaceRHS->set_base_conductor (i);
		aux_rhs.set (0.0);
		m_auxLaplaceAss->adjust_solution (aux_rhs, aux_rhs.dof_distribution ());
	// 2. Solve the auxiliary system:
		phi->set (0.0);
		m_auxLaplaceAss->adjust_solution (*phi, phi->dof_distribution ());
		m_potSolver->apply (*phi, aux_rhs);
	}
}

/**
 * Computes the potential coefficients
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::compute_DVF_potential_coeffs
(
	const SmartPtr<TDomain> & domain, ///< [in] the domain
	const SmartPtr<DoFDistribution> & vertDD ///< [in] the vertex DD
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	
//	Prepare the matrix for the potential coefficients:
	size_t num_b_cond = m_DVF_phi.size ();
	m_potCoeff.resize (num_b_cond, num_b_cond, false);
	if (num_b_cond <= 0) return; // nothing to do: no conductors
	
//	Get the base conductor indices at vertices
	size_t num_vert = vertDD->num_indices ();
	std::vector<int> vert_base_cond (num_vert, -1);
	boost::mpl::for_each<ElemList> (MarkCondVert (this, * vertDD.get (), & (vert_base_cond[0])));
	for (size_t i = 0; i < num_vert; i++)
	{
		int & b_cond = vert_base_cond [i];
		if (b_cond >= 0 && m_DVF_phi [b_cond] == NULL)
			b_cond = -2; // exclude grounded conductors
	}
	
//	Compute the capacity matrix
	m_potCoeff = 0.0;
	boost::mpl::for_each<ElemList> (IntegrateDivDVF
		(this, * domain.get (), * vertDD.get (), & (vert_base_cond[0]), m_potCoeff));
	for (size_t i = 0; i < num_b_cond; i++)
	if (m_DVF_phi [i] == NULL)
		m_potCoeff (i, i) = 1; // set the matrix to identity for the grounded conductors
	else
		for (size_t j = i + 1; j < num_b_cond; j++)
			m_potCoeff (i, j) = m_potCoeff (j, i); // use the symmetry in the lower triangular part
	
//	Invert the capacity matrix to get the potential coefficients
	Invert (m_potCoeff);
}

/**
 * Projects one function, thus performs the main action of the projection.
 * This function assumes that the object is completely initialized, i.e.
 * the solver and the discretization are initialized.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::project_func
(
	const SmartPtr<TDomain> & domain, ///< [in] the domain
	const SmartPtr<DoFDistribution> & edgeDD, ///< [in] the edge DD
	vector_type & u, ///< [in] the vector where to project
	size_t fct, ///< [in] function in u to compute div for
	const SmartPtr<DoFDistribution> & vertDD, ///< [in] the vertex DD
	pot_gf_type & aux_rhs, ///< [in] a grid function for the rhs of the aux. problem
	pot_gf_type & aux_cor ///< [in] a grid function for the sol. of the aux. problem
)
{
//--- Compute the correction due to the divergence:

	aux_cor.set (0.0);
	aux_rhs.set (0.0);
	
	DenseVector<VariableArray1<number> > charge;
	
// 1. Compute the weak div:
	assemble_div (* domain.get(), * edgeDD.get(), u, fct, * vertDD.get(), aux_rhs, charge);
// 2. Correct the divergence in the right-hand side:
	m_auxLaplaceRHS->set_base_conductor (-1);
	m_auxLaplaceAss->adjust_solution (aux_rhs, vertDD);
// 3. Solve the auxiliary system:
	m_auxLaplaceAss->adjust_solution (aux_cor, vertDD);
	m_potSolver->apply (aux_cor, aux_rhs);
// 4. Damp the Dirichlet vector fields:
	if (m_bDampDVFs)
		damp_DVFs (aux_cor, charge);
// 5. Merge the correction into the solution:
	distribute_cor (* domain.get(), * edgeDD.get(), u, fct, * vertDD.get(), aux_cor);
}

/**
 * Assembles the weak divergence operator for all elements of the same type.
 * Remark: Only full-dimensional elements should be considered.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecProject<TDomain, TAlgebra>::weak_div_elem_type
(
	const TDomain & domain, ///< [in] the domain
	const DoFDistribution & edgeDD, ///< [in] the edge DD
	const vector_type & u, ///< [in] the vector to compute div for
	size_t fct, /// [in] function in u to compute div for
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	pot_vector_type & div ///< to update the weak divergence of u
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
//	Get the conductor distribution and the positions of the grid points:
	const std::vector<int> & cond_index = m_spEmMaterial->base_conductor_index ();
	const typename TDomain::position_accessor_type & aaPos = domain.position_accessor ();

//	Arrays for the indices in the vectors:
	std::vector<DoFIndex> vEdgeInd (ref_elem_type::numEdges);
	std::vector<size_t> vVertInd (ref_elem_type::numCorners);
	
//	Loop over all subsets representing insulators:
	for (int si = 0; si < vertDD.num_subsets (); si++)
	if (cond_index [si] == -1)
	{
	//	Loop over all the elements of the given type in the subset
		iterator e_end = vertDD.template end<TElem> (si);
		for (iterator elem_iter = vertDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
			TElem * pElem = *elem_iter;
			MathMatrix<ref_elem_type::numCorners, ref_elem_type::numEdges> B;
			MathVector<ref_elem_type::numEdges> loc_u;
			MathVector<ref_elem_type::numCorners> loc_div;
			
		//	Compute the local weak divergence matrix:
			position_type aCorners [ref_elem_type::numCorners];
			for (size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
				aCorners [co] = aaPos [pElem->vertex (co)];
			NedelecT1_LDisc<TDomain, TElem>::local_div_matrix (&domain, pElem, aCorners,
				(number (*) [ref_elem_type::numEdges]) & (B (0,0)));
			
		//	Compute the local contribution to the weak divergence:
			if (edgeDD.dof_indices (pElem, fct, vEdgeInd) != (size_t) ref_elem_type::numEdges)
				UG_THROW ("NedelecProject: Edge DoF distribution mismatch. Not the Nedelec-Type-1 element?");
			for (size_t i = 0; i < (size_t) ref_elem_type::numEdges; i++)
				loc_u [i] = DoFRef (u, vEdgeInd [i]);
			MatVecMult (loc_div, B, loc_u);
			
		//	Add the local contribution to the global vector:
			if (vertDD.algebra_indices (pElem, vVertInd) != (size_t) ref_elem_type::numCorners)
				UG_THROW ("NedelecProject: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			for (size_t i = 0; i < (size_t) ref_elem_type::numCorners; i++)
				div [vVertInd [i]] += loc_div [i];
		}
	}
}

/**
 * Sets the div operator to 0 in conductors for all elements of the same type
 * and compute the charges of the base conductors. The charges of the base
 * conductors are the sums of the divergences in all the interface vertices
 * between the base conductors and the insulators. As the divergence inside
 * of the conductors is set to zero, we sum over all the nodes in the closure
 * of the conductors.
 * Remark: Only full-dimensional elements should be considered.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecProject<TDomain, TAlgebra>::clear_div_in_conductors
(
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	pot_vector_type & div, ///< [out] to update the weak divergence of u
	DenseVector<VariableArray1<number> > & charge ///< [out] charges of the conductors
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
	const std::vector<int> & base_cond = m_spEmMaterial->base_conductors ();
	if (base_cond.size () == 0) return; // no conductors
	
	std::vector<size_t> vVertInd (ref_elem_type::numCorners);
	const std::vector<int> & base_cond_index = m_spEmMaterial->base_conductor_index ();
	int ci;

//	Loop over all subsets representing conductors:
	for (int si = 0; si < vertDD.num_subsets (); si++)
	if ((ci = base_cond_index [si]) >= 0)
	{
		number & charge_entry = charge [ci];
		
	//	Loop over all the elements of the given type in the subset
		iterator e_end = vertDD.template end<TElem> (si);
		for (iterator elem_iter = vertDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
			TElem * pElem = *elem_iter;
			if (vertDD.algebra_indices (pElem, vVertInd) != (size_t) ref_elem_type::numCorners)
				UG_THROW ("NedelecProject: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			for (size_t i = 0; i < (size_t) ref_elem_type::numCorners; i++)
			{
				number & div_entry = div [vVertInd [i]];
				charge_entry += div_entry;
				div_entry = 0;
			}
		}
	}
}

/**
 * Assembles the weak divergence operator on the whole grid and computes the
 * charges of the base conductors.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::assemble_div
(
	const TDomain & domain, ///< [in] the domain
	const DoFDistribution & edgeDD, ///< [in] the edge DD
	const vector_type & u, ///< [in] the vector to compute div for
	size_t fct, ///< [in] function in u to compute div for
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	pot_vector_type & div, ///< [out] for the weak divergence of u
	DenseVector<VariableArray1<number> > & charge ///< [out] charges of the conductors
)
{
//	The full-dim. grid element types for this dimension:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	
	const std::vector<int> & base_cond = m_spEmMaterial->base_conductors ();
	charge.resize (base_cond.size ());
	if (charge.size () != 0)
		charge = 0.0;
	
	div.set (0.0);
//	Compute the divergence for all the types of the elements:
	boost::mpl::for_each<ElemList> (WeakDiv (this, domain, edgeDD, u, fct, vertDD, div));
//	Clear the entries at all the points in the closure of the conductors:
	boost::mpl::for_each<ElemList> (ClearDivInConductors (this, vertDD, div, charge));
}

/**
 * Damps the Dirichlet vector fields (DVFs)
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::damp_DVFs
(
	pot_vector_type & cor, ///< the potential correction to update
	const DenseVector<VariableArray1<number> > & charge ///< [in] charges of the conductors
)
{
	DenseVector<VariableArray1<number> > factor = m_potCoeff * charge;
	
	for (size_t i = 0; i < m_DVF_phi.size (); i++)
	if (m_DVF_phi [i] != NULL) // skip grounded conductors
		VecScaleAdd (cor, 1.0, cor, factor [i], * (pot_vector_type *) m_DVF_phi [i]);
}

/**
 * Computes the edge-centered correction from the vertex-centered (potential)
 * one by applying the gradient operator.
 *
 * Note that the edge-centered correction is considered only in the closure
 * of insulators. There, it is the gradient of the vertex-centered (potential)
 * correction. But in the conductors (as well as on the Dirichlet boundaries),
 * the vertex-centered correction is constant by construction so that its
 * gradient is 0 any way. (Adjacent conductors are concidered as one conductor
 * in the projection, so that the different factors for the charges in different
 * conductors does not violate the constant value of the vertex-centered
 * correction in any separated piece of conductors.) For this, the function
 * computes the gradient over the whole domain.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::distribute_cor
(
	const TDomain & domain, ///< [in] the domain
	const DoFDistribution & edgeDD, ///< [in] the edge DD
	vector_type & u, ///< [out] the vector to correct (to update)
	size_t fct, ///< [in] function in u to update
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	const pot_vector_type & cor ///< [in] the potential correction for u
)
{
//	Iterator over edges
	typedef DoFDistribution::traits<Edge>::const_iterator t_edge_iterator;

//	Arrays for the indices in the vectors:
	std::vector<size_t> vVertInd (1);
	std::vector<DoFIndex> vEdgeInd (1);
	
//	Loop over edges:
	t_edge_iterator iterEnd = edgeDD.end<Edge> ();
	for (t_edge_iterator iter = edgeDD.begin<Edge> (); iter != iterEnd; iter++)
	{
		number corner_val [2];
		Edge * pEdge = *iter;
		
	//	Get the potential the ends of the edge:
		for (size_t i = 0; i < 2; i++)
		{
		//	Get the multiindices
			if (vertDD.inner_algebra_indices (pEdge->vertex(i), vVertInd) != 1)
				UG_THROW ("NedelecProject:"
					"More than one DoF per vertex for the auxiliary systems.");

		//	Set the Dirichlet entry
			corner_val [i] = cor [vVertInd[0]];
		}
		
	//	Compute the gradient:
		if (edgeDD.inner_dof_indices (pEdge, fct, vEdgeInd) != 1)
			UG_THROW ("NedelecProject:"
				"More than one DoF per edge. Not the Nedelec-Type-1 element?");
		DoFRef (u, vEdgeInd[0])
			-= corner_val [1] - corner_val [0];
	}
}

/*----- Assembling the auxiliary Poisson problems 1: class AuxLaplaceLocAss -----*/

/**
 * Computes the local discretization of the Laplace operator
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecProject<TDomain, TAlgebra>::LocLaplaceA<TElem>::stiffness
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
 * Class constructor that gets the master NedelecProject object and extracts
 * the data necessary for the base class:
 */
template <typename TDomain, typename TAlgebra>
NedelecProject<TDomain, TAlgebra>::AuxLaplaceLocAss::AuxLaplaceLocAss
(
	NedelecProject & master
)
:	IElemDisc<TDomain> (master.m_spVertApproxSpace->name(0).c_str(), master.m_spEmMaterial->subset_names()),
	m_master (master), m_do_assemble_here (false)
{
//	register assemble functions
	register_all_loc_discr_funcs ();
}

// check the basis and the grid
template<typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceLocAss::prepare_setting
(
	const std::vector<LFEID> & vLfeID,
	bool bNonRegular
)
{
	if (bNonRegular)
		UG_THROW ("NedelecProject:"
				" The discretization of the auxiliary systems does not support hanging nodes.\n");

	if (vLfeID.size () != 1)
		UG_THROW ("NedelecProject:"
			" Only scalar grid functions are supported for the potential");

	if (vLfeID[0].type() != LFEID::LAGRANGE || vLfeID[0].order() !=  1)
		UG_THROW ("NedelecProject:"
			" Only Largange-1 functions are supported for the potential");
}

// register the local assembler functions for all the elements and dimensions
template<typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceLocAss::register_all_loc_discr_funcs ()
{
//	get all grid element types in this dimension and below
	typedef typename domain_traits<WDim>::DimElemList ElemList;

//	switch assemble functions
	boost::mpl::for_each<ElemList> (RegisterLocalDiscr (this));
}

// register the local assembler functions for a given element
template<typename TDomain, typename TAlgebra>
template<typename TElem> // the element to register for
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceLocAss::register_loc_discr_func ()
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

// prepares the loop over the elements: check if the subdomain an insulator
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceLocAss::prepare_element_loop
(
	ReferenceObjectID roid, // [in] only elements with this roid are looped over
	int si // [in] and only in this subdomain
)
{
//	We assemble in insulators only
	if (m_master.m_spEmMaterial->base_conductor_index (si) == -1)
		m_do_assemble_here = true;
	else
		m_do_assemble_here = false;
}

/// transfer the local stiffness matrix to the global discretization
template<typename TDomain, typename TAlgebra>
template<typename TElem>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceLocAss::ass_JA_elem
(
	LocalMatrix & J, ///< [in] the matrix to update
	const LocalVector & u, ///< [in] the local solution (not used here)
	GridObject * elem, ///< [in] element to prepare
	const position_type vCornerCoords [] ///< [in] coordinates of the corners of the element
)
{
	if (! m_do_assemble_here) return;
	
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	assemble the local matrix	
	number loc_A [ref_elem_type::numCorners] [ref_elem_type::numCorners];
	NedelecProject<TDomain, TAlgebra>::template LocLaplaceA<TElem>::stiffness (elem, vCornerCoords, loc_A);
	
//	add the local matrix to the global one
	for (size_t from_co = 0; from_co < (size_t) ref_elem_type::numCorners; from_co++)
		for (size_t to_co = 0; to_co < (size_t) ref_elem_type::numCorners; to_co++)
			J (_C_, from_co, _C_, to_co) += loc_A [from_co] [to_co];
}

/*----- Assembling the auxiliary Poisson problems 2: class AuxLaplaceRHS -----*/

/**
 * Class constructor:
 */
template <typename TDomain, typename TAlgebra>
NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::AuxLaplaceRHS
(
	NedelecProject & master
)
:	m_master (master), m_base_cond (-2)
{}

/**
 * Sets zero values on the whole Dirichlet boundary.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::set_zero_Dirichlet
(
	pot_vector_type & u, ///< the vector where to set zeros
	const DoFDistribution * dd ///< the vert.-based DoF distribution
)
{
	if (m_master.m_spDirichlet.invalid ())
		return; // no Dirichlet boundaries specified
	
	std::vector<size_t> vVertInd (1);
	
	SubsetGroup dirichlet_ssgrp (m_master.m_spEmMaterial->subset_handler ());
	m_master.m_spDirichlet->get_dirichlet_subsets (dirichlet_ssgrp);
	
//	Loop over the Dirichlet subsets
	for (size_t j = 0; j < dirichlet_ssgrp.size (); j++)
	{
		int si = dirichlet_ssgrp [j];
		
	//	Loop the edges
		t_edge_iterator iterEnd = dd->end<Edge> (si);
		for (t_edge_iterator iter = dd->begin<Edge> (si); iter != iterEnd; iter++)
		{
			Edge * pEdge = *iter;
			
		//	For both the ends of the edge
			for(size_t i = 0; i < 2; i++)
			{
			//	Get the multiindices
				if (dd->inner_algebra_indices (pEdge->vertex(i), vVertInd) != 1)
					UG_THROW ("NedelecProject:"
						"More than one DoF per vertex for the auxiliary systems.");
	
			//	Set the Dirichlet entry
				u [vVertInd [0]] = 0;
			}
		}
	}
}

/**
 * Sets identity matrix on the whole Dirichlet boundary.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::set_identity_Dirichlet
(
	pot_matrix_type & A, ///< the matrix to set
	const DoFDistribution * dd ///< the vert.-based DoF distribution
)
{
	if (m_master.m_spDirichlet.invalid ())
		return; // no Dirichlet boundaries specified
	
	std::vector<size_t> vVertInd (1);
	
	SubsetGroup dirichlet_ssgrp (m_master.m_spEmMaterial->subset_handler ());
	m_master.m_spDirichlet->get_dirichlet_subsets (dirichlet_ssgrp);
	
//	Loop over the Dirichlet subsets
	for (size_t j = 0; j < dirichlet_ssgrp.size (); j++)
	{
		int si = dirichlet_ssgrp [j];
		
	//	Loop the edges
		t_edge_iterator iterEnd = dd->end<Edge> (si);
		for (t_edge_iterator iter = dd->begin<Edge> (si); iter != iterEnd; iter++)
		{
			Edge * pEdge = *iter;
			
		//	For both the ends of the edge
			for(size_t i = 0; i < 2; i++)
			{
			//	Get the multiindices
				if (dd->inner_algebra_indices (pEdge->vertex(i), vVertInd) != 1)
					UG_THROW ("NedelecProject:"
						"More than one DoF per vertex for the auxiliary systems.");
	
			//	Set the Dirichlet row
				SetDirichletRow (A, vVertInd[0]);
			}
		}
	}
}

/**
 * Sets constant value on the closure of a full-dimensional subset.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::set_value_on_subset
(
	int si, ///< the subset
	number val, ///< the value to set
	pot_vector_type & u, ///< the vector where to set
	const DoFDistribution * dd ///< the vert.-based DoF distribution
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator t_elem_iterator;
	
	std::vector<size_t> vVertInd (1);
	
//	Loop the elements
	t_elem_iterator iterEnd = dd->end<TElem> (si);
	for (t_elem_iterator iter = dd->begin<TElem> (si); iter != iterEnd; iter++)
	{
		TElem * pElem = *iter;
		
	//	Loop the corners
		for(size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
		{
		//	Get the multiindices
			if (dd->inner_algebra_indices (pElem->vertex(co), vVertInd) != 1)
				UG_THROW ("NedelecProject:"
					"More than one DoF per vertex for the auxiliary systems.");

		//	Set the value
			u [vVertInd [0]] = val;
		}
	}
}

/**
 * Sets identity matrix on the closure of a full-dimensional subset.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::set_identity_on_subset
(
	int si, ///< the subset
	pot_matrix_type & A, ///< the matrix to set
	const DoFDistribution * dd ///< the vert.-based DoF distribution
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator t_elem_iterator;
	
	std::vector<size_t> vVertInd (1);
	
//	Loop the elements
	t_elem_iterator iterEnd = dd->end<TElem> (si);
	for (t_elem_iterator iter = dd->begin<TElem> (si); iter != iterEnd; iter++)
	{
		TElem * pElem = *iter;
		
	//	Loop the corners
		for(size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
		{
		//	Get the multiindices
			if (dd->inner_algebra_indices (pElem->vertex(co), vVertInd) != 1)
				UG_THROW ("NedelecProject:"
					"More than one DoF per vertex for the auxiliary systems.");

		//	Set the row
			SetDirichletRow (A, vVertInd[0]);
		}
	}
}

/**
 * Sets the conductor-rows of a given Jacobian to the identity.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::adjust_jacobian
(
	pot_matrix_type & J,
	const pot_vector_type & u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time,
	ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol,
	const number s_a0
)
{
//	Set all matrix rows at Dirichlet boundaries to identity:
	set_identity_Dirichlet (J, dd.get());
	
// Set all matrix rows in conductors to identity:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	const std::vector<int> & base_cond = m_master.m_spEmMaterial->base_conductor_index ();
	
	for (size_t si = 0; si < base_cond.size (); si++)
	if (base_cond [si] >= 0)
		boost::mpl::for_each<ElemList> (SetIdentityOnSubset (this, si, J, dd.get ()));
}

/**
 * Sets the conductor-entries of a given defect vector to 0.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::adjust_defect
(
	pot_vector_type & d,
	const pot_vector_type & u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time,
	ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol,
	const std::vector<number> * vScaleMass,
	const std::vector<number> * vScaleStiff
)
{
//	Set all entries at Dirichlet boundaries to zero:
	set_zero_Dirichlet (d, dd.get ());
	
// Set all entries in conductors to zero:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	const std::vector<int> & base_cond = m_master.m_spEmMaterial->base_conductor_index ();
	
	for (size_t si = 0; si < base_cond.size (); si++)
	if (base_cond [si] >= 0)
		boost::mpl::for_each<ElemList> (SetValueOnSubset (this, si, 0, d, dd.get ()));
}

/**
 * Sets the conductor-entries of a given solution vector to the boundary values.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProject<TDomain, TAlgebra>::AuxLaplaceRHS::adjust_solution
(
	pot_vector_type	& u,
	ConstSmartPtr<DoFDistribution> dd,
	int type,
	number time
)
{
//	Set all entries at Dirichlet boundaries to identity:
	set_zero_Dirichlet (u, dd.get ());
	
//	Set all entries in conductors to 0 or 1:
	typedef typename domain_traits<WDim>::DimElemList ElemList;
	const std::vector<int> & base_cond = m_master.m_spEmMaterial->base_conductor_index ();
	
	for (size_t si = 0; si < base_cond.size (); si++)
	{
		int the_base_cond = base_cond [si];
		if (the_base_cond < 0) continue;
		number val = (the_base_cond == m_base_cond)? 1 : 0;
		boost::mpl::for_each<ElemList> (SetValueOnSubset (this, si, val, u, dd.get ()));
	}
}

/**
 * Initializes 'vert_base_cond' with the indices of the base conductors
 * according to the closure of the conductive subsets.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecProject<TDomain, TAlgebra>::mark_cond_vert_elem_type
(
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	int * vert_base_cond ///< [out] indices of the base conductors for every vertex
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
	int base_cond_ind;
	
//	Get the conductor distribution and the positions of the grid points:
	const std::vector<int> & cond_index = m_spEmMaterial->base_conductor_index ();

//	Array for the indices in the vectors:
	std::vector<size_t> vVertInd (ref_elem_type::numCorners);
	
//	Loop over all subsets representing conductors:
	for (int si = 0; si < vertDD.num_subsets (); si++)
	if ((base_cond_ind = cond_index [si]) >= 0)
	{
	//	Loop over all the elements of the given type in the subset
		iterator e_end = vertDD.template end<TElem> (si);
		for (iterator elem_iter = vertDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
		//	Add the local contribution to the global vector:
			if (vertDD.algebra_indices (*elem_iter, vVertInd) != (size_t) ref_elem_type::numCorners)
				UG_THROW ("NedelecProject: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			for (size_t i = 0; i < (size_t) ref_elem_type::numCorners; i++)
				vert_base_cond [vVertInd [i]] = base_cond_ind;
		}
	}
}

/**
 * Integration of div E over boundaries of conductors (for one type of elements)
 * to get the capacity matrix. We compute only the lower triangular part of
 * the capacity matrix: This matrix is symmetric.
 */
template <typename TDomain, typename TAlgebra>
template <typename TElem>
void NedelecProject<TDomain, TAlgebra>::integrate_div_DVF_elem_type
(
	const TDomain & domain, ///< [in] the domain
	const DoFDistribution & vertDD, ///< [in] the vertex DD
	const int * vert_base_cond, ///< [in] indices of the base conductors for every vertex
	DenseMatrix<VariableArray2<number> > & C ///< [out] the capacity matrix to update
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	
//	Get the positions of the grid points:
	const typename TDomain::position_accessor_type & aaPos = domain.position_accessor ();

//	Get the conductor distribution and the positions of the grid points:
	const std::vector<int> & cond_index = m_spEmMaterial->base_conductor_index ();

//	Array for the indices in the vectors:
	std::vector<size_t> vVertInd (ref_elem_type::numCorners);
	
//	Loop over all subsets representing insulators:
	for (int si = 0; si < vertDD.num_subsets (); si++)
	if (cond_index [si] == -1)
	{
	//	Loop over all the elements of the given type in the subset
		iterator e_end = vertDD.template end<TElem> (si);
		for (iterator elem_iter = vertDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
			TElem * pElem = *elem_iter;
			size_t corner_cond [ref_elem_type::numCorners];
			bool cond_bnd_flag;
			
		//	get the indices in the vectors:
			if (vertDD.algebra_indices (pElem, vVertInd) != (size_t) ref_elem_type::numCorners)
				UG_THROW ("NedelecProject: Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			cond_bnd_flag = false;
			for (size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
				if ((corner_cond [co] = vert_base_cond [vVertInd [co]]) >= 0)
					cond_bnd_flag = true;
			if (! cond_bnd_flag) continue; // not at a boundary of a non-grounded conductor
			
		//	get the corner positions:
			position_type aCorners [ref_elem_type::numCorners];
			for (size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
				aCorners [co] = aaPos [pElem->vertex (co)];
			
		//	assemble the local matrix	
			number loc_A [ref_elem_type::numCorners] [ref_elem_type::numCorners];
			LocLaplaceA<TElem>::stiffness (pElem, aCorners, loc_A);
			
		//	add the contributions to the capacity matrix
			for (int to_cond = 0; to_cond < (int) m_DVF_phi.size (); to_cond++)
			{
			//	In this computation, C(i,j) is the charge induced by phi(j)
			//	in conductor i. The conductor i, where the charge is induced,
			//	is called 'the from-conductor', and the conductor j, whose
			//	potential induces the charge, is called 'the to-conductor'.
				pot_vector_type * phi = m_DVF_phi [to_cond];
				if (phi == NULL) continue; // this is a grounded conductor
				
				for (size_t to_co = 0; to_co < (size_t) ref_elem_type::numCorners; to_co++)
				{
					int from_cond;
					number phi_val = (* phi) [vVertInd [to_co]];
					for (size_t from_co = 0; from_co < (size_t) ref_elem_type::numCorners; from_co++)
					// Exclude inner vertices of insulators and grounded conductors, as well as the upper triangle
					if ((from_cond = corner_cond [from_co]) >= to_cond)
						C (from_cond, to_cond) += loc_A [from_co] [to_co] * phi_val;
				}
			}
		}
	}
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
