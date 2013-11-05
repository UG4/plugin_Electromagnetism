/**
 * nedelec_transfer_impl.h - implementation of class members of the transfer
 * operators for the Nelelec (Whitney-1) elements.
 * 
 * Created on: 22.03.2013
 * Author: D. Logashenko
 */
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "nedelec_local_ass.h"

namespace ug{
namespace Electromagnetism{

/**
 * Computes the local coordinates of a vertex according to the assumption
 * of the regular refinement. The local coordinates are computed w.r.t.
 * a given basis element that should be associated with the parent of
 * the vertex whose all corners must be corners of the base element.
 * (Otherwise an exception is thrown.) The vector of the local coordinates
 * of the given vertex is the arithmetical average of the vectors of the
 * local coordinates of its parent w.r.t. the base element.
 *
 * Acknowledgements to S. Reiter.
 */
template <typename TDomain, typename TAlgebra, typename TElem>
void NedelecProlongationMatrixHelper<TDomain, TAlgebra, TElem>::GetRegularLocalCoordinate
(
	const MultiGrid * mg, ///< [in] the grid hierarchy
	VertexBase * v, ///< [in] the vertex
	TElem * base, ///< [in] the base element (for the local coordinates)
	MathVector<TElem::dim> & local ///< [out] to save the local coordinates
)
{
//	Get the parent of the vertex (note that this is typically not the 'base')
	GeometricObject * parent = mg->get_parent (v);
	
//	Get the vertices of the parent
	Grid::vertex_traits::secure_container vrts;
	((MultiGrid *) mg)->associated_elements (vrts, parent);
	size_t n_co = vrts.size ();
	UG_ASSERT (n_co != 0, "GetRegularLocalCoordinate: No associated vertices.")
	
//	Get the reference element for the base
	const DimReferenceElement<TElem::dim> & rRefElem = ReferenceElementProvider::get<TElem::dim>
		(geometry_traits<TElem>::REFERENCE_OBJECT_ID);
	
//	Average the local coordinates of the vertices
	local = 0.0;
	for (size_t co = 0; co < n_co; co++)
	{
		int i = GetVertexIndex (base, vrts[co]);
		if (i < 0)
			UG_THROW ("GetRegularLocalCoordinate: No (implicit) parent-child relation between the vertex and the base element.");
		local += rRefElem.corner (i);
	}
	local /= n_co;
}

/**
 * Assembles the prolongation matrix for the DoFs between the edges
 * for one type of the grid elements. (This instantiation works with edges.)
 * To get the total prolongation matrix, the functions for all the types should
 * be called.
 *
 * Note that we assume here that the father edge is subdivided into two
 * edges with equal length (or is not subdivided at all). At least, the
 * prolongation coefficients are computed for the regular refinement.
 */
template <typename TDomain, typename TAlgebra>
void NedelecProlongationMatrixHelper<TDomain, TAlgebra, Edge>::assemble_prolongation_matrix
(
	const domain_type & domain, ///< [in] the domain
	const DoFDistribution & coarseDD, ///< [in] the coarse grid DD
	const DoFDistribution & fineDD, ///< [in] the fine grid DD
	matrix_type & mat, ///< [out] the interpolation matrix
	std::vector<bool> & vIsRestricted ///< [out] whether a coarse grid DoF has children
)
{
	typedef DoFDistribution::traits<EdgeBase>::const_iterator iterator;
	
// Multiindices to access the components
	std::vector<DoFIndex> c_ind (1), f_ind (1);

// Get the multigrid:
	const MultiGrid * grid = coarseDD.multi_grid().get ();
	
// Get the number of the functions:
	size_t num_fct = coarseDD.num_fct ();
	UG_ASSERT (num_fct == fineDD.num_fct (), "NedelecTransfer: Coarse and find DD mismatch.");
	
// Loop over all subsets:
	for (int si = 0; si < coarseDD.num_subsets (); si++)
	{
	// Loop over all edges in the subset
		iterator e_end = coarseDD.template end<EdgeBase> (si);
		for (iterator edge_iter = coarseDD.template begin<EdgeBase> (si);
			edge_iter != e_end; ++edge_iter)
		{
			number coef;
			EdgeBase * c_edge = * edge_iter;
			const size_t n_children = grid->num_children<EdgeBase, EdgeBase> (c_edge);
			
			if (n_children == 0)
				continue;
			else if (n_children == 1)
			{
				EdgeBase * f_edge = grid->get_child<EdgeBase, EdgeBase> (c_edge,  0);
				
			//	Check the edge orientation:
				GeometricObject * corner_0 = grid->get_parent (f_edge->vertex(0));
				if (corner_0 == (GeometricObject *) (c_edge->vertex (0)))
					coef = 1; // the edges should have the same orientation
				else if (corner_0 == (GeometricObject *) (c_edge->vertex (1)))
					coef = -1; // the edges have the inverse orientation
				else
					UG_THROW ("NedelecTransfer: Cannot find out the edge orientation.");
				
			//	Add the values:
				for (size_t fct = 0; fct < num_fct; fct++)
				{
				//	Get the DoFs:
					if (coarseDD.inner_dof_indices (c_edge, fct, c_ind) != 1
							|| fineDD.inner_dof_indices (f_edge, fct, f_ind) != 1)
						UG_THROW ("NedelecTransfer:"
							"More than one DoF per edge. Not the Nedelec-type-1 element?");
				//	Assemble the + or - identity matrix here:
					DoFRef (mat, f_ind[0], c_ind[0]) += coef;
					vIsRestricted [c_ind[0][0]] = true;
				}
			}
			else if (n_children == 2)
			{
				for (size_t child = 0; child < 2; child++)
				{
					EdgeBase * f_edge = grid->get_child<EdgeBase, EdgeBase> (c_edge,  child);
					
				//	Check the edge orientation:
					GeometricObject * corner_0 = grid->get_parent (f_edge->vertex(0));
					GeometricObject * corner_1 = grid->get_parent (f_edge->vertex(1));
					if (corner_0 == (GeometricObject *) (c_edge->vertex (0))
						|| corner_1 == (GeometricObject *) (c_edge->vertex (1)))
						coef = 0.5; // the edges should have the same orientation
					else if (corner_0 == (GeometricObject *) (c_edge->vertex (1))
						|| corner_1 == (GeometricObject *) (c_edge->vertex (0)))
						coef = -0.5; // the edges have the inverse orientation
					else
						UG_THROW ("NedelecTransfer: Cannot find out the edge orientation.");
					
				//	Add the values:
					for (size_t fct = 0; fct < num_fct; fct++)
					{
					//	Get the DoFs:
						if (coarseDD.inner_dof_indices (c_edge, fct, c_ind) != 1
								|| fineDD.inner_dof_indices (f_edge, fct, f_ind) != 1)
							UG_THROW ("NedelecTransfer:"
								"More than one DoF per edge. Not the Nedelec-type-1 element?");
					//	Assemble the matrix entry:
						DoFRef (mat, f_ind[0], c_ind[0]) += coef;
						vIsRestricted [c_ind[0][0]] = true;
					}
				}
			}
			else
				UG_THROW ("NedelecTransfer: Only regular refinement is supported,"
							" but some edge is subdivided into " << n_children << " ones.");
		}
	}
}

/**
 * Assembles the prolongation matrix for the DoFs between the edges
 * for one type of the grid elements. (This instantiation assumes that the
 * dimensionality of the element is greater than 1.) To get the total
 * prolongation matrix, the functions for all the types should be called.
 */
template <typename TDomain, typename TAlgebra, typename TElem>
void NedelecProlongationMatrixHelper<TDomain, TAlgebra, TElem>::assemble_prolongation_matrix
(
	const domain_type & domain, ///< [in] the domain
	const DoFDistribution & coarseDD, ///< [in] the coarse grid DD
	const DoFDistribution & fineDD, ///< [in] the fine grid DD
	matrix_type & mat, ///< [out] the interpolation matrix
	std::vector<bool> & vIsRestricted ///< [out] whether a coarse grid DoF has children
)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator iterator;
	const ReferenceObjectID roid = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	
// Multiindices to access the components
	std::vector<DoFIndex> c_ind (ref_elem_type::numEdges), f_ind (1);
	
// Get the multigrid:
	const MultiGrid * grid = coarseDD.multi_grid().get ();
	
// Get position accessor for the integration:
	const typename TDomain::position_accessor_type & aaPos = domain.position_accessor ();

// Get the number of the functions:
	size_t num_fct = coarseDD.num_fct ();
	UG_ASSERT (num_fct == fineDD.num_fct (), "NedelecTransfer: Coarse and find DD mismatch.");
	
// Loop over all subsets:
	for (int si = 0; si < coarseDD.num_subsets (); si++)
	{
	// Loop over all the elements of the given type in the subset
		iterator e_end = coarseDD.template end<TElem> (si);
		for (iterator elem_iter = coarseDD.template begin<TElem> (si);
			elem_iter != e_end; ++elem_iter)
		{
			TElem * c_elem = * elem_iter;
			const size_t n_children = grid->num_children<EdgeBase, TElem> (c_elem);
			
			if (n_children == 0)
				continue;
			
		// Get the corner positions:
			position_type aCorners [ref_elem_type::numCorners];
			for (size_t co = 0; co < (size_t) ref_elem_type::numCorners; co++)
				aCorners [co] = aaPos [c_elem->vertex (co)];
			
		// Loop over the child elements:
			for (size_t child = 0; child < n_children; child++)
			{
			// Get the fine grid edge:
				EdgeBase * edge = grid->get_child<EdgeBase, TElem> (c_elem,  child);
				
			// Get the local coordinates of the edge center w.r.t. the parent element:
				MathVector<TElem::dim> loc_0, loc_1, loc_center;
				GetRegularLocalCoordinate (grid, edge->vertex (0), c_elem, loc_0);
				GetRegularLocalCoordinate (grid, edge->vertex (1), c_elem, loc_1);
				VecAdd (loc_center, loc_0, loc_1);
				loc_center /= 2;
				
			// Get the length of the edge according to the local coordinates
				MathVector<WDim> corner_0, corner_1, edge_vec;
				DimReferenceMapping<dim, WDim> & map = ReferenceMappingProvider::get<dim, WDim> (roid);
				map.update (aCorners);
				map.local_to_global (corner_0, loc_0);
				map.local_to_global (corner_1, loc_1);
				VecSubtract (edge_vec, corner_1, corner_0);
				
			// Get the shapes and assemble the contribution to the interpolation matrix:
				MathVector<WDim> shape [ref_elem_type::numEdges];
				NedelecT1_LDisc<TDomain, TElem>::get_shapes (domain, c_elem, aCorners, loc_center, shape);
				for (size_t fct = 0; fct < num_fct; fct++)
				{
				//	Get the DoFs:
					if (coarseDD.dof_indices (c_elem, fct, c_ind) != (size_t) ref_elem_type::numEdges
							|| fineDD.inner_dof_indices (edge, fct, f_ind) != 1)
						UG_THROW ("NedelecTransfer:"
							"More than one DoF per edge. Not the Nedelec-type-1 element?");
				//	Add the contributions to the matrix entry:
					for (size_t i = 0; i < (size_t) ref_elem_type::numEdges; i++)
					{
						number coef = VecDot (shape [i], edge_vec);
						DoFRef (mat, f_ind [0], c_ind [i]) += coef;
						vIsRestricted [c_ind [0] [0]] = true;
					}
				}
			}
		}
	}
}

/**
 * Initializes the object: Verifies the data and computes the prolongation
 * matrix.
 */
template <typename TDomain, typename TAlgebra>
void NedelecTransfer<TDomain, TAlgebra>::init ()
{
//	The grid element types for this dimension:
	typedef typename domain_traits<WDim>::AllElemList ElemList;
	
//	Verify the approximation space:
	check_approximation_space ();
	
//	Get the DoF distributions:
	const DoFDistribution & coarseDD = * m_spApproxSpace->dof_distribution (m_coarseLevel);
	const DoFDistribution & fineDD = * m_spApproxSpace->dof_distribution (m_fineLevel);
	
//  Get number of dofs on the grid levels:
	const size_t numFineDoFs = fineDD.num_indices ();
	const size_t numCoarseDoFs = coarseDD.num_indices ();

//	Check if grid distribution has dofs, otherwise skip creation since father
//	elements may not exist in parallel:
	if(numFineDoFs == 0 || numCoarseDoFs == 0)
		return;

//	Assemble the prolongation matrix:
	m_vIsRestricted.clear (); m_vIsRestricted.resize (numCoarseDoFs, false);
	try
	{
		m_prolongation_matrix.resize_and_clear (numFineDoFs, numCoarseDoFs);
	}
	UG_CATCH_THROW ("NedelecTransfer: Failed to allocate the prolongation matrix.");
	boost::mpl::for_each<ElemList> (AssembleProlongationMatrix (this, coarseDD, fineDD));

//	Done:
	m_bInit = true;
}

/**
 * Performs prolongation by applying the prolongation matrix to the given
 * vector.
 */
template <typename TDomain, typename TAlgebra>
void NedelecTransfer<TDomain, TAlgebra>::prolongate
(
	vector_type & uFine, ///< [out] the prolongated vector
	const vector_type & uCoarse ///< [in] the vector to prolongate
)
{
	if (! m_bInit)
		UG_THROW("NedelecTransfer: Operator not initialized.");
	
// Some assertions
	UG_ASSERT(uFine.size () >= m_prolongation_matrix.num_rows (), "Vector ["
		<< uFine.size() << "] must be >= Row size " << m_prolongation_matrix.num_rows ());
	UG_ASSERT(uCoarse.size () >= m_prolongation_matrix.num_cols (), "Vector ["
		<< uCoarse.size() << "] must be >= Col size " << m_prolongation_matrix.num_cols ());

// Apply the prolongation matrix
	if (! m_prolongation_matrix.apply (uFine, uCoarse))
	{
		std::stringstream ss;
		ss << "NedelecTransfer::prolongate: Cannot apply matrix. ";
#ifdef UG_PARALLEL
		ss << "(Type uCoarse = " << uCoarse.get_storage_mask ();
#endif
		UG_THROW (ss.str().c_str ());
	}

// Set Dirichlet nodes to zero again
	try
	{
		for (size_t i = 0; i < m_vConstraint.size (); ++i)
			if (m_vConstraint[i]->type () & CT_DIRICHLET)
				m_vConstraint[i]->adjust_defect (uFine, uFine,
					m_spApproxSpace->dof_distribution (m_fineLevel));
	}
	UG_CATCH_THROW("NedelecTransfer::prolongate: Error while setting dirichlet defect to zero.");

// Call further constraints constraints (= adjust_restrict, member of class constraint)
	try
	{
		for(size_t i = 0; i < m_vConstraint.size (); ++i)
			m_vConstraint[i]->adjust_prolongation (uFine, m_fineLevel, uCoarse, m_coarseLevel);
	}
	UG_CATCH_THROW("NedelecTransfer::prolongate: Error while setting dirichlet defect to zero.");
}

/**
 * Performs prolongation by applying the prolongation matrix to the given
 * vector.
 */
template <typename TDomain, typename TAlgebra>
void NedelecTransfer<TDomain, TAlgebra>::do_restrict
(
	vector_type & uCoarse, ///< [out] the restricted vector
	const vector_type & uFine ///< [in] the vector to restrict
)
{
	if (! m_bInit)
		UG_THROW("NedelecTransfer: Operator not initialized.");
	
	vector_type	uTmp; uTmp.resize(uCoarse.size());

// Some assertions
	UG_ASSERT (uFine.size () >= m_prolongation_matrix.num_rows (), "Vector ["
		<< uFine.size () << "] must be >= Row size " << m_prolongation_matrix.num_rows ());
	UG_ASSERT (uCoarse.size () >= m_prolongation_matrix.num_cols (), "Vector ["
		<< uCoarse.size() << "] must be >= Col size " << m_prolongation_matrix.num_cols ());

// Apply the transposed prolongation matrix
	if (! m_prolongation_matrix.apply_transposed (uTmp, uFine))
		UG_THROW ("NedelecTransfer::restrict: Cannot apply transposed matrix.");

// Copy only restricted values:
//	This is needed for adaptive grids, where the defect that has been
//	projected from the surface should remain and only hidden (i.e.
//	indices with children) should be changed.
	for (size_t i = 0; i < uTmp.size (); ++i)
		if (m_vIsRestricted [i])
			uCoarse[i] = uTmp[i];

// Set dirichlet nodes to zero again
	try
	{
		for(size_t i = 0; i < m_vConstraint.size (); ++i)
			if (m_vConstraint[i]->type () & CT_DIRICHLET)
				m_vConstraint[i]->adjust_defect (uCoarse, uCoarse,
					m_spApproxSpace->dof_distribution (m_coarseLevel));
	}
	UG_CATCH_THROW("NedelecTransfer::restrict: Error while setting dirichlet defect to zero.");

// call restrictions due to added constraints (= adjust_restrict, member of class constraint)
	try
	{
		for (size_t i = 0; i < m_vConstraint.size (); ++i)
			m_vConstraint[i]->adjust_restriction (uCoarse, m_coarseLevel, uFine, m_fineLevel);
	}
	UG_CATCH_THROW("NedelecTransfer::restrict: Error while setting dirichlet defect to zero.");
}

/**
 * Checks the approximation space
 */
template <typename TDomain, typename TAlgebra>
void NedelecTransfer<TDomain, TAlgebra>::check_approximation_space ()
{
	if (m_spApproxSpace.invalid ())
		UG_THROW ("NedelecTransfer: Approximation space not set.");
	
	for (size_t fct = 0; fct < m_spApproxSpace->num_fct (); fct++)
		if (m_spApproxSpace->local_finite_element_id (fct). type() != LFEID::NEDELEC)
			UG_THROW ("NedelecTransfer: "
						"All the functions should be based on the Nedelec element "
						"but function #" << fct << " is not.");
}

/**
 * Sets and checks the grid levels to work on:
 */
template <typename TDomain, typename TAlgebra>
void NedelecTransfer<TDomain, TAlgebra>::set_levels
(
	GridLevel coarseLevel, ///< [in] the coarse grid level
	GridLevel fineLevel ///< [in] the fine grid level
)
{
	m_bInit = false;
	m_fineLevel = fineLevel;
	m_coarseLevel = coarseLevel;

	if(m_fineLevel.level () - m_coarseLevel.level () != 1)
		UG_THROW("NedelecTransfer: Can only project between successive level.");

	if(m_fineLevel.type () != GridLevel::LEVEL || m_coarseLevel.type () != GridLevel::LEVEL)
		UG_THROW("NedelecTransfer: Can only project between level dof distributions, but fine="
				<< m_fineLevel << ", coarse=" << m_coarseLevel);
}

/**
 * Copies parameters from another object
 */
template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferOperator<TAlgebra> > NedelecTransfer<TDomain, TAlgebra>::clone ()
{
	SmartPtr<NedelecTransfer> op (new NedelecTransfer (m_spApproxSpace));
	
	for(size_t i = 0; i < m_vConstraint.size (); ++i)
		op->add_constraint (m_vConstraint[i]);
	
	return op;
}

/**
 * Adds a new constraint to the transfer operators
 */
template <typename TDomain, typename TAlgebra>
void NedelecTransfer<TDomain, TAlgebra>::add_constraint
(
	SmartPtr<IConstraint<TAlgebra> > pp
)
{
//	add only once
	if (std::find (m_vConstraint.begin(), m_vConstraint.end(), pp) != m_vConstraint.end ())
		return;
	m_vConstraint.push_back (pp);
}

/**
 * Removes a constraint from the transfer operators
 */
template <typename TDomain, typename TAlgebra>
void NedelecTransfer<TDomain, TAlgebra>::remove_constraint
(
	SmartPtr<IConstraint<TAlgebra> > pp
)
{
	m_vConstraint.erase (m_vConstraint.begin (),
	                     std::remove (m_vConstraint.begin (), m_vConstraint.end (), pp));
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
