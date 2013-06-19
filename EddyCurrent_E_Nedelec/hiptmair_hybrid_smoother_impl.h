/*
 * Implementation of the class functions templates of the class of the
 * hybrid smoother by Hiptmair.
 */

namespace ug{
namespace Electromagnetism{

//==== Initialization ====

/**
 * Gets the correspondence between the edge DoFs and the vertex DoFs
 * (i.e. fills the tEdgeInfo-structures in m_vEdgeInfo). Furthermore,
 * this function marks edges at Dirichlet boundaries.
 *
 * \tparam TDomain	type of domain in the class
 * \tparam TAlgebra	type of algebra in the class
 *
 * \param[in] pEdgeDD	the current edge-centered DoF distribution of the grid functions
 * \param[in] pVertDD	the currend vertex-centered DoF distribution of the grid functions
 */
template <typename TDomain, typename TAlgebra>
void TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra>::get_edge_vert_correspondence
(
	const DoFDistribution * pEdgeDD,
	const DoFDistribution * pVertDD
)
{
	typedef DoFDistribution::traits<Edge>::const_iterator t_edge_iter;
	
	m_vEdgeInfo.resize (pEdgeDD->num_indices (), false);
	std::vector<size_t> vEdgeInd (1), vVertInd (1);
	
//	Get the edge-to-vertex correspondence
	t_edge_iter edgeIterEnd = pEdgeDD->end<Edge> ();
	for (t_edge_iter edgeIter = pEdgeDD->begin<Edge> (); edgeIter != edgeIterEnd; ++edgeIter)
	{
		EdgeBase * pEdge = *edgeIter;
		if (pEdgeDD->inner_algebra_indices (pEdge, vEdgeInd) != 1)
			UG_THROW (name () << ": Edge DoF distribution mismatch. Not the Nedelec-Type-1 element?");
		tEdgeInfo & EdgeInfo = m_vEdgeInfo [vEdgeInd [0]];
		
		EdgeInfo.Dirichlet = false;
		
		for (size_t i = 0; i < 2; i++)
		{
			VertexBase * pVertex = pEdge->vertex (i);
			if (pVertDD->inner_algebra_indices (pVertex, vVertInd) != 1)
				UG_THROW (name () << ": Vertex DoF distribution mismatch. Not the Lagrange-Order-1 element?");
			EdgeInfo.vrt_index [i] = vVertInd [0];
		}
	}
	
	if (m_spDirichlet.invalid ()) return; // no Dirichlet boundaries specified
	
//	Mark Dirichlet edges
	SubsetGroup dirichlet_ssgrp (pEdgeDD->subset_handler ());
	m_spDirichlet->get_dirichlet_subsets (dirichlet_ssgrp);
	for (size_t j = 0; j < dirichlet_ssgrp.size (); j++)
	{
		int si = dirichlet_ssgrp [j];
	//	Loop the edges in the subset
		t_edge_iter iterEnd = pEdgeDD->end<Edge> (si);
		for (t_edge_iter iter = pEdgeDD->begin<Edge> (si); iter != iterEnd; iter++)
		{
			EdgeBase * pEdge = *iter;
			if (pEdgeDD->inner_algebra_indices (pEdge, vEdgeInd) != 1)
				UG_THROW (name () << ": Edge DoF distribution mismatch. Not the Nedelec-Type-1 element?");
			m_vEdgeInfo[vEdgeInd[0]].Dirichlet = true;
		}
	}
}

/**
 * Performs the multiplication \f$ G^T M^{(1)}_h G \f$. Furthermore,
 * marks the 'conductive' vertices, i.e. the vertices that lie on the closure
 * of the conductive subset (i.e. fills the m_vConductiveVertex array).
 * 'Conductive' vertices are ends of the edges whose entries of the mass matrix
 * \f$ M^{(1)}_h\f$ are non-zero. Any other vertex, which is not an end of
 * such an edge, is considered as 'non-conductive'.
 *
 * For the non-conductive nodes, the function sets the row of the potential
 * stiffness matrix to the identity.
 *
 * This function uses the m_vEdgeInfo-array to get the correspondence between
 * the edges.
 */
template <typename TDomain, typename TAlgebra>
void TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra>::compute_GtMG ()
{
	size_t N_edges = m_spSysMat->num_rows ();
	size_t N_verts = m_spPotMat->num_rows ();
	
	m_spPotMat->set (0.0);
	m_vConductiveVertex.resize (N_verts, false);
	memset (& (m_vConductiveVertex.at (0)), 0, N_verts * sizeof (bool));
	
	// 1. Compute G^T M^{(1)}_h G and mark the 'conductive' vertices.
	// Loop over the rows:
	for (size_t row = 0; row < N_edges; row++)
	{
		bool conductive_edge = false;
		// Remark: For 'conductive edges', the mass matrix is non-zero. We check that.
		
		size_t startVert [2];
		memcpy (startVert, m_vEdgeInfo[row].vrt_index, 2 * sizeof (size_t));
		
		// Loop over the connections:
		for (typename matrix_type::row_iterator row_it = m_spSysMat->begin_row(row);
				row_it != m_spSysMat->end_row(row); ++row_it)
		{
			double M_entry = BlockRef (row_it.value (), 1, 0); // the entry of M^{(1)}_h
			
			// Skip zero entries
			if (zero_mass_entry (M_entry))
				continue; // to prevent extra connections in the potential stiffness matrix
			else
				conductive_edge = true;
			
			size_t destVert [2];
			memcpy (destVert, m_vEdgeInfo[row_it.index()].vrt_index, 2 * sizeof (size_t));
			
			// Update the entries of the potential stiffness matrix:
			(* m_spPotMat) (startVert[0], destVert[0]) += M_entry;
			(* m_spPotMat) (startVert[0], destVert[1]) -= M_entry;
			(* m_spPotMat) (startVert[1], destVert[0]) -= M_entry;
			(* m_spPotMat) (startVert[1], destVert[1]) += M_entry;
		}
		
		// Mark the 'conductive' nodes:
		if (conductive_edge)
			m_vConductiveVertex[startVert[0]] = m_vConductiveVertex[startVert[1]] = true;
	}
	
	// 2. Unmark vertices at the Dirichlet boundary:
	// Loop over the rows:
	for (size_t row = 0; row < N_edges; row++)
	{
		tEdgeInfo & EdgeInfo = m_vEdgeInfo [row];
		if (EdgeInfo.Dirichlet)
			m_vConductiveVertex[EdgeInfo.vrt_index[0]]
			 = m_vConductiveVertex[EdgeInfo.vrt_index[1]] = false;
	}
	
	// 2. Assemble the identity matrix for the 'non-conductive' vertices
	// (otherwise we would have zero lines there).
	// Loop over the vertices:
	for (size_t row = 0; row < N_verts; row++)
		if (! m_vConductiveVertex [row])
			SetDirichletRow (* m_spPotMat.get (), row);
}

/**
 * Computes the matrix \f$A_{pot}\f$ for the smoother in the potential space
 * and marks the "conductive nodes" for the computation of the correction.
 * \f{eqnarray*}{
 *  A_{pot} = G^T M^{(1)}_h G,
 * \f}
 * where \f$G\f$ is the note-to-edge gradient incidence matrix and \f$M^{(1)}_h\f$
 * is the mass matrix of the Whitney-1-forms, taken from the (2, 1)-position of
 * the block-form of the system matrix of the time-harmonic system:
 * \f{eqnarray*}{
 * \left (
 *  \begin{array} {cc}
 *   - S_h (\mu^{-1}) & \omega M^{(1)}_h (\sigma) \\
 *   \omega M^{(1)}_h (\sigma) & S_h (\mu^{-1})
 *  \end{array}
 * \right )
 * \f}
 * (this is the edge-centered matrix).
 *
 * \tparam TDomain	type of domain in the class
 * \tparam TAlgebra	type of algebra in the class
 *
 * \param[in] pEdgeDD	the current edge-centered DoF distribution of the grid functions
 * \param[in] pVertDD	the currend vertex-centered DoF distribution of the grid functions
 */
template <typename TDomain, typename TAlgebra>
void TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra>::compute_potential_matrix
(
	const DoFDistribution * pEdgeDD,
	const DoFDistribution * pVertDD
)
{
	size_t N_edges = pEdgeDD->num_indices ();
	size_t N_verts = pVertDD->num_indices ();
	
	// To be on the safe side:
	if (N_edges != m_spSysMat->num_rows () || N_edges != m_spSysMat->num_cols ())
		UG_THROW (name () << ": smoother not init. or illegal matrix type");
	
	// Allocate the auxiliary matrices and vectors:
	m_spPotMat->resize_and_clear (N_verts, N_verts);
	
	// Get the edge-vertex dof correspondence:
	get_edge_vert_correspondence (pEdgeDD, pVertDD);
	
	// Compute the potential stiffness matrix:
	compute_GtMG ();
}

/**
 *	This function initializes the smoother with the matrix J and the geometry
 *	read from the GridFunction u. Note that u should be a GridFunction, not merely
 *	a vector.
 *
 * \param[in]	J		the matrix of the time-harmonic Nedelec discretization
 * \param[in]	u		the grid function containing the grid info (the values are not used)
 * \returns		bool	success flag
 */
template <typename TDomain, typename TAlgebra>
bool TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra>::init
(
	SmartPtr<ILinearOperator<vector_type> > J,
	const vector_type & u
)
{
	try
	{
	// Check the vertex approx. space:
		if (m_spVertApproxSpace.invalid())
			UG_THROW (name() << ", init: No approx. space for the potential specified.");
		
	// Cast and check if a suitable operator
		m_spSysMat = J.template cast_dynamic<MatrixOperator<matrix_type, vector_type> >();
		if(m_spSysMat.invalid())
			UG_THROW (name() << ", init: No or illegal lin. operator passed (not a matrix?).");
		
	// Try to cast the vector to the grid function and get the edge DoF distribution:
		const TGridFunc * pGridF;
		if ((pGridF = dynamic_cast<const TGridFunc *> (&u)) == 0)
			UG_THROW (name() << ", init: No DoFDistribution specified.");
		m_spEdgeDD = ((TGridFunc *) pGridF)->dof_distribution ();
		m_GridLevel = m_spEdgeDD->grid_level ();
		m_spVertDD = m_spVertApproxSpace->dof_distribution (m_GridLevel);
		
	// Compute the matrix of the potential equation:
		compute_potential_matrix (m_spEdgeDD.get (), m_spVertDD.get ());
	
	// Allocate the grid functions for the defect in the potential space:
		delete m_pPotCorRe; m_pPotCorRe = 0;
		delete m_pPotCorIm; m_pPotCorIm = 0;
		m_pPotCorRe = new TPotGridFunc (m_spVertApproxSpace, m_spVertDD);
		m_pPotCorIm = new TPotGridFunc (m_spVertApproxSpace, m_spVertDD);
	
	// Initialize the subordinated smoother for the vertex dof:
		*m_pPotCorRe = 0;
		if (! m_spVertSmoother->init (m_spPotMat, *m_pPotCorRe))
			UG_THROW (name() << ", init: Failed to initialize the subordinated vertex-based smoother.");
		
	// Initialize the subordinated smoother for the edge dofs:
		if (! m_spEdgeSmoother->init (J, u))
			UG_THROW (name() << ", init: Failed to initialize the subordinated edge-based smoother.");
	}
	catch (...)
	{
		m_bInit = false;
		delete m_pPotCorRe; m_pPotCorRe = 0;
		delete m_pPotCorIm; m_pPotCorIm = 0;
		m_vEdgeInfo.resize (0, false);
		m_vConductiveVertex.resize (0, false);
		m_spPotMat->resize_and_clear (0, 0);
		throw;
	}
	
	m_bInit = true;
	return true;
}

//==== Computation of the correction ====

/**
 * Computes the defect in the potential space: \f$ d_{pot} = G^T d \f$, where
 * \f$ d \f$ is the defect of the original system. This is done only for the
 * 'conductive' vertices. At the other vertices, \f$ d_{pot} \f$ is set to 0.
 *
 * \param[in]	d			the original (edge-centered) defect
 * \param[out]	potDefRe	real part of \f$ d_{pot} \f$
 * \param[out]	potDefIm	imaginary part of \f$ d_{pot} \f$
 */
template <typename TDomain, typename TAlgebra>
void TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra>::collect_edge_defect
(
	const vector_type & d,
	pot_vector_type & potDefRe,
	pot_vector_type & potDefIm
)
{
	size_t N_edges = m_vEdgeInfo.size ();
	
	potDefRe = 0; potDefIm = 0;
	
	// Loop over the edges:
	for (size_t edge = 0; edge < N_edges; edge++)
	{
		size_t i;
		tEdgeInfo & EdgeInfo = m_vEdgeInfo [edge];
		number Re_d = BlockRef (d [edge], 0);
		number Im_d = BlockRef (d [edge], 1);
		
		if (m_vConductiveVertex [i = EdgeInfo.vrt_index[0]])
		{
			potDefRe [i] -= Re_d;
			potDefIm [i] -= Im_d;
		}
		
		if (m_vConductiveVertex [i = EdgeInfo.vrt_index[1]])
		{
			potDefRe [i] += Re_d;
			potDefIm [i] += Im_d;
		}
	}
}

/**
 * Updates the edge-centered correction: \f$ c := c_{edge} + G c_{pot} \f$.
 *
 * \param[in]	potCorRe	real part of \f$ c_{pot} \f$
 * \param[in]	potCorIm	imaginary part of \f$ c_{pot} \f$
 * \param[out]	c			the final, edge-centered correction
 */
template <typename TDomain, typename TAlgebra>
void TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra>::distribute_vertex_correction
(
	pot_vector_type & potCorRe,
	pot_vector_type & potCorIm,
	vector_type & c
)
{
	size_t N_edges = m_vEdgeInfo.size ();
	
	// Loop over the edges:
	for (size_t edge = 0; edge < N_edges; edge++)
	{
		tEdgeInfo & EdgeInfo = m_vEdgeInfo [edge];
		
		BlockRef (c [edge], 0) += potCorRe [EdgeInfo.vrt_index[1]]
									- potCorRe [EdgeInfo.vrt_index[0]];
		
		BlockRef (c [edge], 1) += potCorIm [EdgeInfo.vrt_index[1]]
									- potCorIm [EdgeInfo.vrt_index[0]];
	}
}

/**
 * This method applies the iterator operator, i.e. invertes the approximate
 * matrix. The defect d remains unchanged.
 *
 * The idea is to apply first the edge-centered smoother (specified as the
 * argument of the constructor of the class). To the correction obtained in
 * this way, the vertex-centered correction is added. The latter is obtained
 * by the application of the vertex-centered smoother to the system
 * \f{eqnarray*}{
 * \left (
 *  \begin{array} {cc}
 *   0 & A_{pot} \\
 *   A_{pot} & 0
 *  \end{array}
 * \right )
 * \left (
 *  \begin{array} {c}
 *   \mathbf{Re} c_{pot} \\
 *   \mathbf{Im} c_{pot}
 *  \end{array}
 * \right )
 * =
 * \left (
 *  \begin{array} {c}
 *   - \mathbf{Re} d_{pot} \\
 *   \mathbf{Im} d_{pot}
 *  \end{array}
 * \right ),
 * \f}
 * where \f$ A_{pot} = G^T M^{(1)}_h G \f$ and \f$ d_{pot} = G^T d \f$. To
 * add the vertex-centered correction, obtained in this way, to the edge-centered
 * one, one sets: \f$ c := c_{edge} + G c_{pot} \f$.
 *
 * \param[in]	d		defect
 * \param[out]	c		correction
 * \returns		bool	success flag
 */
template <typename TDomain, typename TAlgebra>
bool TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra>::apply
(
	vector_type & c,
	const vector_type & d
)
{
//---- Protection:

//	Check that object is initialized:
	if(! m_bInit)
	{
		UG_LOG("ERROR in '"<<name()<<"::apply': Object not initialized.\n");
		return false;
	}

//	Check sizes:
	if (d.size() != m_spSysMat->num_rows())
		UG_THROW("ERROR in '"<<name()<<"::apply': Vector [size= "<<d.size()<<
					"] and Row [size= "<<m_spSysMat->num_rows()<<"] sizes have to match!");
	if (c.size() != m_spSysMat->num_cols())
		UG_THROW("ERROR in '"<<name()<<"::apply': Vector [size= "<<c.size()<<
					"] and Column [size= "<<m_spSysMat->num_cols()<<"] sizes have to match!");
	if (d.size() != c.size())
		UG_THROW("ERROR in '"<<name()<<"::apply': Vector [d size= "<<d.size()<<
					", c size = " << c.size() << "] sizes have to match!");

//	Temporary vectors:
	vector_type auxEdgeDef (d.size ());
	pot_vector_type potDefRe (m_pPotCorRe->size ());
	pot_vector_type potDefIm (m_pPotCorIm->size ());
	
//---- Numerics:

	auxEdgeDef = d;
	
//	Apply the edge-based smoother:
	if (! m_bSkipEdge)
	{
		if (! m_spEdgeSmoother->apply_update_defect (c, auxEdgeDef))
			return false;
	}
	else c = 0.0;
	
	if (! m_bSkipVertex)
	{
	//	Compute the 'vertex defect' of the potential:
		collect_edge_defect (auxEdgeDef, potDefRe, potDefIm);
		
	//	Apply the vertex-centered smoother:
		if (! m_spVertSmoother->apply (*m_pPotCorRe, potDefIm))
			return false;
		if (! m_spVertSmoother->apply (*m_pPotCorIm, potDefRe))
			return false;	
		
	//	Add the vertex-centered correction into the edge-centered one:
		distribute_vertex_correction (*m_pPotCorRe, *m_pPotCorIm, c);
	}
	
//	Damping (the standard way, like for IPreconditioner):
	number kappa = this->damping()->damping (c, d, m_spSysMat);
	if (kappa != 1.0)
		c *= kappa;

	return true;
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
