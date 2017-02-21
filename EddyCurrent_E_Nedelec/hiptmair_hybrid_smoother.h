/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

/*
 * The hybrid smoother proposed by R. Hiptmair for the discretizations
 * based on the Nedelec element. The smoother is implemented for the
 * discretizations of the time-harmonic equations.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__HIPTMAIR_HYBRID_SMOOTHER__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__HIPTMAIR_HYBRID_SMOOTHER__

// basic ug4 headers
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library-specific ug4 headers
#include "lib_algebra/operator/interface/linear_iterator.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/damping.h"
#include "lib_algebra/operator/debug_writer.h"

// plugin headers
#include "../em_material.h"

namespace ug{
namespace Electromagnetism{

/// The hybrid smoother by R. Hiptmair
/**
 * This class implements the hybrid smoother proposed by R. Hiptmair for the
 * discretizations based on the Nedelec element. The smoother is based on the
 * approximate Helmholtz decomponstion of the operator. This implementation is
 * constructed for the discretizations of the time-harmonic equations. The class
 * computes the matrix for the discrete potential and calls two further smoothers
 * specified by the user: one in the space of the edge elements (Whitney-1 or
 * Nedelec forms) and one in the vertex-centered potential space (i.e. the Whitney-0
 * or Lagrange-1 elements). Accordingly, this smoother requires 2 approximation
 * spaces: one for the edge dofs and one for the vertex dofs.
 *
 * References:
 * <ul>
 * <li> R. Hiptmair. Multigrid method for Maxwell’s equations. SIAM J. Numer. Anal., 36(1): pp. 204–225 (1998), DOI 10.1137/S0036142997326203.
 * <li> R. Beck, P. Deuflhard, R. Hiptmair, R. Hoppe und B. Wohlmuth. Adaptive multilevel methods for edge element discretizations of Maxwell’s equations. Surveys on Mathematics for Industry, 8(3-4): pp. 271–312 (1999).
 * </ul>
 *
 * \tparam		TDomain			Type of Domain
 * \tparam		TAlgebra		Type of Algebra
 */
template <typename TDomain, typename TAlgebra>
class TimeHarmonicNedelecHybridSmoother :
	public ILinearIterator<typename TAlgebra::vector_type>,
	public DebugWritingObject<TAlgebra>
{
	typedef TimeHarmonicNedelecHybridSmoother<TDomain,TAlgebra> this_type;
	
public:
///	Vector type
	typedef typename TAlgebra::vector_type vector_type;
///	Matrix type
	typedef typename TAlgebra::matrix_type matrix_type;
///	Matrix Operator type
	typedef MatrixOperator<matrix_type, vector_type> matrix_operator_type;

/// The auxiliary algebra type for the space of the potential. (Note: It should be scalar.)
	typedef CPUAlgebra TPotAlgebra;
///	Vector type for the potential space
	typedef typename TPotAlgebra::vector_type pot_vector_type;
/// Matrix type for the potential space
	typedef typename TPotAlgebra::matrix_type pot_matrix_type;
///	Matrix Operator type for the potential space
	typedef MatrixOperator<pot_matrix_type, pot_vector_type> pot_matrix_operator_type;
	
/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> TGridFunc;
/// Grid function type for the potential
	typedef GridFunction<TDomain, TPotAlgebra> TPotGridFunc;

private:

/// Measure for numerically zero entries in the mass matrix
	bool zero_mass_entry (double val) {return fabs (val) < 1e-64;} // this could be merely 0

/// Edge-centered matrix of the original equation
	SmartPtr<matrix_operator_type> m_spSysMat;
	
/// Vertex-centered matrix for the potential
	SmartPtr<pot_matrix_operator_type> m_spPotMat;
	
/// Vertex-centered grid function for the Re-part of potential corrections (this is a GridFunction to allow geometry-dependent smoothers)
	pot_vector_type * m_pPotCorRe;
/// Vertex-centered grid function for the Im-part of potential corrections (this is a GridFunction to allow geometry-dependent smoothers)
	pot_vector_type * m_pPotCorIm;
	
/// DoF distribution for the Nedelec elements
	SmartPtr<DoFDistribution> m_spEdgeDD;
	
/// Approximation space for the potential (vertex-centered) space
	SmartPtr<ApproximationSpace<TDomain> > m_spVertApproxSpace;
/// Level DoF distribution for the vertex centered grid func.
	SmartPtr<DoFDistribution> m_spVertDD;
	
/// Grid level (or surface) where to smooth
	GridLevel m_GridLevel;
	
/// Smoother for the edge dofs
	SmartPtr<ILinearIterator<vector_type> > m_spEdgeSmoother;
	
/// Smoother for the vertex dofs
	SmartPtr<ILinearIterator<pot_vector_type> > m_spVertSmoother;

///	Dirichlet boundary
	SmartPtr<EMDirichlet<TDomain, TAlgebra> > m_spDirichlet;
	
/// Structure of the storage for the information about the edge-vertex interconnections
	struct tEdgeInfo
	{
		size_t vrt_index [2]; ///< vertex dof's of the beginning and the end of the edge
		
		char flags; ///< flags: (0) Dirichlet, (1) conductive vertex 0, (2) conductive vertex 1, (3) init-ed
		
		bool is_Dirichlet () {return (flags & 1) != 0;}
		bool conductive_vrt_0 () {return (flags & 2) != 0;}
		bool conductive_vrt_1 () {return (flags & 4) != 0;}
		bool is_init () {return (flags & 8) != 0;}
		
		void clear_flags () {flags = 0;}
		void set_Dirichlet () {flags |= 1;}
		void set_conductive_vrt_0 () {flags |= 2;}
		void set_conductive_vrt_1 () {flags |= 4;}
		void set_init () {flags |= 8;}
	};
/// Storage for the information about the edge-vertex interconnections
	VariableArray1<tEdgeInfo> m_vEdgeInfo;

/// Whether initialized
	bool m_bInit;
	
/// Needed mainly for debugging: Whether to skip one of the stages
	bool m_bSkipEdge, m_bSkipVertex;

public:
/// Constructor setting the approx. spaces and the default subsmoothers
	TimeHarmonicNedelecHybridSmoother
	(
		SmartPtr<ApproximationSpace<TDomain> > vertApproxSpace, ///< vertex-centered approx. space
		SmartPtr<ILinearIterator<vector_type> > edgeSmoother, ///< the edge-centered smoother
		SmartPtr<ILinearIterator<pot_vector_type> > vertSmoother ///< the vertex-centered smoother
	)
	: m_spPotMat (new pot_matrix_operator_type),
	  m_pPotCorRe (NULL), m_pPotCorIm (NULL),
	  m_spVertApproxSpace (vertApproxSpace),
	  m_spEdgeSmoother (edgeSmoother), m_spVertSmoother (vertSmoother),
	  m_bInit (false),
	  m_bSkipEdge (false), m_bSkipVertex (false)
	{
		if (m_spVertApproxSpace.invalid ())
			UG_THROW (name() << ": illegal vert.-centered approx. space");
		if (m_spEdgeSmoother.invalid ())
			UG_THROW (name() << ": illegal edge-centered smoother");
		if (m_spVertSmoother.invalid ())
			UG_THROW (name() << ": illegal vert.-centered smoother");
	}
	
/// Destructor
	~TimeHarmonicNedelecHybridSmoother ()
	{
		delete m_pPotCorRe; delete m_pPotCorIm;
	}

/// Returns the name
	const char* name() const {return "Hiptmair hybrid smoother for Whitney-1 elements";}
	
///	Currently returns false because the computation of the potentials is not purely parallel up to now
	bool supports_parallel() const {return true;}
	
/// Initialization using a matrix and a GridFunction (not merely a vector!)
	bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

/// We cannot initialize without the geometry (this version of init cannot not work)
	bool init(SmartPtr<ILinearOperator<vector_type> > L)
	{
		UG_THROW(name () << ": Cannot initialize the hybrid smoother without the geometry. Specify the 2nd argument for init!");
		return false;
	}

///	Computes the correction
	bool apply (vector_type & c, const vector_type & d);

///	Computes the correction and updates the defect d := d - L*c
	bool apply_update_defect (vector_type & c, vector_type & d)
	{
	//	compute the correction
		if(! apply (c, d)) return false;

	// 	update the defect d := d - A*c
		if(! m_spSysMat->matmul_minus (d, c))
			UG_THROW (name() << "::apply_update_defect: failed to execute matmul_minus");

		return true;
	}

///	Sets the Dirichlet boundary
	void set_Dirichlet
	(
		SmartPtr<EMDirichlet<TDomain, TAlgebra> > spDirichlet ///< the Dirichlet BC object
	)
	{
		m_spDirichlet = spDirichlet;
	}
	
///	Skip flag the edge smoother
	void set_skip_edge_smoother (bool skip_edge) {m_bSkipEdge = skip_edge;}
///	Skip flag the vertex smoother
	void set_skip_vertex_smoother (bool skip_vertex) {m_bSkipVertex = skip_vertex;}

///	Clone the smoother by copying the data
	SmartPtr<ILinearIterator<vector_type, vector_type> > clone ()
	{
		SmartPtr<this_type> newInst
			(new this_type (m_spVertApproxSpace,
				m_spEdgeSmoother->clone (), m_spVertSmoother->clone ()));
		newInst->set_Dirichlet (m_spDirichlet);
		newInst->set_skip_edge_smoother (m_bSkipEdge);
		newInst->set_skip_vertex_smoother (m_bSkipVertex);
		newInst->set_debug (this->debug_writer ());
		newInst->set_damp (this->damping ());
		return newInst;
	}
	
private: // Auxiliary functions:

///	Computes the matrix for the smoother in the potential space and marks the "conductive nodes"
	void compute_potential_matrix
	(
		const DoFDistribution * pEdgeDD, ///< edge-centered DoF distribution of the grid functions
		const DoFDistribution * pVertDD ///< vertex-centered DoF distribution of the grid functions
	);

/// Gets the correspondence between the edges and the vertices
	void get_edge_vert_correspondence
	(
		const DoFDistribution * pEdgeDD, ///< edge-centered DoF distribution of the grid functions
		const DoFDistribution * pVertDD ///< vertex-centered DoF distribution of the grid functions
	);

/// Computes the product \f$ G^T M^{(1)}_h G \f$
	void compute_GtMG ();

/// Computes the vertex-centered defect \f$ d_{pot} = G^T d \f$
	void collect_edge_defect
	(
		const vector_type & d, ///< original (edge-centered) defect
		pot_vector_type & potDefRe, ///< real part of \f$ d_{pot} \f$
		pot_vector_type & potDefIm ///< imaginary part of \f$ d_{pot} \f$
	);

/// Adds the vertex-centered correction to the edge-centered one:
	void distribute_vertex_correction
	(
		pot_vector_type & potCorRe, ///< real part of the potential correction \f$ c_{pot} \f$
		pot_vector_type & potCorIm, ///< imaginary part of the potential correction \f$ c_{pot} \f$
		vector_type & c ///< final (edge-centered) correction
	);
	
#ifdef UG_PARALLEL
///	"or" reduction operation class for the conductivity condition
	struct t_red_op_or
	{
		static inline bool op (bool a, bool b) {return a || b;}
	};
///	"and" reduction operation class for the conductivity condition
	struct t_red_op_and
	{
		static inline bool op (bool a, bool b) {return a && b;}
	};
#endif

}; // class TimeHarmonicNedelecHybridSmoother

} // end namespace Electromagnetism
} // end namespace ug

// Implementation of the functions:
#include "hiptmair_hybrid_smoother_impl.h"

#endif /* __H__UG__PLUGINS__ELECTROMAGNETISM__HIPTMAIR_HYBRID_SMOOTHER__ */

/* End of File */
