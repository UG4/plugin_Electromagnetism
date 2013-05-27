/**
 * The hybrid smoother proposed by R. Hiptmair for the discretizations
 * based on the Nedelec element. The smoother is implemented for the
 * discretizations of the time-harmonic equations.
 *
 * Created on: 04.12.2012
 * Author: D. Logashenko
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

/// Vertex-centered grid functions (these are GridFunction to allow geometry-dependent smoothers)
///\{
	pot_vector_type * m_pPotDefRe;
	pot_vector_type * m_pPotDefIm;
	
	pot_vector_type m_potCorRe;
	pot_vector_type m_potCorIm;
///\}
	
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

/// Storage for the information about the edge-vertex interconnections
/// \{
	struct tEdgeInfo
	{
		size_t vrt_index [2]; ///< vertix dof's of the beginning and the end of the edge
	};
	VariableArray1<tEdgeInfo> m_vEdgeInfo;
/// \}

/// Flags of the 'conductive' vertices:
	VariableArray1<bool> m_vConductiveVertex;
	
/// Whether initialized
	bool m_bInit;

public:
/// Constructor setting the approx. spaces and the default subsmoothers
	TimeHarmonicNedelecHybridSmoother
	(
		SmartPtr<ApproximationSpace<TDomain> > vertApproxSpace, ///< vertex-centered approx. space
		SmartPtr<ILinearIterator<vector_type> > edgeSmoother, ///< the edge-centered smoother
		SmartPtr<ILinearIterator<pot_vector_type> > vertSmoother ///< the vertex-centered smoother
	)
	: m_spPotMat (new pot_matrix_operator_type),
	  m_pPotDefRe (0), m_pPotDefIm (0),
	  m_spVertApproxSpace (vertApproxSpace),
	  m_spEdgeSmoother (edgeSmoother), m_spVertSmoother (vertSmoother),
	  m_bInit(false)
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
		delete m_pPotDefRe; delete m_pPotDefIm;
	}

/// Returns the name
	const char* name() const {return "Hiptmair hybrid smoother for Whitney-1 elements";}
	
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

	///	Clone the smoother by copying the data
	SmartPtr<ILinearIterator<vector_type, vector_type> > clone ()
	{
		SmartPtr<this_type> newInst
			(new this_type (m_spVertApproxSpace,
				m_spEdgeSmoother->clone (), m_spVertSmoother->clone ()));
		newInst->set_debug (this->debug_writer ());
		newInst->set_damp (this->damping ());
		return newInst;
	}

private: // Auxiliary functions:

///	Computes the matrix for the smoother in the potential space and marks the "conductive nodes"
	void compute_potential_matrix
	(
		const DoFDistribution * pEdgeDD, ///< ptr to the DoF distribution for the edge-centered grid func.
		const DoFDistribution * pVertDD ///< ptr to the DoF distribution for the vertex-centered grid func.
	);

/// Gets the correspondence between the edges and the vertices
	void get_edge_vert_correspondence
	(
		const DoFDistribution * pEdgeDD, ///< ptr to the DoF distribution for the edge-centered grid func.
		const DoFDistribution * pVertDD ///< ptr to the DoF distribution for the vertex-centered grid func.
	);

/// Computes the product \f$ G^T M^{(1)}_h G \f$
	void compute_GtMG ();

/// Computes the vertex-centered defect \f$ d_{pot} G^T d \f$
	void collect_edge_defect
	(
		const vector_type & d ///< the original (edge-centered) defect
	);

/// Adds the vertex-centered correction to the edge-centered one:
	void distribute_vertex_correction
	(
		vector_type & c ///< the final (edge-centered) correction
	);

}; // class TimeHarmonicNedelecHybridSmoother

} // end namespace Electromagnetism
} // end namespace ug

// Implementation of the functions:
#include "hiptmair_hybrid_smoother_impl.h"

#endif /* __H__UG__PLUGINS__ELECTROMAGNETISM__HIPTMAIR_HYBRID_SMOOTHER__ */

/* End of File */
