/**
 * nedelec_transfer.h - class of the transfer operators for the Nelelec
 * (Whitney-1) elements.
 * 
 * Created on: 22.03.2013
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_TRANSFER__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_TRANSFER__

#include "common/common.h"
#include "lib_disc/operator/linear_operator/transfer_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{
namespace Electromagnetism{

/**
 * Helper class for assembling the prolongation matrix
 */
template <typename TDomain, typename TAlgebra, typename TElem>
class NedelecProlongationMatrixHelper
{
///	Type of Domain
	typedef TDomain domain_type;
	
///	Type of algebra
	typedef TAlgebra algebra_type;

///	Type of Vector
	typedef typename TAlgebra::matrix_type matrix_type;

/// world dimention
	static const int WDim = TDomain::dim;

/// reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	
/// dimensionality of the reference element
	static const int dim = ref_elem_type::dim;

/// position type in the domain
	typedef typename TDomain::position_type position_type;

///	computes the local coordinates of a vertex according to the assumption of the regular refinement
	static void GetRegularLocalCoordinate
	(
		const MultiGrid * mg, ///< [in] the grid hierarchy
		VertexBase * vrt, ///< [in] the vertex
		TElem * base, ///< [in] the base element (for the local coordinates)
		MathVector<TElem::dim> & local ///< [out] to save the local coordinates
	);
	
public:

/// assembles the prolongation matrix for one type of the grid elements
	static void assemble_prolongation_matrix
	(
		const domain_type & domain, ///< [in] the domain
		const DoFDistribution & coarseDD, ///< [in] the coarse grid DD
		const DoFDistribution & fineDD, ///< [in] the fine grid DD
		matrix_type & mat, ///< [out] the interpolation matrix
		std::vector<bool> & vIsRestricted ///< [out] whether a coarse grid DoF has children
	);
};

/**
 * Helper class for assembling the prolongation matrix: Instantiation for edges.
 */
template <typename TDomain, typename TAlgebra>
class NedelecProlongationMatrixHelper<TDomain, TAlgebra, Edge>
{
///	Type of Domain
	typedef TDomain domain_type;
	
///	Type of algebra
	typedef TAlgebra algebra_type;

///	Type of Vector
	typedef typename TAlgebra::matrix_type matrix_type;

/// world dimention
	static const int WDim = TDomain::dim;

public:

/// assembles the prolongation matrix for one type of the grid elements
	static void assemble_prolongation_matrix
	(
		const domain_type & domain, ///< [in] the domain
		const DoFDistribution & coarseDD, ///< [in] the coarse grid DD
		const DoFDistribution & fineDD, ///< [in] the fine grid DD
		matrix_type & mat, ///< [out] the interpolation matrix
		std::vector<bool> & vIsRestricted ///< [out] whether a coarse grid DoF has children
	);
};

/**
 * \brief Class of the prolongation and the restriction of the Nedelec DoFs.
 *
 * The prolongation is based on the transformation of the Nelelec DoF values
 * into the vector-valued functions using the basis Nedelec (Whitney-1) shape
 * functions on the coarse grid elements and the computation of the Nedelec DoF
 * values of these vector functions on the fine grid (son) elements. Note that
 * the DoF values are merely the circulations of the vector functions over the
 * edges.
 *
 * The restriction is the transposed prolongation.
 *
 * Note the special case of the "projected boundaries", i.e. when after the
 * regular refinement, new vertices on the boundary sides of the coarse grid
 * elements (or anywhere else) are moved to the actual geometric boundary. In
 * this case, this implementation of the transfer operators works as if it
 * would be applied before moving the new fine grid vertices, i.e. just after
 * the regular refinement. So, the procedure can be considered as "refine
 * regularily - transfer - project the boundary". Due to this approach, the
 * restricted defect belongs to the image of the coarse grid matrix, so that
 * all the systems in the multigrid hierarchy are solvable. Otherwise this
 * may not be the case, cf. O. Sterz. Modellierung und Numerik zeitharmonischer
 * Wirbelstromprobleme in 3D. PhD thesis, 2003.
 *
 * To this end, the evaluations of the Nedelec (Whitney-1) shape functions are
 * performed not according to the geometrical (global) coordinates of the
 * new vertices, but according to their local coordinates in the full-dim.
 * parent element, whereas these local coordinates are computed using averaging
 * of the local coordinates of the corresponding corners. Thus this local
 * coordinates correspond to the positions "before the moving" and not to
 * the actual positions in the space.
 *
 * For grid functions having several Nedelec DoFs on every edge, the
 * prolongation and the restriction is performed separately for every
 * subfunction.
 */
template <typename TDomain, typename TAlgebra>
class NedelecTransfer: public ITransferOperator<TAlgebra>
{
public:
/// This type
	typedef NedelecTransfer<TDomain, TAlgebra> this_type;
	
///	Type of Domain
	typedef TDomain domain_type;
	
///	Type of algebra
	typedef TAlgebra algebra_type;

///	Type of Vector
	typedef typename TAlgebra::vector_type vector_type;

///	Type of Vector
	typedef typename TAlgebra::matrix_type matrix_type;

/// world dimention
	static const int WDim = TDomain::dim;
	
public:

///	Constructor setting approximation space
	NedelecTransfer (SmartPtr<ApproximationSpace<TDomain> > approxSpace)
	:	m_bInit (false), m_spApproxSpace (approxSpace)
		{clear_constraints ();};
	
///	Set levels
	virtual void set_levels (GridLevel coarseLevel, GridLevel fineLevel);
	
///	clears dirichlet post processes
	void clear_constraints () {m_vConstraint.clear();};
	
///	adds a dirichlet post process (not added if already registered)
	void add_constraint (SmartPtr<IConstraint<TAlgebra> > pp);
	
///	removes a post process
	void remove_constraint (SmartPtr<IConstraint<TAlgebra> > pp);
	
///	initializes the operator (computes the prolongation matrix etc)
	void init ();
	
/// applies the prolongation
	void prolongate (vector_type & uFineOut, const vector_type & uCoarse);
	
/// apples the restriction = transposed prolongation
	void do_restrict (vector_type & uCoarse, const vector_type & uFine);
	
///	returns new instance with same setting
	SmartPtr<ITransferOperator<TAlgebra> > clone ();
	
private:
	
/// checks the approximation space
	void check_approximation_space ();
	
/// assembles the prolongation matrix for one type of the grid elements
	template <typename TElem>
	void assemble_prolongation_matrix
	(
		matrix_type & mat, ///< [out] the interpolation matrix
		std::vector<bool> & vIsRestricted, ///< [out] whether a coarse grid DoF has children
		const DoFDistribution & coarseDD, ///< [in] the coarse grid DD
		const DoFDistribution & fineDD ///< [in] the fine grid DD
	);
	
/// a helper class to call all the type-dependent assembling functions
	struct AssembleProlongationMatrix
	{
		this_type * m_pThis;
		const domain_type & m_domain;
		const DoFDistribution & m_coarseDD;
		const DoFDistribution & m_fineDD;
		
		AssembleProlongationMatrix
		(
			this_type* pThis,
			const DoFDistribution & coarseDD,
			const DoFDistribution & fineDD
		)
		:	m_pThis (pThis), m_domain (* m_pThis->m_spApproxSpace->domain().get ()),
			m_coarseDD (coarseDD), m_fineDD (fineDD) {}
		
		template <typename TElem> void operator() (TElem &)
		{
			NedelecProlongationMatrixHelper<TDomain, TAlgebra, TElem>::assemble_prolongation_matrix
				(m_domain, m_coarseDD, m_fineDD,
					m_pThis->m_prolongation_matrix, m_pThis->m_vIsRestricted);
		}
	};
	
private:

///	initialization flag
	bool m_bInit;
	
///	matrix to store prolongation
	matrix_type m_prolongation_matrix;

///	restriction flag
	std::vector<bool> m_vIsRestricted;

///	list of post processes
	std::vector<SmartPtr<IConstraint<TAlgebra> > > m_vConstraint;
	
///	approximation space
	SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;
	
///	fine grid level
	GridLevel m_fineLevel;
	
///	coarse grid level
	GridLevel m_coarseLevel;
};

} // end namespace Electromagnetism
} // end namespace ug

#include "nedelec_transfer_impl.h"

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_TRANSFER__

/* End of File */
