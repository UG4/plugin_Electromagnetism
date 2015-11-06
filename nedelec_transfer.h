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
 * nedelec_transfer.h - class of the transfer operators for the Nelelec
 * (Whitney-1) elements.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_TRANSFER__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_TRANSFER__

#include "common/common.h"
#include "lib_disc/operator/linear_operator/transfer_interface.h"

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
		Vertex * vrt, ///< [in] the vertex
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
class NedelecProlongationMatrixHelper<TDomain, TAlgebra, RegularEdge>
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
class NedelecTransfer: public ITransferOperator<TDomain, TAlgebra>
{
public:
/// This type
	typedef NedelecTransfer<TDomain, TAlgebra> this_type;

///	Type of base class
	typedef ITransferOperator<TDomain, TAlgebra> base_type;

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
	:	ITransferOperator<TDomain, TAlgebra>(), m_bInit (false), m_spApproxSpace (approxSpace)
		{};
	
///	Set levels
	virtual void set_levels (GridLevel coarseLevel, GridLevel fineLevel);
	
///	initializes the operator (computes the prolongation matrix etc)
	void init ();
	
/// applies the prolongation
	void prolongate (vector_type & uFineOut, const vector_type & uCoarse);
	
/// apples the restriction = transposed prolongation
	void do_restrict (vector_type & uCoarse, const vector_type & uFine);
	
///	returns new instance with same setting
	SmartPtr<ITransferOperator<TDomain, TAlgebra> > clone ();
	
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
	using base_type::m_vConstraint;
	
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
