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
 * nedelec_source.h - classes for computation of divergence-free sources
 * in form of Nedelec-type-1 element DoFs.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_SOURCE__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_SOURCE__

#include "common/common.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/operator/linear_operator/transfer_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/debug_writer.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{
namespace Electromagnetism{

/// Class for computation of loop currents
/**
 * Class for computation of weakly divergence-free loop currents that are
 * non-zero only in a given subdomain \f$T\f$ (i.e. in a coil).
 *
 * The range (image) of the stiffness matrix of the discrete
 * \f$\mathbf{rot} \, \mathbf{rot}\f$ operator is the subspace of all
 * functions othrogonal to all the gradients \f$ \mathbf{G} \psi \f$, where
 * \f$ \mathbf{G} \f$ is the incidence matrix between the vertices and the edges,
 * and \f$ \psi \f$ is any vertex-centered grid function. Thus, for the
 * solvability of the discretized system (for ex. of the E-based formulation of
 * the eddy-current model), its right-hand side \f$M^{(1)} J \f$ (\f$M^{(1)}\f$
 * being the mass matrix of the Whitney-1 shapes) must satisfy
 * \f{eqnarray*}{
 *  ( \mathbf{G}^T M^{(1)} J, \psi ) = ( M^{(1)} J, \mathbf{G} \psi ) = 0
 * \f}
 * This is really true only in insulators. In conductors, a mass matrix is added
 * to the \f$\mathbf{rot} \, \mathbf{rot}\f$ operator, so that the kernel is
 * trivial. We assume that a coil is modelled as a zero-conductivity medium with
 * the source, so that \f$T\f$ is completely contained in insulators.
 * Thus the former condition implies
 * \f{eqnarray*}{
 *  \mathbf{G}_{i}^T M^{(1)} J = 0,
 * \f}
 * (with \f$ J = 0 \f$ out of \f$ T \f$) i.e. the weak divergence of the rhs
 * must be zero.
 *
 * This class implements a computation of such a weakly divergence-free
 * current \f$J\f$ in a given set of insulators \f$T\f$ constituting (topologically)
 * a ring (torus) by specifying a cut with the jump of the potential and a
 * part adjacent to this cut that specifies the positive direction of the
 * current. The current is represented as \f$J = \hat{\mathbf{G}} \phi\f$ and the
 * potential \f$\phi\f$ is computed as the solution of a discretized Laplace
 * equation with the special boundary conditions. Here, \f$\hat{\mathbf{G}}\f$
 * is the matrix obtained from \f$\mathbf{G}\f$ by replacing the rows
 *
 * The computed \f$J\f$ is rescaled so that the current over the cut is
 * equal to the given value. Different values of the currents can be used for
 * different function names: cf. the 'set' functions of this class.
 *
 * Note that the mass matrix \f$M^{(1)}\f$ is not the same as the mass matrix
 * at the off-diagonal positions of the stiffness matrix of the time-harmonic
 * problem. The matrix \f$M^{(1)}\f$ mentioned here for the computation of the
 * right-hand side is assembled only of the local mass matrices for the grid
 * elements from \f$T\f$. The resulting discretized system for \f$\phi\f$ consists
 * the discrete Laplace operator on \f$T\f$, the discrete Neumann-0 (natural)
 * boundary condition on \f$\partial{T}\f$ and the identity matrix for all other
 * DoFs. Note furthermore that the same \f$M^{(1)}\f$ should be used in the
 * assembling of the right-hand side of the discretized Maxwell equations.
 *
 * References:
 * <ul>
 * <li> O. Sterz. Modellierung und Numerik zeitharmonischer Wirbelstromprobleme in 3D. PhD thesis, 2003.
 * </ul>
 *
 * \tparam	TDomain		Domain type
 * \tparam	TAlgebra	Algebra type
 */
template <typename TDomain, typename TAlgebra>
class NedelecLoopCurrent
{
/// This type
	typedef NedelecLoopCurrent<TDomain, TAlgebra> this_type;
	
public:
///	Type of Domain
	typedef TDomain domain_type;
	
///	Type of Grid:
	typedef typename TDomain::grid_type grid_type;
	
/// Type of subset handler
	typedef typename domain_type::subset_handler_type subset_handler_type;
	
///	Type of algebra (for the Nedelec-element-based grid functions)
	typedef TAlgebra algebra_type;
///	Type of Vector (for the Nedelec-element-based grid functions)
	typedef typename TAlgebra::vector_type vector_type;
///	Type of Vector (for the Nedelec-element-based grid functions)
	typedef typename TAlgebra::matrix_type matrix_type;

/// The auxiliary algebra type for the space of the potential. (Note: It should be scalar.)
	typedef CPUAlgebra TPotAlgebra;
///	Vector type for the potential space
	typedef typename TPotAlgebra::vector_type pot_vector_type;
/// Matrix type for the potential space
	typedef typename TPotAlgebra::matrix_type pot_matrix_type;
///	Grid function type for the potential space
	typedef GridFunction<TDomain, TPotAlgebra> pot_gf_type;

/// world dimention
	static const int WDim = TDomain::dim;
	
///	position type
	typedef typename TDomain::position_type position_type;
	
private:

///	Type of the attachment and its accessor for flags
	typedef AChar a_edge_flag_type;
	typedef Grid::EdgeAttachmentAccessor<a_edge_flag_type> aa_edge_flag_type;
	typedef ABool a_vert_flag_type;
	typedef Grid::VertexAttachmentAccessor<a_vert_flag_type> aa_vert_flag_type;
	
public:

///	Constructor
	NedelecLoopCurrent
	(
		const char * ssNames, ///< names of the subsets of the source (up to the subset of the pos. dir.)
		const char * posSsNames, ///< names of the subsets of the positive direction
		const char * cutSsNames, ///< names of the surfaces on the cut of the loop
		SmartPtr<ApproximationSpace<TDomain> > vertApproxSpace, ///< vertex-centered approx. space
		SmartPtr<ILinearOperatorInverse<pot_vector_type> > potSolver ///< linear solver for the potential
	);
	
///	Setting the electric current value
	void set
	(
		const char * fctNames, ///< names of the components
		number I ///< electric current for the functions
	);

///	Computation of the source (only if 'init'ialized): distributes the potential
	void compute
	(
		SmartPtr<GridFunction<TDomain, TAlgebra> > sp_u ///< the grid function for the source
	);
	
///	Returns the subsets where the source is defined
	std::string subsets () {return m_allSsNames;};
	
private:
	
	/// Marks edges that belong to the loop source domain (for one type of elements)
	template <typename TElem>
	void mark_source_edges
	(
		const DoFDistribution & edgeDD, ///< [in] the edge DD
		aa_edge_flag_type & in_source ///< [out] the flags to update
	);
	/// Helper class for marking the edges in the source
	struct MarkSourceEdges
	{
		this_type * m_pThis;
		const DoFDistribution & m_edgeDD;
		aa_edge_flag_type & m_in_source;
		
		MarkSourceEdges
		(
			this_type * pThis, ///< [in] pointer to the master class
			const DoFDistribution & edgeDD, ///< [in] the edge DD
			aa_edge_flag_type & in_source ///< [out] the flags to update
		)
		: m_pThis (pThis), m_edgeDD (edgeDD), m_in_source (in_source) {}
		
		template <typename TElem> void operator() (TElem &)
		{
			m_pThis->template mark_source_edges<TElem> (m_edgeDD, m_in_source);
		}
	};
	
	/// Computes the potential of the source
	void compute_potential
	(
		pot_gf_type & pot_u ///< a grid function for the potential
	);
	
	/// Computes the gradient of the potential
	void distribute_source_potential
	(
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		pot_vector_type & src_pot, ///< [in] potential to distribute
		DoFDistribution & edgeDD, ///< [in] the edge DD
		size_t func, ///< [in] index of the function
		number value, ///< [in] value of the source
		vector_type & src_field ///< [out] the computed source field
	);
	
	/// Helper class for assembling the local stiffness matrix of Laplacian in various routines
	template <typename TElem>
	class LocLaplaceA
	{
	public:
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
		static const int numCorners = ref_elem_type::numCorners;
		
		/// Computes the local discretization of the Laplace operator
		static void stiffness
		(
			GridObject * elem, ///< [in] element to compute for
			const position_type vCornerCoords [], ///< [in] coordinates of the corners of the element
			number loc_A [numCorners] [numCorners] ///< [out] the local stiffness matrix
		);
		
		/// Checks whether corners are on boundary (returns true if yes)
		static bool bnd_corners
		(
			GridObject * elem, ///< [in] element to compute for
			const SubsetGroup & bndGrp, ///< [in] boundary subsets
			bool bnd_flag [numCorners] ///< [out] true for bnd nodes, false for others
		);
	};
	
	/// Class for local assembling of the auxiliary Laplace operator
	/**
	 * This class assembles the local discretization of the Laplace operator
	 * in the closure of the subdomain of the source. (Out of this closure,
	 * the constraint of class OutOfSource is used.)
	 */
	class AuxLaplaceLocAss : public IElemDisc<TDomain>
	{
	private:
	///	a name for the component in the grid functions
		static const size_t _C_ = 0;
		
	/// maximum number of corners of an element
		static const size_t maxCorners = domain_traits<WDim>::MaxNumVerticesOfElem;
		
	public:
	///	constructor
		AuxLaplaceLocAss
		(
			NedelecLoopCurrent & master ///< [in] the master object
		);
		
	///	switches between non-regular and regular grids
		virtual void prepare_setting
		(
			const std::vector<LFEID> & vLfeID,
			bool bNonRegular
		);
	
	private:
	/// assembling functions
	/// \{
		template <typename TElem>
		void prepare_element_loop (ReferenceObjectID roid, int si);

		template <typename TElem>
		void ass_JA_elem (LocalMatrix & J, const LocalVector & u, GridObject * elem, const position_type vCornerCoords []);
		
		template <typename TElem>
		void ass_rhs_elem (LocalVector & d, GridObject * elem, const position_type vCornerCoords []);

		template <typename TElem>
		void prepare_element (const LocalVector & u, GridObject * elem, ReferenceObjectID roid, const position_type vCornerCoords [])
			{};
		template <typename TElem>
		void finish_element_loop ()
			{};
		template <typename TElem>
		void ass_JM_elem (LocalMatrix & J, const LocalVector & u, GridObject * elem, const position_type vCornerCoords [])
			{
				UG_THROW ("NedelecProject: Attempt to use nonstationary discretization for the potential.");
			};
		template <typename TElem>
		void ass_dA_elem (LocalVector & d, const LocalVector & u, GridObject * elem, const position_type vCornerCoords [])
			{
				UG_THROW ("NedelecProject: Attempt to assemble non-linear system for the potential.");
			};
		template <typename TElem>
		void ass_dM_elem (LocalVector & d, const LocalVector & u, GridObject * elem, const position_type vCornerCoords [])
			{
				UG_THROW ("NedelecProject: Attempt to use nonstationary discretization for the potential.");
			};
	/// \}
	
	/// registration of the assembling functions
	/// \{
		void register_all_loc_discr_funcs();

		struct RegisterLocalDiscr {
				RegisterLocalDiscr (AuxLaplaceLocAss* pThis) : m_pThis (pThis) {}
				AuxLaplaceLocAss* m_pThis;
				template< typename TElem > void operator () (TElem&)
				{m_pThis->register_loc_discr_func<TElem> ();}
		};

		template <typename TElem>
		void register_loc_discr_func();
	/// \}

	private:
		
	///	the master object of the projection
		NedelecLoopCurrent & m_master;
		
	/// group of the positive direction subsets
		SubsetGroup m_posSsGrp;
	/// group of the surfaces at the cut
		SubsetGroup m_cutSsGrp;
	
	///	whether we assemble in the subset specifying the positive direction
		bool m_in_pos_subset;
	};
	
///	Constraint that sets the problem to 0-identity out of the source
/**
 *	This class restricts the discretized Laplace equation to the subsets
 *	of the source: In the parts of the domain that do not belong to the
 *	closure of the source, the discretization matrix is set to identity,
 *	and the right-hand side is set to zero.
 */
	class OutOfSource : public IDomainConstraint<TDomain, TPotAlgebra>
	{
	private:
	/// maximum number of corners of an element
		static const size_t maxCorners = (size_t) element_list_traits<typename domain_traits<WDim>::DimElemList>::maxCorners;
		
	/// Iterator over vertices
		typedef DoFDistribution::traits<Vertex>::const_iterator t_vert_iterator;
	
	private:
	
	/// marks the vertices belonging to elements in the source (for one element type)
		template <typename TElem>
		void mark_source_vertices_elem_type
		(
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			aa_vert_flag_type & in_source ///< [out] the flags to update
		);
	/// Helper class for marking the vertices in the closure of the source subdomain
		struct MarkSourceVertices
		{
			OutOfSource * m_pThis;
			const DoFDistribution & m_vertDD;
			aa_vert_flag_type & m_in_source;
			
			MarkSourceVertices
			(
				OutOfSource * pThis, ///< [in] pointer to the master class
				const DoFDistribution & vertDD, ///< [in] the vertex DD
				aa_vert_flag_type & in_source ///< [out] theflags to update
			)
			: m_pThis (pThis), m_vertDD (vertDD), m_in_source (in_source) {}
			
			template <typename TElem> void operator() (TElem &)
			{
				m_pThis->template mark_source_vertices_elem_type<TElem> (m_vertDD, m_in_source);
			}
		};
	/// marks the vertices belonging to elements in the source (for all element types)
		void mark_source_vertices
		(
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			aa_vert_flag_type & in_source ///< [out] the flags
		);
	
	///	sets to identity all the matrix rows that do not belong to the closure of the source subdomain
		void adjust_matrix
		(
			const DoFDistribution & vertDD, ///< the vertex DD
			aa_vert_flag_type & in_source, ///< the flags
			pot_matrix_type & A ///< the matrix to adjust
		);
		
	///	sets to 0 all the entries of a vector that do not belong to the closure of the source subdomain
		void adjust_vector
		(
			const DoFDistribution & vertDD, ///< the vertex DD
			aa_vert_flag_type & in_source, ///< the flags
			pot_vector_type & u ///< the vector to adjust
		);
		
	public:
	///	constructor
		OutOfSource
		(
			NedelecLoopCurrent & master ///< [in] the master object
		);
		
	/// sets a unity row for all conductor indices
		void adjust_jacobian
		(
			pot_matrix_type & J,
			const pot_vector_type & u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol = SPNULL,
			const number s_a0 = 1.0
		)
		{
			MultiGrid * mg = (MultiGrid *) dd->multi_grid().get ();
			a_vert_flag_type a_in_source;
			mg->attach_to_vertices (a_in_source);
			aa_vert_flag_type aa_in_source (*mg, a_in_source);
			mark_source_vertices (* dd.get (), aa_in_source);
			adjust_matrix (* dd.get (), aa_in_source, J);
			mg->detach_from_edges (a_in_source);
		}
	
	/// sets a zero value in the defect for all conductor indices
		void adjust_defect
		(
			pot_vector_type & d,
			const pot_vector_type & u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol = SPNULL,
			const std::vector<number> * vScaleMass = NULL,
			const std::vector<number> * vScaleStiff = NULL
		)
		{
			MultiGrid * mg = (MultiGrid *) dd->multi_grid().get ();
			a_vert_flag_type a_in_source;
			mg->attach_to_vertices (a_in_source);
			aa_vert_flag_type aa_in_source (*mg, a_in_source);
			mark_source_vertices (* dd.get (), aa_in_source);
			adjust_vector (* dd.get (), aa_in_source, d);
			mg->detach_from_edges (a_in_source);
		}
	
	/// sets the value in the solution for all conductor indices
		void adjust_solution
		(
			pot_vector_type	& u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		)
		{
			MultiGrid * mg = (MultiGrid *) dd->multi_grid().get ();
			a_vert_flag_type a_in_source;
			mg->attach_to_vertices (a_in_source);
			aa_vert_flag_type aa_in_source (*mg, a_in_source);
			mark_source_vertices (* dd.get (), aa_in_source);
			adjust_vector (* dd.get (), aa_in_source, u);
			mg->detach_from_edges (a_in_source);
		}
	
	///	sets unity rows in A and dirichlet values in right-hand side b
		void adjust_linear
		(
			pot_matrix_type & A,
			pot_vector_type & b,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		)
		{
			MultiGrid * mg = (MultiGrid *) dd->multi_grid().get ();
			a_vert_flag_type a_in_source;
			mg->attach_to_vertices (a_in_source);
			aa_vert_flag_type aa_in_source (*mg, a_in_source);
			mark_source_vertices (* dd.get (), aa_in_source);
			adjust_matrix (* dd.get (), aa_in_source, A);
			adjust_vector (* dd.get (), aa_in_source, b);
			mg->detach_from_edges (a_in_source);
		}
	
	///	sets the dirichlet value in the right-hand side
		void adjust_rhs
		(
			pot_vector_type & b,
			const pot_vector_type & u,
			ConstSmartPtr<DoFDistribution> dd,
			int type,
			number time = 0.0
		)
		{
			OutOfSource::adjust_solution (b, dd, type, time);
		}

	///	returns the type of the constraints
		int type () const {return CT_DIRICHLET;}
	
	private:
		
	///	the master object of the projection
		NedelecLoopCurrent & m_master;
	};
	
///	Computes the flux of of the gradient of the potential over the cut (for one type of elements)
	template <typename TElem>
	void get_flux_of_pot
	(
		const TDomain & domain, ///< [in] the domain
		const pot_vector_type & pot, ///< [in] the potential field
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		number & flux ///< [out] the flux to update
	);
	/// Helper class for computation of the flux of potential over the cut
	struct GetFluxOfPotential
	{
		this_type * m_pThis;
		const domain_type & m_domain;
		const pot_vector_type & m_pot;
		const DoFDistribution & m_vertDD;
		number & m_flux;
		
		GetFluxOfPotential
		(
			this_type * pThis, ///< [in] pointer to the master class
			const TDomain & domain, ///< [in] the domain
			const pot_vector_type & pot, ///< [in] the potential field
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			number & flux ///< [out] the flux to update
		)
		: m_pThis (pThis), m_domain (domain), m_pot (pot), m_vertDD (vertDD), m_flux (flux)
		{
			m_flux = 0;
		}
		
		template <typename TElem> void operator() (TElem &)
		{
			m_pThis->template get_flux_of_pot<TElem> (m_domain, m_pot, m_vertDD, m_flux);
		}
	};
	
private:
	
///	Structure for keeping electric current data
	struct TSrcData
	{
		std::string fctNames; ///< function names for the value
		number I; ///< value of the electric current
		
		/// Constructor
		TSrcData (const char * the_fctNames, number the_I)
			: fctNames (the_fctNames), I (the_I) {};
	};
	
///	Array of the electric current data
	std::vector<TSrcData> m_vSrcData;

///	Names of all the subsets of the source
	std::string m_allSsNames;
///	Names of the subsets of the positive direction
	std::string m_posSsNames;
///	Names of the surfaces of the cut of the loop
	std::string m_cutSsNames;
	
///	Group of all the subsets of the source
	SubsetGroup m_allSsGrp;
///	Group of the subsets of the positive direction
	SubsetGroup m_posSsGrp;
///	Group of the surfaces of the cut of the loop
	SubsetGroup m_cutSsGrp;
	
///	Vertex-centered approximation space
	SmartPtr<ApproximationSpace<TDomain> > m_spVertApproxSpace;
	
///	Local discretization of the auxiliary equations
	SmartPtr<AuxLaplaceLocAss> m_auxLocLaplace;
///	Extension of the matrices and vectors to the whole domain
	SmartPtr<OutOfSource> m_outOfSource;
///	Global discretization of the auxiliary equations
/**
 *	This discretization object uses m_auxLocLaplace and m_outOfSource
 *	to assemble the Laplace equation and the jump right-hand side in the
 *	subsets of the source, as well as the identity equations in the other
 *	parts of the domain.
 */
	SmartPtr<DomainDiscretization<TDomain, TPotAlgebra> > m_auxLaplaceAss;
///	Matrix of the discretization
	SmartPtr<AssembledLinearOperator<TPotAlgebra> > m_auxLaplaceOp;
///	Solver for the auxliliary equations
	SmartPtr<ILinearOperatorInverse<pot_vector_type> > m_potSolver;
	
///	Computed scaling of the potential
	number m_pot_scaling;
};

} // end namespace Electromagnetism
} // end namespace ug

#include "nedelec_source_impl.h"

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_PROJECT__

/* End of File */
