/**
 * nedelec_source.h - classes for computation of divergence-free sources
 * in form of Nedelec-type-1 element DoFs.
 *
 * Created on: Jun. 30, 2013
 * Author: D. Logashenko
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

/**
 * Class for computation of loop currents
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

/// world dimention
	static const int WDim = TDomain::dim;
	
///	position type
	typedef typename TDomain::position_type position_type;

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

///	Computation of the source
	void compute
	(
		SmartPtr<GridFunction<TDomain, TAlgebra> > sp_u, ///< the grid function for the source
		const char * fct_names ///< the function name
	);
	
private:
	
	/// Marks edges that belong to the loop source domain (for one type of elements)
	template <typename TElem>
	void mark_source_edges
	(
		const DoFDistribution & edgeDD, ///< [in] the edge DD
		VariableArray1<char> & in_source ///< [out] the arrays of flags to update
	);
	/// Helper class for marking the edges in the source
	struct MarkSourceEdges
	{
		this_type * m_pThis;
		const DoFDistribution & m_edgeDD;
		VariableArray1<char> & m_in_source;
		
		MarkSourceEdges
		(
			this_type * pThis, ///< [in] pointer to the master class
			const DoFDistribution & edgeDD, ///< [in] the edge DD
			VariableArray1<char> & in_source ///< [out] the arrays of flags to update
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
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		const GridLevel & grid_lev, ///< [in] grid level of the vector to correct
		pot_vector_type & pot_u ///< [out] a vector for the potential
	);
	
	/// Computes the gradient of the potential
	void distribute_source_potential
	(
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		pot_vector_type & src_pot, ///< [in] potential to distribute
		const DoFDistribution & edgeDD, ///< [in] the edge DD
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
			GeometricObject * elem, ///< [in] element to compute for
			const position_type vCornerCoords [], ///< [in] coordinates of the corners of the element
			number loc_A [numCorners] [numCorners] ///< [out] the local stiffness matrix
		);
		
		/// Checks whether corners are on boundary (returns true if yes)
		static bool bnd_corners
		(
			GeometricObject * elem, ///< [in] element to compute for
			const SubsetGroup & bndGrp, ///< [in] boundary subsets
			bool bnd_flag [numCorners] ///< [out] true for bnd nodes, false for others
		);
	};
	
	/// Class for local assembling of the auxiliary Laplace operator
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
		void ass_JA_elem (LocalMatrix & J, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords []);
		
		template <typename TElem>
		void ass_rhs_elem (LocalVector & d, GeometricObject * elem, const position_type vCornerCoords []);

		template <typename TElem>
		void prepare_element (const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords [])
			{};
		template <typename TElem>
		void finish_element_loop ()
			{};
		template <typename TElem>
		void ass_JM_elem (LocalMatrix & J, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords [])
			{
				UG_THROW ("NedelecProject: Attempt to use nonstationary discretization for the potential.");
			};
		template <typename TElem>
		void ass_dA_elem (LocalVector & d, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords [])
			{
				UG_THROW ("NedelecProject: Attempt to assemble non-linear system for the potential.");
			};
		template <typename TElem>
		void ass_dM_elem (LocalVector & d, const LocalVector & u, GeometricObject * elem, const position_type vCornerCoords [])
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
	class OutOfSource : public IDomainConstraint<TDomain, TPotAlgebra>
	{
	private:
	/// maximum number of corners of an element
		static const size_t maxCorners = DimReferenceElement<WDim>::MAXCORNERS;
		
	/// Iterator over vertices
		typedef DoFDistribution::traits<Vertex>::const_iterator t_vert_iterator;
	
	private:
	
	/// marks the vertices belonging to elements in the source (for one element type)
		template <typename TElem>
		void mark_source_vertices_elem_type
		(
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			VariableArray1<bool> & in_source ///< [out] the arrays of flags to update
		);
	/// Helper class for marking the vertices in the source
		struct MarkSourceVertices
		{
			OutOfSource * m_pThis;
			const DoFDistribution & m_vertDD;
			VariableArray1<bool> & m_in_source;
			
			MarkSourceVertices
			(
				OutOfSource * pThis, ///< [in] pointer to the master class
				const DoFDistribution & vertDD, ///< [in] the vertex DD
				VariableArray1<bool> & in_source ///< [out] the arrays of flags to update
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
			VariableArray1<bool> & in_source ///< [out] the arrays of flags
		);
	
	///	sets to identity all the matrix rows that do not belong to the closure of the source domain
		void adjust_matrix
		(
			VariableArray1<bool> & in_source, ///< the arrays of flags
			pot_matrix_type & A ///< the matrix to adjust
		);
		
	///	sets to 0 all the entries of a vector that do not belong to the closure of the source domain
		void adjust_vector
		(
			VariableArray1<bool> & in_source, ///< the arrays of flags
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
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol = NULL,
			const number s_a0 = 1.0
		)
		{
			VariableArray1<bool> in_source (dd->num_indices ());
			mark_source_vertices (* dd.get (), in_source);
			adjust_matrix (in_source, J);
		}
	
	/// sets a zero value in the defect for all conductor indices
		void adjust_defect
		(
			pot_vector_type & d,
			const pot_vector_type & u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol = NULL,
			const std::vector<number> * vScaleMass = NULL,
			const std::vector<number> * vScaleStiff = NULL
		)
		{
			VariableArray1<bool> in_source (dd->num_indices ());
			mark_source_vertices (* dd.get (), in_source);
			adjust_vector (in_source, d);
		}
	
	/// sets the value in the solution for all conductor indices
		void adjust_solution
		(
			pot_vector_type	& u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		)
		{
			VariableArray1<bool> in_source (dd->num_indices ());
			mark_source_vertices (* dd.get (), in_source);
			adjust_vector (in_source, u);
		}
	
	///	sets unity rows in A and dirichlet values in right-hand side b
		void adjust_linear
		(
			pot_matrix_type & A,
			pot_vector_type & b,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		)
		{
			VariableArray1<bool> in_source (dd->num_indices ());
			mark_source_vertices (* dd.get (), in_source);
			adjust_matrix (in_source, A);
			adjust_vector (in_source, b);
		}
	
	///	sets the dirichlet value in the right-hand side
		void adjust_rhs
		(
			pot_vector_type & b,
			const pot_vector_type & u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		)
		{
			OutOfSource::adjust_solution (b, dd, time);
		}

	///	returns the type of the constraints
		int type () const {return CT_DIRICHLET;}
	
	private:
		
	///	the master object of the projection
		NedelecLoopCurrent & m_master;
	};
	
private:

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
	
	SmartPtr<ApproximationSpace<TDomain> > m_spVertApproxSpace;
	
///	Local discretization of the auxiliary equations
	SmartPtr<AuxLaplaceLocAss> m_auxLocLaplace;
///	Extension of the matrices and vectors to the whole domain
	SmartPtr<OutOfSource> m_outOfSource;
///	Global discretization of the auxiliary equations
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
