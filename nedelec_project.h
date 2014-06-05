/*
 * nedelec_project.h - classes for the projection of the functions based
 * on the Nedelec element to the space of the divergence-free functions.
 * 
 * Created on: 07.04.2013
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_PROJECT__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_PROJECT__

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

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

#include "em_material.h"

namespace ug{
namespace Electromagnetism{

/// Projection procedure class
/**
 * This class implements the procedure for projection of vector fields
 * in the Nedelec-element representation to the divergence-free space
 * and elimination of the non-zero potentials of conductors.
 */
template <typename TDomain, typename TAlgebra>
class NedelecProject
{
/// This type
	typedef NedelecProject<TDomain, TAlgebra> this_type;
	
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

public:

///	Constructor
	NedelecProject
	(
		SmartPtr<EMaterial<TDomain> > emMatherial, ///< properties of the materials
		SmartPtr<ApproximationSpace<TDomain> > vertApproxSpace, ///< vertex-centered approx. space
		SmartPtr<ILinearOperatorInverse<pot_vector_type> > potSolver ///< linear solver for the potential
	);
	
///	Destructor
	~NedelecProject ()
	{
		for (size_t i = 0; i < m_DVF_phi.size (); i++) delete m_DVF_phi [i];
	}
	
///	Sets the Dirichlet boundary
	void set_Dirichlet
	(
		SmartPtr<EMDirichlet<TDomain, TAlgebra> > spDirichlet ///< the Dirichlet BC object
	)
	{
		m_spDirichlet = spDirichlet;
	}
	
///	Sets whether to damp Dirichlet vector fields
	void set_dampDVFs (bool val)
	{
		m_bDampDVFs = val;
	}
	
///	Performs the projection
	void apply
	(
		SmartPtr<GridFunction<TDomain, TAlgebra> > sp_u, ///< [in] the grid function to project
		const char * fct_name ///< [in] the function name
	);
	
///	Computes the weak divergence in insulators
	void compute_div
	(
		SmartPtr<GridFunction<TDomain, TAlgebra> > sp_u, ///< [in] the vector field grid function
		const char * u_fct_name, ///< [in] the function name of the Nedelec DoFs
		SmartPtr<GridFunction<TDomain, TPotAlgebra> > sp_div ///< [out] the grid function for the divergence
	);
	
private:
	
	/// Allocates memory for the DVFs associated with the ungrounded conductors
	void alloc_DVFs
	(
		const SmartPtr<TDomain> & domain, ///< [in] the domain
		pot_gf_type & aux_rhs ///< [in] grid function for the auxiliary rhs
	);
	
	///	Computes the Dirichlet vector fields (DVFs)
	void compute_DVFs
	(
		pot_gf_type & aux_rhs ///< grid function for the auxiliary rhs
	);
	
	/// Computes the potential coefficients
	void compute_DVF_potential_coeffs
	(
		const SmartPtr<TDomain> & domain, ///< [in] the domain
		const SmartPtr<DoFDistribution> & vertDD ///< [in] the vertex DD
	);
	
	/// Projects one function (i.e. performs the main action):
	inline void project_func
	(
		const SmartPtr<TDomain> & domain, ///< [in] the domain
		const SmartPtr<DoFDistribution> & edgeDD, ///< [in] the edge DD
		vector_type & u, ///< [in] the vector where to project
		size_t fct, ///< [in] function in u to compute div for
		const SmartPtr<DoFDistribution> & vertDD, ///< [in] the vertex DD
		pot_gf_type & aux_rhs, ///< [in] a grid function for the rhs of the aux. problem
		pot_gf_type & aux_cor ///< [in] a grid function for the sol. of the aux. problem
	);
	
//	The weak divergence and gradient operators:
	
	/// Computes the weak divergence (for all elements of the same type)
	template <typename TElem>
	void weak_div_elem_type
	(
		const TDomain & domain, ///< [in] the domain
		const DoFDistribution & edgeDD, ///< [in] the edge DD
		const vector_type & u, ///< [in] the vector to compute div for
		size_t fct, ///< [in] function in u to compute div for
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		pot_vector_type & div ///< to update the weak divergence of u
	);
	
	/// Helper class for the computation the weak divergence
	struct WeakDiv
	{
		this_type * m_pThis;
		const TDomain & m_domain;
		const DoFDistribution & m_edgeDD;
		const vector_type & m_u;
		size_t m_fct;
		const DoFDistribution & m_vertDD;
		pot_vector_type & m_div;
		
		WeakDiv
		(
			this_type * pThis,
			const TDomain & domain, ///< [in] the domain
			const DoFDistribution & edgeDD, ///< [in] the edge DD
			const vector_type & u, ///< [in] the vector to compute div for
			size_t fct, ///< [in] function in u to compute div for
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			pot_vector_type & div ///< to update the weak divergence of u
		)
		:	m_pThis (pThis), m_domain (domain), m_edgeDD (edgeDD), m_u (u),
			m_fct (fct), m_vertDD (vertDD), m_div (div) {}
		
		template <typename TElem> void operator() (TElem &)
		{
			m_pThis->template weak_div_elem_type<TElem>
				(m_domain, m_edgeDD, m_u, m_fct, m_vertDD, m_div);
		}
	};
	
	/// Sets the div operator to 0 in conductors (for all elements of the same type)
	template <typename TElem>
	void clear_div_in_conductors
	(
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		pot_vector_type & div, ///< [out] to update the weak divergence of u
		DenseVector<VariableArray1<number> > & charge ///< [out] charges of the conductors
	);
	
	/// Helper class for assembling the weak divergence operator
	struct ClearDivInConductors
	{
		this_type * m_pThis;
		const DoFDistribution & m_vertDD;
		pot_vector_type & m_div;
		DenseVector<VariableArray1<number> > & m_charge;
		
		ClearDivInConductors
		(
			this_type * pThis,
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			pot_vector_type & div, ///< [out] to update the weak divergence of u
			DenseVector<VariableArray1<number> > & charge ///< [out] charges of the conductors
		)
		: m_pThis (pThis), m_vertDD (vertDD), m_div (div), m_charge (charge) {}
		
		template <typename TElem> void operator() (TElem &)
		{
			m_pThis->template clear_div_in_conductors<TElem> (m_vertDD, m_div, m_charge);
		}
	};
	
	/// Assembles the weak divergence operator on the whole grid
	void assemble_div
	(
		const TDomain & domain, ///< [in] the domain
		const DoFDistribution & edgeDD, ///< [in] the edge DD
		const vector_type & u, ///< [in] the vector to compute div for
		size_t fct, ///< [in] function in u to compute div for
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		pot_vector_type & div, ///< [out] for the weak divergence of u
		DenseVector<VariableArray1<number> > & charge ///< [out] charges of the conductors
	);
	
	/// Damps the Dirichlet vector fields (DVFs)
	void damp_DVFs
	(
		pot_vector_type & cor, ///< the potential correction to update
		const DenseVector<VariableArray1<number> > & charge ///< [in] charges of the conductors
	);
	
	/// Updates the grid function by the potential correction
	void distribute_cor
	(
		const TDomain & domain, ///< [in] the domain
		const DoFDistribution & edgeDD, ///< [in] the edge DD
		vector_type & u, ///< [out] the vector to correct (to update)
		size_t fct, ///< [in] function in u to update
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		const pot_vector_type & cor ///< [in] the potential correction for u
	);
	
//	Laplace operators for the auxiliary problems
	
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
			NedelecProject & master ///< [in] the master object
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
		void prepare_element (const LocalVector & u, GridObject * elem, const position_type vCornerCoords [])
			{};
		template <typename TElem>
		void finish_element_loop ()
			{};
		template <typename TElem>
		void ass_rhs_elem (LocalVector & d, GridObject * elem, const position_type vCornerCoords [])
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
		NedelecProject & m_master;
		
	///	flag whether to assemble on a given element
		bool m_do_assemble_here;
	};
	
///	Constraint that assembles the rhs und the bc for the auxiliary problems
	class AuxLaplaceRHS : public IDomainConstraint<TDomain, TPotAlgebra>
	{
	private:
	/// Iterator over edges
		typedef DoFDistribution::traits<Edge>::const_iterator t_edge_iterator;
	
	/// maximum number of corners of an element
		static const size_t maxCorners = (size_t) element_list_traits<typename domain_traits<WDim>::DimElemList>::maxCorners;
		
	/// Iterator over vertices
		typedef DoFDistribution::traits<Vertex>::const_iterator t_vert_iterator;
	
	public:
	///	constructor
		AuxLaplaceRHS
		(
			NedelecProject & master ///< [in] the master object
		);
		
	///	sets the base conductor index (or -1 for the divergence correction)
		void set_base_conductor
		(
			int base_cond = -1 ///< index of the base conductor
		)
		{
			m_base_cond = base_cond;
		}
		
	/// sets a unity row for all conductor indices
		void adjust_jacobian
		(
			pot_matrix_type & J,
			const pot_vector_type & u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol = SPNULL,
			const number s_a0 = 1.0
		);
	
	/// sets a zero value in the defect for all conductor indices
		void adjust_defect
		(
			pot_vector_type & d,
			const pot_vector_type & u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0,
			ConstSmartPtr<VectorTimeSeries<pot_vector_type> > vSol = SPNULL,
			const std::vector<number> * vScaleMass = NULL,
			const std::vector<number> * vScaleStiff = NULL
		);
	
	/// sets the value in the solution for all conductor indices
		void adjust_solution
		(
			pot_vector_type	& u,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		);
	
	///	sets unity rows in A and dirichlet values in right-hand side b
		void adjust_linear
		(
			pot_matrix_type & A,
			pot_vector_type & b,
			ConstSmartPtr<DoFDistribution> dd,
			number time = 0.0
		)
		{ // Note that this function is not really used, so it needs not to be optimal.
			AuxLaplaceRHS::adjust_jacobian (A, b, dd, time); // the 2nd arg. is dummy: A does not depend on u
			AuxLaplaceRHS::adjust_solution (b, dd, time);
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
			AuxLaplaceRHS::adjust_solution (b, dd, time);
		}

	///	returns the type of the constraints
		int type () const {return CT_DIRICHLET;}
	
	private:
		
	///	sets zero values on the whole Dirichlet boundary
		void set_zero_Dirichlet
		(
			pot_vector_type & u, ///< the vector where to set zeros
			const DoFDistribution * dd ///< the vert.-based DoF distribution
		);
	
	///	sets identity matrix on the whole Dirichlet boundary
		void set_identity_Dirichlet
		(
			pot_matrix_type & A, ///< the matrix to set
			const DoFDistribution * dd ///< the vert.-based DoF distribution
		);
	
	/// sets constant value on the closure of a full-dimensional subset
		template <typename TElem>
		void set_value_on_subset
		(
			int si, ///< the subset
			number val, ///< the value to set
			pot_vector_type & u, ///< the vector where to set
			const DoFDistribution * dd ///< the vert.-based DoF distribution
		);
		
	/// Helper class for 'set_value_on_subset'
		struct SetValueOnSubset
		{
			AuxLaplaceRHS * m_pThis;
			int m_si;
			number m_val;
			pot_vector_type & m_u;
			const DoFDistribution * m_dd;
			
			SetValueOnSubset
			(
				AuxLaplaceRHS * pThis,
				int si, ///< the subdomain where to set the values
				number val, ///< the value to set
				pot_vector_type & u, ///< the vector where to set the values
				const DoFDistribution * dd ///< the edge DD
			)
			: m_pThis (pThis), m_si (si), m_val (val), m_u (u), m_dd (dd) {};
			
			template <typename TElem> void operator() (TElem &)
			{
				m_pThis->template set_value_on_subset<TElem> (m_si, m_val, m_u, m_dd);
			}
		};
	
	/// sets identity matrix on the closure of a full-dimensional subset
		template <typename TElem>
		void set_identity_on_subset
		(
			int si, ///< the subset
			pot_matrix_type & A, ///< the matrix to set
			const DoFDistribution * dd ///< the vert.-based DoF distribution
		);
	
	/// Helper class for 'set_value_on_subset'
		struct SetIdentityOnSubset
		{
			AuxLaplaceRHS * m_pThis;
			int m_si;
			pot_matrix_type & m_A;
			const DoFDistribution * m_dd;
			
			SetIdentityOnSubset
			(
				AuxLaplaceRHS * pThis,
				int si, ///< the subdomain where to set the values
				pot_matrix_type & A, ///< the matrix to set
				const DoFDistribution * dd ///< the edge DD
			)
			: m_pThis (pThis), m_si (si), m_A (A), m_dd (dd) {};
			
			template <typename TElem> void operator() (TElem &)
			{
				m_pThis->template set_identity_on_subset<TElem> (m_si, m_A, m_dd);
			}
		};
	
	private:
		
	///	the master object of the projection
		NedelecProject & m_master;
	
	///	the base conductor index (or -1 for insulators)
		int m_base_cond;
	};
	
//	Computation of the charges of the DVFs
	
	/// Sets the base conductor indices for every vertex (for all elements of the same type)
	template <typename TElem>
	void mark_cond_vert_elem_type
	(
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		int * vert_base_cond ///< [out] indices of the base conductors for every vertex
	);
	
	/// Helper class for setting the base conductor indices to vertices
	struct MarkCondVert
	{
		this_type * m_pThis;
		const DoFDistribution & m_vertDD;
		int * m_vert_base_cond;
		
		MarkCondVert
		(
			this_type * pThis,
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			int * vert_base_cond ///< [out] indices of the base conductors for every vertex
		)
		:	m_pThis (pThis), m_vertDD (vertDD), m_vert_base_cond (vert_base_cond) {}
		
		template <typename TElem> void operator() (TElem &)
		{
			m_pThis->template mark_cond_vert_elem_type<TElem> (m_vertDD, m_vert_base_cond);
		}
	};
	
	/// Integrates div E over boundaries of conductors for elements of the same type
	template <typename TElem>
	void integrate_div_DVF_elem_type
	(
		const TDomain & domain, ///< [in] the domain
		const DoFDistribution & vertDD, ///< [in] the vertex DD
		const int * vert_base_cond, ///< [in] indices of the base conductors for every vertex
		DenseMatrix<VariableArray2<number> > & C ///< [out] the matrix to update
	);
	
	/// Helper class for computation of the charges of the DVFs
	struct IntegrateDivDVF
	{
		this_type * m_pThis;
		const TDomain & m_domain;
		const DoFDistribution & m_vertDD;
		const int * m_vert_base_cond;
		DenseMatrix<VariableArray2<number> > & m_C;
		
		IntegrateDivDVF
		(
			this_type * pThis,
			const TDomain & domain, ///< [in] the domain
			const DoFDistribution & vertDD, ///< [in] the vertex DD
			const int * vert_base_cond, ///< [out] indices of the base conductors for every vertex
			DenseMatrix<VariableArray2<number> > & C ///< [out] the matrix to update
		)
		:	m_pThis (pThis), m_domain (domain), m_vertDD (vertDD),
			m_vert_base_cond (vert_base_cond), m_C (C) {}
		
		template <typename TElem> void operator() (TElem &)
		{
			m_pThis->template integrate_div_DVF_elem_type<TElem>
				(m_domain, m_vertDD, m_vert_base_cond, m_C);
		}
	};
	
private:

///	Properties of the matherials in the domain
	SmartPtr<EMaterial<TDomain> > m_spEmMaterial;
	
/// Approximation space for the potential (vertex-centered) space
	SmartPtr<ApproximationSpace<TDomain> > m_spVertApproxSpace;
	
///	Whether to damp the Dirichlet vector fields (DVFs)
	bool m_bDampDVFs;
	
///	Dirichlet boundary
	SmartPtr<EMDirichlet<TDomain, TAlgebra> > m_spDirichlet;
	
///	Local discretization of the auxiliary equations
	SmartPtr<AuxLaplaceLocAss> m_auxLocLaplace;
///	Rhs and BC for the auxiliary problem
	SmartPtr<AuxLaplaceRHS> m_auxLaplaceRHS;
///	Global discretization of the auxiliary equations
	SmartPtr<DomainDiscretization<TDomain, TPotAlgebra> > m_auxLaplaceAss;
///	Matrix of the discretization
	SmartPtr<AssembledLinearOperator<TPotAlgebra> > m_auxLaplaceOp;
///	Solver for the auxliliary equations
	SmartPtr<ILinearOperatorInverse<pot_vector_type> > m_potSolver;
///	Dirichlet vector fields (DVFs)
	std::vector<pot_gf_type *> m_DVF_phi;
///	Potential coefficients:
	DenseMatrix<VariableArray2<number> > m_potCoeff;
};

} // end namespace Electromagnetism
} // end namespace ug

#include "nedelec_project_impl.h"

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_PROJECT__

/* End of File */
