/*
 * nedelec_local_ass_impl.cpp
 * Implementation of the functions for assembling the local stiffness, mass and
 * weak divergence matrices of the Nedelec-type-1 (Whitney-1) based discretization
 * of the Maxwell equations, as well as the shape functions and their curls
 * of the Nedelec element.
 *
 * Created on: 18.09.2012
 * Author: D. Logashenko
 */

/* UG4 headers: */
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/reference_element/reference_mapping.h"

namespace ug{
namespace Electromagnetism{

/** computes the gradients of the Whitney-0 (Lagrange P1) shape functions
 *
 * Note that we assume that the reference mapping is linear: The Nedelec
 * elements are only implemented for simplices.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc<TDomain, TElem>::compute_W0_grads
(
	const position_type * corners, /**< [in] array of the global corner coordinates */
	MathVector<WDim> grad [numCorners] /**< [out] gradients of the Whitney-0 functions */
)
{
// We work with simplices only:
	if (! ReferenceMapping<ref_elem_type, WDim>::isLinear)
		UG_THROW ("The Nedelec elements are implemented for triangles and tetrahedra only.");

// Whitney-0 shapes:
	const W0_shapes_type & W0_shapes = Provider<W0_shapes_type>::get ();
	
///	reference mapping
	ReferenceMapping<ref_elem_type, WDim> ref_mapping;
	ref_mapping.update (corners);

/// compute the gradients

	MathVector<dim> local_pos;
	MathVector<dim> local_grad [numCorners];
	MathMatrix<WDim,dim> JtInv;
	
	local_pos = 0; // we compute the gradients at this point (remember: the mapping is linear)
	W0_shapes.grads (local_grad, local_pos);
	ref_mapping.jacobian_transposed_inverse (JtInv, local_pos);
	for (size_t co = 0; co < numCorners; co++)
		MatVecMult (grad[co], JtInv, local_grad[co]);
}

/** gets the correspondence between the edge dof indices and the corners of the element
 *
 * The edge dof indices correspond to the order of the edges in the reference
 * element. To every edge, two vertices (the corners of the element) are
 * assigned. We need them to be ordered independently on the element. This
 * function assignes to every edge (i.e. every edge dof index) two indices
 * of the corners of the element in order that corresponds to the global
 * order of the vertices associated with these corners.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc<TDomain, TElem>::get_edge_corners
(
	const TDomain& domain, /** [in] the domain */
	TElem * elem, /**< [in] the element */
	size_t edge_corner [numEdges] [2] /**< [out] edge dof -> corner of the element */
)
{
	const grid_type * grid = domain.grid().get ();
	Grid::edge_traits::secure_container edge_list;
	
	UG_ASSERT ((grid != 0), "No grid in the domain.");
	((grid_type *) grid)->associated_elements_sorted (edge_list, elem);
	UG_ASSERT ((edge_list.size () == numEdges), "Mismatch of numbers of corners and vertices of an element");
	
	for (size_t e = 0; e < numEdges; e++)
	{
		const EdgeBase * edge = edge_list[e];
		for (size_t i = 0; i < 2; i++)
		{
			const VertexBase * vert = edge->vertex (i); size_t co = 0;
			while (elem->vertex (co) != vert)
				if ((++co) == numCorners)
					UG_THROW ("Internal error in Nedelec disc.: vertex-edge mismatch");
			edge_corner[e][i] = co;
		}
	}
}

/** assembles the stiffness and mass matrices of the rot-rot operator
 *
 * The matrices are not premultiplied by any physical factors:
 * these factors are assumed to be constant over the element and must
 * be applied later.
 * Note that we assume that the reference mapping is linear: The Nedelec
 * elements are only implemented for simplices.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc<TDomain, TElem>::local_stiffness_and_mass
(
	const TDomain & domain, /** [in] the domain */
	TElem * elem, /**< [in] element */
	const position_type * corners, /**< [in] array of the global corner coordinates */
	number S [maxNumEdges][maxNumEdges], /**< [out] local stiffness matrix */
	number M [maxNumEdges][maxNumEdges] /**< [out] local mass matrix */
)
{
// clear the local matrices:
	memset (S, 0, maxNumEdges * maxNumEdges * sizeof (number));
	memset (M, 0, maxNumEdges * maxNumEdges * sizeof (number));
	
// get volume of the grid element
	number V = ElementSize<ref_elem_type, WDim> (corners);
	
/* 1. Computation of the stiffness matrix: */
	
// get the gradients of the Whitney-0 elements
	MathVector<WDim> grad_w0 [numCorners];
	compute_W0_grads (corners, grad_w0);

// get the correspondence of the edges and the corners:
	size_t edge_corner [numEdges] [2];
	get_edge_corners (domain, elem, edge_corner);

// compute \f$ \mathbf{rot} w^{(1)}_e / 2 = \mathbf{grad} w^{(0)}_m \times \mathbf{grad} w^{(0)}_n \f$
// for every edge \f$ e = (m, n) \f$
	MathVector<WDim> half_rot_w1 [numEdges];
	for (size_t e = 0; e < numEdges; e++)
		GenVecCross (half_rot_w1[e], grad_w0[edge_corner[e][0]], grad_w0[edge_corner[e][1]]);

// compute the entries of the local stiffness matrix using its symmetry
	for (size_t e_1 = 0; e_1 < numEdges; e_1++)
		for (size_t e_2 = 0; e_2 <= e_1; e_2++)
			S[e_1][e_2] = VecDot (half_rot_w1[e_1], half_rot_w1[e_2]) * 4 * V;
	for (size_t e_1 = 0; e_1 < numEdges-1; e_1++)
		for (size_t e_2 = e_1 + 1; e_2 < numEdges; e_2++) S[e_1][e_2] = S[e_2][e_1];
	
/* 2. Computation of the mass matrix: */
	
// REMARK: Below we compute the integrals of (w^{(1)}_e1,  w^{(1)}_e2) over the element
// by the MIDPOINT rule. The approximation error should be consistent with the one
// of the whole discretization. (Is it true?)

// compute the values of the w^{(1)}_e-functions at the center of the element
	MathVector<WDim> w1_at_center[numEdges];
	for (size_t e = 0; e < numEdges; e++)
	{
	// All the w^{(0)} (i.e. Lagrange) functions have the same value at the center
	// of the element. (This can be considered as a definition of the 'center'.)
	// This value is 1.0 / numCorners:
		VecSubtract (w1_at_center[e], grad_w0[edge_corner[e][0]], grad_w0[edge_corner[e][1]]);
		w1_at_center[e] /= numCorners;
	}
// assemble the mass matrix
	for (size_t e_1 = 0; e_1 < numEdges; e_1++)
		for (size_t e_2 = 0; e_2 <= e_1; e_2++)
			M[e_1][e_2] = VecDot (w1_at_center[e_1], w1_at_center[e_2]) * V;
	for (size_t e_1 = 0; e_1 < numEdges-1; e_1++)
		for (size_t e_2 = e_1 + 1; e_2 < numEdges; e_2++) M[e_1][e_2] = M[e_2][e_1];
}

/** assembles the discrete weak div operator
 *
 * This function assembles the local weak divergence operator \f$ \mathbf{B} \f$:
 * \f[
 *  \mathbf{B}_{v,e} = - \int_T w^{(1)}_e \cdot \mathbf{grad} w^{(0)}_v \, dx,
 * \f]
 * where \f$T\f$ is the element, \f$v\f$ its corner and \f$e\f$ its edge. The
 * integration is perfomed by the midpoint rule which is exact because the
 * integrand is linear. Thus, for \f$e = [v_1, v_2]\f$
 * \f[
 *  \mathbf{B}_{v,e} = - \frac{|T|}{N} \mathbf{grad} w^{(0)}_v
 *		\cdot (\mathbf{grad} w^{(0)}_{v_2} - \mathbf{grad} w^{(0)}_{v_1}) \, dx.
 * \f]
 * (with \f$N\f$ being the number of corners of \f$T\f$).
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc<TDomain, TElem>::local_div_matrix
(
	const TDomain & domain, /**< [in] the domain */
	TElem * elem, /**< [in] element */
	const position_type * corners, /**< [in] array of the global corner coordinates */
	number B [][numEdges] /**< [out] local weak divergence operator matrix */
)
{
// clear the local matrix:
	memset (B, 0, numCorners * numEdges * sizeof (number));
	
// get the gradients of the Whitney-0 elements
	MathVector<WDim> grad_w0 [numCorners];
	compute_W0_grads (corners, grad_w0);

// get the correspondence of the edges and the corners:
	size_t edge_corner [numEdges] [2];
	get_edge_corners (domain, elem, edge_corner);

// compute the values of the (- |T| w^{(1)}_e)-functions at the center of the element
	number V = - ElementSize<ref_elem_type, WDim> (corners) / numCorners;
	MathVector<WDim> T_w1_at_center[numEdges];
	for (size_t e = 0; e < numEdges; e++)
	{
		VecSubtract (T_w1_at_center[e], grad_w0[edge_corner[e][0]], grad_w0[edge_corner[e][1]]);
		T_w1_at_center[e] *= V;
	}
	
// compute the weak div matrix
	for (size_t v = 0; v < numCorners; v++)
		for (size_t e = 0; e < numEdges; e++)
			B[v][e] = VecDot (T_w1_at_center[e], grad_w0[v]);
}

/** computes the Nedelec shapes at a given point
 * 
 * Note that if the dimensionality of the given element is lower than that of
 * the world then the shapes are projected to the element.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc<TDomain, TElem>::get_shapes
(
	const TDomain & domain, /**< [in] the domain */
	TElem * elem, /**< [in] element */
	const position_type * corners, /**< [in] array of the global corner coordinates */
	const MathVector<dim> local, /**< [in] local coordinates of the point where to compute */
	MathVector<WDim> shapes [] /**< [out] array for the shapes */
)
{
// get the gradients of the Whitney-0 elements
	MathVector<WDim> grad_w0 [numCorners];
	compute_W0_grads (corners, grad_w0);

// get the correspondence of the edges and the corners:
	size_t edge_corner [numEdges] [2];
	get_edge_corners (domain, elem, edge_corner);

// Whitney-0 shapes:
	const W0_shapes_type& W0_shapes = Provider<W0_shapes_type>::get ();
	
// compute the shapes:
	number w0 [numCorners];
	W0_shapes.shapes (w0, local);
	for (size_t e = 0; e < numEdges; e++)
	{
		size_t m = edge_corner[e][0], n = edge_corner[e][1];
		VecScaleAdd (shapes[e], w0[m], grad_w0[n], - w0[n], grad_w0[m]);
	}
}

/** computes of the values of the grid functions
 *
 * The values of the DoFs are scalar, but the resulting function is a vector field.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc<TDomain, TElem>::interpolate
(
	const TDomain & domain, /** [in] the domain */
	TElem * elem, /**< [in] element */
	const position_type * corners, /**< [in] array of the global corner coordinates */
	const number dof [], /**< [in] arrays of values of the Nedelec degrees of freedom */
	const MathVector<dim> local [], /**< [in] local coordinates of the points where to compute */
	const size_t n_pnt, /**< [in] number of the points where to compute */
	MathVector<WDim> values [] /**< [out] where to store the computed n_pnt values */
)
{
// get the gradients of the Whitney-0 elements
	MathVector<WDim> grad_w0 [numCorners];
	compute_W0_grads (corners, grad_w0);

// get the correspondence of the edges and the corners:
	size_t edge_corner [numEdges] [2];
	get_edge_corners (domain, elem, edge_corner);

// Whitney-0 shapes:
	const W0_shapes_type & W0_shapes = Provider<W0_shapes_type>::get ();
	
// compute the values at all the points
	for (size_t pnt = 0; pnt < n_pnt; pnt++)
	{
		// get Whitney-0 shapes
		number w0 [numCorners];
		W0_shapes.shapes (w0, local[pnt]);
		
		// loop the shape functions
		MathVector<WDim>& value = values[pnt];
		value = 0.0;
		for (size_t e = 0; e < numEdges; e++)
		{
			MathVector<WDim> shape;
			size_t m = edge_corner[e][0], n = edge_corner[e][1];
			VecScaleAdd (shape, w0[m], grad_w0[n], - w0[n], grad_w0[m]);
			VecScaleAppend (value, dof[e], shape);
		}
	}
}

/** computes the curl of the grid functions
 *
 * Note that curl is constant over the element.
 * 
 * In 2d, the value of the curl should have the form \f$(0, 0, z)\f$. But
 * to keep the dimensionality of the output, it is represented as \f$(z, 0)\f$.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc<TDomain, TElem>::curl
(
	const TDomain & domain, /** [in] the domain */
	TElem * elem, /**< [in] element */
	const position_type * corners, /**< [in] array of the global corner coordinates */
	const number dof [], /**< [in] arrays of values of the Nedelec degrees of freedom */
	MathVector<WDim> & curl /**< [out] where to store the computed curl */
)
{
// get the gradients of the Whitney-0 elements
	MathVector<WDim> grad_w0 [numCorners];
	compute_W0_grads (corners, grad_w0);

// get the correspondence of the edges and the corners:
	size_t edge_corner [numEdges] [2];
	get_edge_corners (domain, elem, edge_corner);

// compute \f$ \mathbf{rot} w^{(1)}_e / 2 = \mathbf{grad} w^{(0)}_m \times \mathbf{grad} w^{(0)}_n \f$
// for every edge \f$ e = (m, n) \f$
	MathVector<WDim> half_rot_w1 [numEdges];
	for (size_t e = 0; e < numEdges; e++)
		GenVecCross (half_rot_w1[e], grad_w0[edge_corner[e][0]], grad_w0[edge_corner[e][1]]);
	
// compute the curl
	curl = 0.0;
	for (size_t e = 0; e < numEdges; e++)
		VecScaleAppend (curl, dof[e], half_rot_w1[e]);
	curl *= 2;
}

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
