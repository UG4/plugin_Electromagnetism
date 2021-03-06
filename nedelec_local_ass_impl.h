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

/*
 * nedelec_local_ass_impl.cpp
 * Implementation of the functions for assembling the local stiffness, mass and
 * weak divergence matrices of the Nedelec-type-1 (Whitney-1) based discretization
 * of the Maxwell equations, as well as the shape functions and their curls
 * of the Nedelec element.
 */

/* UG4 headers: */
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/reference_element/reference_mapping.h"

namespace ug{
namespace Electromagnetism{

/**
 * Computes the gradients of the Whitney-0 (Lagrange P1) shape functions.
 *
 * \remark We assume that the reference mapping is linear: Here, the Nedelec
 * elements are only implemented for simplices.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::compute_W0_grads
(
	const position_type * corners, /**< [in] array of the global corner coordinates */
	MathVector<WDim> grad [numCorners] /**< [out] gradients of the Whitney-0 functions */
)
{
// We assume here that TElem is a simplex (triangle or tetrahedron):
	UG_ASSERT ((ReferenceMapping<ref_elem_type, WDim>::isLinear),
		"Not a simplex (" << ref_elem_type::REFERENCE_OBJECT_ID
			<< ") as a template param. of NedelecT1_LDisc_forSimplex: Only triangles and tetrahedra accepted.");

// Whitney-0 shapes:
	const W0_shapes_type & W0_shapes = Provider<W0_shapes_type>::get ();
	
// Reference mapping
	ReferenceMapping<ref_elem_type, WDim> ref_mapping;
	ref_mapping.update (corners);

// Compute the gradients

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
 * assigned. We need them to be ordered independently of the element. This
 * function assignes to every edge (i.e. every edge dof index) two indices
 * of the corners of the element in order that corresponds to the global
 * order of the vertices associated with these corners.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::get_edge_corners
(
	const TDomain * domain, /**< [in] the domain */
	TElem * elem, /**< [in] the element */
	size_t edge_corner [numEdges] [2] /**< [out] edge dof -> corner of the element */
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	const ref_elem_type & rRefElem = Provider<ref_elem_type>::get ();
	
	const grid_type * grid = domain->grid().get ();
	UG_ASSERT ((grid != 0), "No grid in the domain.");
	
	for (size_t e = 0; e < numEdges; e++)
	{
		int co_0 = rRefElem.id (1, e, 0, 0);
		int co_1 = rRefElem.id (1, e, 0, 1);
		UG_ASSERT ((co_0 >= 0 && co_1 >= 0), "NedelecT1_LDisc_forSimplex::get_edge_corners: Internal error.");
		
		const Edge * edge = ((grid_type *) grid)->get_edge (elem, e);
		
		if (elem->vertex (co_0) == edge->vertex (0))
		{
			UG_ASSERT ((elem->vertex (co_0) == edge->vertex (0)), "NedelecT1_LDisc_forSimplex::get_edge_corners: Internal error.")
			edge_corner[e][0] = co_0;
			edge_corner[e][1] = co_1;
		}
		else
		{
			UG_ASSERT ((elem->vertex (co_0) == edge->vertex (1) && elem->vertex (co_1) == edge->vertex (0)),
				"NedelecT1_LDisc_forSimplex::get_edge_corners: Internal error.")
			edge_corner[e][0] = co_1;
			edge_corner[e][1] = co_0;
		}
	}
}

/**
 * Assembles the local stiffness and mass matrices of the rot-rot operator
 *
 * The matrices are not premultiplied by any physical factors:
 * these factors are assumed to be constant over the element and must
 * be applied later.
 *
 * \remark We assume that the reference mapping is linear: The Nedelec
 * elements are only implemented for simplices.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::local_stiffness_and_mass
(
	const TDomain * domain, /**< [in] the domain */
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
	
// REMARK: Below we compute the integrals of (w^{(1)}_e1,  w^{(1)}_e2) over the
// element by a generalization of the Simpson's quadrature rule for simplexes, cf.
// A. Horwitz, A version of Simpson’s rule for multiple integrals,
// Journal of Computational and Applied Mathematics 134 (2001), pp. 1-11,
// DOI: 10.1016/S0377-0427(00)00444-1

	static const number lambda = ((number) (WDim + 1)) / (WDim + 2);
	
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
// compute the values of w^{(1)}_e-functions at the corners
	MathVector<WDim> w1_at_co[numEdges][numCorners];
	for (size_t e = 0; e < numEdges; e++)
	{
	// Every w^{(1)} is non-zero at only two corners of the element: the ends
	// of the edge: At all the other corners the corresponding w^{(0)} are zero.
		for (size_t co = 0; co < numCorners; co++) w1_at_co[e][co] = 0;
		w1_at_co[e][edge_corner[e][0]] = grad_w0[edge_corner[e][1]];
		w1_at_co[e][edge_corner[e][1]] -= grad_w0[edge_corner[e][0]];
	}
// assemble the mass matrix
	for (size_t e_1 = 0; e_1 < numEdges; e_1++)
		for (size_t e_2 = 0; e_2 <= e_1; e_2++)
		{
			number t = 0;
			for (size_t co = 0; co < numCorners; co++)
				t += VecDot (w1_at_co[e_1][co], w1_at_co[e_2][co]);
			t /= numCorners;
			
			M[e_1][e_2] = (lambda * VecDot (w1_at_center[e_1], w1_at_center[e_2])
				+ (1 - lambda) * t) * V;
		}
	for (size_t e_1 = 0; e_1 < numEdges-1; e_1++)
		for (size_t e_2 = e_1 + 1; e_2 < numEdges; e_2++) M[e_1][e_2] = M[e_2][e_1];
}

/**
 * Assembles the local mass matrix of the Nedelec element (i.e. performs a part
 * of the task of local_stiffness_and_mass).
 *
 * \remark We assume that the reference mapping is linear: The Nedelec
 * elements are only implemented for simplices.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::local_mass
(
	const TDomain * domain, /**< [in] the domain */
	TElem * elem, /**< [in] element */
	const position_type * corners, /**< [in] array of the global corner coordinates */
	number M [maxNumEdges][maxNumEdges] /**< [out] local mass matrix */
)
{
// clear the local matrix:
	memset (M, 0, maxNumEdges * maxNumEdges * sizeof (number));
	
// get volume of the grid element
	number V = ElementSize<ref_elem_type, WDim> (corners);
	
/* Computation of the mass matrix: */
	
// get the gradients of the Whitney-0 elements
	MathVector<WDim> grad_w0 [numCorners];
	compute_W0_grads (corners, grad_w0);

// get the correspondence of the edges and the corners:
	size_t edge_corner [numEdges] [2];
	get_edge_corners (domain, elem, edge_corner);

// REMARK: Below we compute the integrals of (w^{(1)}_e1,  w^{(1)}_e2) over the
// element by a generalization of the Simpson's quadrature rule for simplexes, cf.
// A. Horwitz, A version of Simpson’s rule for multiple integrals,
// Journal of Computational and Applied Mathematics 134 (2001), pp. 1-11,
// DOI: 10.1016/S0377-0427(00)00444-1

	static const number lambda = ((number) (WDim + 1)) / (WDim + 2);
	
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
// compute the values of w^{(1)}_e-functions at the corners
	MathVector<WDim> w1_at_co[numEdges][numCorners];
	for (size_t e = 0; e < numEdges; e++)
	{
	// Every w^{(1)} is non-zero at only two corners of the element: the ends
	// of the edge: At all the other corners the corresponding w^{(0)} are zero.
		for (size_t co = 0; co < numCorners; co++) w1_at_co[e][co] = 0;
		w1_at_co[e][edge_corner[e][0]] = grad_w0[edge_corner[e][1]];
		w1_at_co[e][edge_corner[e][1]] -= grad_w0[edge_corner[e][0]];
	}
// assemble the mass matrix
	for (size_t e_1 = 0; e_1 < numEdges; e_1++)
		for (size_t e_2 = 0; e_2 <= e_1; e_2++)
		{
			number t = 0;
			for (size_t co = 0; co < numCorners; co++)
				t += VecDot (w1_at_co[e_1][co], w1_at_co[e_2][co]);
			t /= numCorners;
			
			M[e_1][e_2] = (lambda * VecDot (w1_at_center[e_1], w1_at_center[e_2])
				+ (1 - lambda) * t) * V;
		}
	for (size_t e_1 = 0; e_1 < numEdges-1; e_1++)
		for (size_t e_2 = e_1 + 1; e_2 < numEdges; e_2++) M[e_1][e_2] = M[e_2][e_1];
}

/**
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
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::local_div_matrix
(
	const TDomain * domain, /**< [in] the domain */
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

/**
 * Computes the Nedelec shapes at a given point
 * 
 * \remark If the dimensionality of the given element is lower than that of
 * the world then the shapes are projected to the element.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::get_shapes
(
	const TDomain * domain, /**< [in] the domain */
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

/**
 * Computes of the values of the grid functions.
 *
 * The values of the DoFs are scalar, but the resulting function is a vector field.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::interpolate
(
	const TDomain * domain, /**< [in] the domain */
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

/**
 * Computes the curl of the grid functions
 *
 * \remark Note that curl is constant over the element.
 * 
 * In 2d, the value of the curl should have the form \f$(0, 0, z)\f$. But
 * to keep the dimensionality of the output, it is represented as \f$(z, 0)\f$.
 */
template <typename TDomain, typename TElem>
void NedelecT1_LDisc_forSimplex<TDomain, TElem>::curl
(
	const TDomain * domain, /**< [in] the domain */
	TElem * elem, /**< [in] element */
	const position_type * corners, /**< [in] array of the global corner coordinates */
	const number dof [], /**< [in] arrays of values of the Nedelec degrees of freedom */
	MathVector<WDim> & curl_vec /**< [out] where to store the computed curl */
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
	curl_vec = 0.0;
	for (size_t e = 0; e < numEdges; e++)
		VecScaleAppend (curl_vec, dof[e], half_rot_w1[e]);
	curl_vec *= 2;
}

/**
 * This function computes the flux of the curl through a given plane inside
 * of a given grid element. The plane is identified by a point on it and the
 * normal. The flux is multiplied by the norm of the given normal vector (i.e.
 * specify the unit normal to get the standard flux). The function returns 
 * the area of the intersection if the element is intersected by the
 * plane. Otherwise the function returns exactly 0.0.
 */
template <typename TDomain>
number NedelecInterpolation<TDomain, 3, 3>::curl_flux
(
	const TDomain * domain, /**< [in] the domain */
	GridObject * elem, /**< [in] element */
	const position_type corners [], /**< [in] array of the global corner coordinates */
	const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
	const MathVector<3> & normal, /**< [in] normal to the plane */
	const position_type pnt, /**< [in] point on the plane (identifying the plane) */
	number & flux /**< [out] the flux */
)
{
//	current implementation considers only tetrahedra
	if (elem->reference_object_id () != ROID_TETRAHEDRON)
		UG_THROW ("NedelecInterpolation::curl_flux:"
					" No implementation of the Nedelec elements for"
					" Reference Object " << elem->reference_object_id () <<
					" of reference dim. 3 in a 3d domain. (This must be a tetrahedron.)");
	
//	compute the flux
	number area;
	MathVector<3> sect_corners [4];
	size_t n_intersect = IntersectPlaneWithTetrahedron (sect_corners, pnt, normal, corners);
	switch (n_intersect)
	{
	case 0:
		flux = 0.0;
		return 0.0;
		
	case 3:
		area = TriangleArea (sect_corners[0], sect_corners[1], sect_corners[2]);
		break;
	case 4:
		area = TriangleArea (sect_corners[0], sect_corners[1], sect_corners[2])
				+ TriangleArea (sect_corners[0], sect_corners[2], sect_corners[3]);
		break;
	
	default:
		UG_THROW ("NedelecInterpolation::curl_flux:"
			" Illegal number of intersections of a plane with a tetrahedron.");
	}
	
//	compute the curl (as a vector); note that it is constant in the element
	MathVector<3> curl_vec;
	NedelecT1_LDisc<TDomain, Tetrahedron>::curl (domain, (Tetrahedron *) elem,
		corners, dofs, curl_vec);
	flux = VecDot (curl_vec, normal) * area;
	
	return area;
};

} // end namespace Electromagnetism
} // end namespace ug

/* End of File */
