/*
 * nedelec_local_ass.h
 * Declarations of the functions for assembling the local stiffness, mass and
 * weak divergence matrices of the Nedelec-type-1 (Whitney-1) based discretization
 * of the Maxwell equations, as well as the shape functions and their curls
 * of the Nedelec element.
 *
 * Created on: 18.09.2012
 * Author: D. Logashenko
 */

#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__ROT_ROT_ASS__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__ROT_ROT_ASS__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/element_list_traits.h"
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"

namespace ug{
namespace Electromagnetism{

/// \ingroup lib_disc_elem_disc
/// @{

/// Tool kit for the Whitney-1 (Nedelec) based FE discretization of the rot-rot operators
/**
 * Class for the local discretization of the rot-rot operator using
 * the Nedelec-type-1 (Whitney-1) elements. This class may not contain
 * any member except for static functions, so no instances of this class
 * need to be created.
 *
 * This basic class template does not contain any numerical algorithms as they
 * are specific for different types of the elements. Cf. the specializations
 * for simplices (\see NedelecT1_LDisc_forSimplex).
 *
 * \remark
 * Note that the definition of the Whitney-1 function depends on the ordering of
 * the nodes in the element. More precisely, for a different ordering, the shape
 * function \f$w^{(1)}_e\f$ may have a different sign. Thus, the ordering should
 * be consistent in all the elements: If two elements \f$T_1\f$ and \f$T_2\f$
 * share the edge \f$e\f$ with the ends at corners \f$i\f$ and \f$j\f$
 * (that are shared, too) then this edge should be \f$e = (i, j)\f$ for both
 * the elements (and not for example \f$e = (i, j)\f$ for \f$T_1\f$ but
 * \f$e = (j, i)\f$ for \f$T_2\f$). In all the other respects, the ordering
 * does not play any essential role. In particular, the singes of the degrees of
 * freedom in the solution depend on the global ordering of the corners (but
 * the absolute values of these dofs are invariant). For this, the solution
 * of the discretized system may be only considered in connection with the
 * given ordering of the corners (in the edges).
 */
template <typename TDomain, typename TElem>
class NedelecT1_LDisc
{
public:

/// world dimention
	static const int WDim = TDomain::dim;
	
/// type of the grid
	typedef typename TDomain::grid_type grid_type;

/// type of the geometric positions (WDim-vectors)
	typedef typename TDomain::position_type position_type;
	
///	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

/// shapes and derivatives of the Whitney-0 (Lagrange P1) shape functions
	typedef LagrangeP1<ref_elem_type> W0_shapes_type;

///	dimension of reference element
	static const int dim = ref_elem_type::dim;

/// total number of the corners
	static const size_t numCorners = ref_elem_type::numCorners;
	
/// total number of the edges
	static const size_t numEdges = ref_elem_type::numEdges;
	
/// max. number of the edges of the full-dimensional elements in the domain
	static const size_t maxNumEdges = (size_t) element_list_traits<typename domain_traits<WDim>::DimElemList>::maxEdges;

public:

/// assembles the stiffness matrix of the rot-rot operator
	static void local_stiffness_and_mass
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		number S [maxNumEdges][maxNumEdges], /**< [out] local stiffness matrix */
		number M [maxNumEdges][maxNumEdges] /**< [out] local mass matrix */
	)
	{
		UG_THROW ("Whitney-1 (Nedelec) shapes not implemented for roid " << ref_elem_type::REFERENCE_OBJECT_ID << ".");
	};
	
///	assembles the discrete weak div operator
	static void local_div_matrix
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		number B [][numEdges] /**< [out] local weak divergence operator matrix */
	)
	{
		UG_THROW ("Whitney-1 (Nedelec) shapes not implemented for roid " << ref_elem_type::REFERENCE_OBJECT_ID << ".");
	};
	
///	computes the Nedelec shapes at a given point
	static void get_shapes
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		const MathVector<dim> local, /**< [in] local coordinates of the point where to compute */
		MathVector<WDim> shapes [] /**< [out] array for the shapes */
	)
	{
		UG_THROW ("Whitney-1 (Nedelec) shapes not implemented for roid " << ref_elem_type::REFERENCE_OBJECT_ID << ".");
	};
	
///	computes of the values of the grid functions
	static void interpolate
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<dim> local [], /**< [in] local coordinates of the points where to compute */
		const size_t n_pnt, /**< [in] number of the points where to compute */
		MathVector<WDim> values [] /**< [out] where to store the computed n_pnt values */
	)
	{
		UG_THROW ("Whitney-1 (Nedelec) shapes not implemented for roid " << ref_elem_type::REFERENCE_OBJECT_ID << ".");
	};
	
/// computes the curl of the grid functions (in 2d represented as \f$(z, 0)\f$ instead of \f$(0, 0, z)\f$)
	static void curl
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		MathVector<WDim> & curl_vec /**< [out] where to store the computed curl */
	)
	{
		UG_THROW ("Whitney-1 (Nedelec) shapes not implemented for roid " << ref_elem_type::REFERENCE_OBJECT_ID << ".");
	};
};

/// Helper class for the specialization of NedelecT1_LDisc for simplices (triangles and tetrahedrons)
/**
 * This class template implements the Whitney-1 (Nedelec) elements for triangles
 * and tetrahedrons (i.e. implements the classical variant of these elements),
 * \see NedelecT1_LDisc.
 *
 * Consider a grid element \f$T\f$. Let \f$w^{(0)}_i\f$ denote the linear (Lagrangian)
 * shape functions equal to 1 at corner i and to 0 at the other corners of \f$T\f$.
 * These are the Whitney-0 shape functions. The Whitney-1-shape functions \f$w^{(1)}_e\f$
 * (constituting the Nedelec element) are assigned to edges \f$e\f$ of \f$T\f$.
 * For the edge \f$e = (i, j)\f$, this function is given by formula
 * \f[
 *  w^{(1)}_e = w^{(0)}_i \mathbf{grad} w^{(0)}_j - w^{(0)}_j \mathbf{grad} w^{(0)}_i.
 * \f]
 * These are vector-valued shape functions. They are well-defined for simplices
 * (i.e. for triangles in 2d and tetrahedra in 3d), whereas the corresponding
 * physical problems are typically formulated only in 3d.
 *
 * Functions \f$w^{(1)}_e\f$ are used for discretization of the differential
 * operators of type \f$\mathbf{rot} \mu^{-1} \mathbf{rot}\f$ with the discretization
 * error \f$ O(h) \f$. The physical parameters (like \f$ \mu \f$) are then
 * considered constant over \f$T\f$, so that this tool kit sets the physical
 * factors to 1. (The local matrices can be multiplied by these factors later.)
 * The local stiffness matrix computed by this tool kit is
 * \f[
 *  (S_{h,T})_{e_1,e_2} = \int_T \mathbf{rot} w^{(1)}_{e_1} \cdot \mathbf{rot} w^{(1)}_{e_2} dx,
 * \f]
 * and the mass matrix is
 * \f[
 *  (M^{(1)}_{h,T})_{e_1,e_2} = \int \sigma w^{(1)}_{e_1} \cdot w^{(1)}_{e_2} dx.
 * \f]
 * The computation of the stiffness matrix is exact. (The integrand is a linear function.)
 * But the computation of the mass matrix is implemented by an inexact quadrature formula.
 *
 * \remark
 * Note that the definition of the Whitney-1 function depends on the ordering of
 * the nodes in the element. More precisely, for a different ordering, the shape
 * function \f$w^{(1)}_e\f$ may have a different sign. Thus, the ordering should
 * be consistent in all the elements: If two elements \f$T_1\f$ and \f$T_2\f$
 * share the edge \f$e\f$ with the ends at corners \f$i\f$ and \f$j\f$
 * (that are shared, too) then this edge should be \f$e = (i, j)\f$ for both
 * the elements (and not for example \f$e = (i, j)\f$ for \f$T_1\f$ but
 * \f$e = (j, i)\f$ for \f$T_2\f$). In all the other respects, the ordering
 * does not play any essential role. In particular, the singes of the degrees of
 * freedom in the solution depend on the global ordering of the corners (but
 * the absolute values of these dofs are invariant). For this, the solution
 * of the discretized system may be only considered in connection with the
 * given ordering of the corners (in the edges).
 *
 * References:
 * <ul>
 * <li> O. Sterz. Modellierung und Numerik zeitharmonischer Wirbelstromprobleme in 3D. PhD thesis, 2003.
 * <li> A. Bossavit. Computational Electromagnetism. Academic Press (Boston), 1998 (available in the Internet)
 * </ul>
 *
 * \tparam	TDomain		Domain type
 * \tparam	TElem		Element type
 */
template <typename TDomain, typename TElem>
class NedelecT1_LDisc_forSimplex
{
public:

/// world dimention
	static const int WDim = TDomain::dim;
	
/// type of the grid
	typedef typename TDomain::grid_type grid_type;

/// type of the geometric positions (WDim-vectors)
	typedef typename TDomain::position_type position_type;
	
///	type of reference element
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

/// shapes and derivatives of the Whitney-0 (Lagrange P1) shape functions
	typedef LagrangeP1<ref_elem_type> W0_shapes_type;

///	dimension of reference element
	static const int dim = ref_elem_type::dim;

/// total number of the corners
	static const size_t numCorners = ref_elem_type::numCorners;
	
/// total number of the edges
	static const size_t numEdges = ref_elem_type::numEdges;
	
/// max. number of the edges of the full-dimensional elements in the domain
	static const size_t maxNumEdges = (size_t) element_list_traits<typename domain_traits<WDim>::DimElemList>::maxEdges;

private:
	
/// computes the gradients of the Whitney-0 (Lagrange P1) shape functions
	static void compute_W0_grads
	(
		const position_type * corners, /**< [in] array of the global corner coordinates */
		MathVector<WDim> grad [numCorners] /**< [out] gradients of the Whitney-0 functions */
	);
	
/// gets the correspondence between the edge dof indices and the corners of the element
	static void get_edge_corners
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] the element */
		size_t edge_corner [numEdges] [2] /**< [out] edge dof -> corner of the element */
	);
	
public:
	
/// assembles the stiffness matrix of the rot-rot operator
	static void local_stiffness_and_mass
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		number S [maxNumEdges][maxNumEdges], /**< [out] local stiffness matrix */
		number M [maxNumEdges][maxNumEdges] /**< [out] local mass matrix */
	);
	
///	assembles the discrete weak div operator
	static void local_div_matrix
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		number B [][numEdges] /**< [out] local weak divergence operator matrix */
	);
	
///	computes the Nedelec shapes at a given point
	static void get_shapes
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		const MathVector<dim> local, /**< [in] local coordinates of the point where to compute */
		MathVector<WDim> shapes [] /**< [out] array for the shapes */
	);
	
///	computes of the values of the grid functions
	static void interpolate
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<dim> local [], /**< [in] local coordinates of the points where to compute */
		const size_t n_pnt, /**< [in] number of the points where to compute */
		MathVector<WDim> values [] /**< [out] where to store the computed n_pnt values */
	);
	
/// computes the curl of the grid functions (in 2d represented as \f$(z, 0)\f$ instead of \f$(0, 0, z)\f$)
	static void curl
	(
		const TDomain * domain, /**< [in] the domain */
		TElem * elem, /**< [in] element */
		const position_type * corners, /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		MathVector<WDim> & curl_vec /**< [out] where to store the computed curl */
	);
};

/// Specialization of NedelecT1_LDisc for triangles
/**
 * Triangle is considered as a simplex.
 * \see NedelecT1_LDisc and \see NedelecT1_LDisc_forSimplex
 */
template <typename TDomain>
class NedelecT1_LDisc<TDomain, Triangle> : public NedelecT1_LDisc_forSimplex<TDomain, Triangle>
{
};

/// Specialization of NedelecT1_LDisc for tetrahedra
/**
 * Tetrahedron is considered as a simplex.
 * \see NedelecT1_LDisc and \see NedelecT1_LDisc_forSimplex
 */
template <typename TDomain>
class NedelecT1_LDisc<TDomain, Tetrahedron> : public NedelecT1_LDisc_forSimplex<TDomain, Tetrahedron>
{
};

///@}

/// Interpolation of the Nedelec dofs and their curls
/**
 * This class implements the transformation of the Nedelec element into
 * the vector field.
 *
 * \tparam TDomain	type of the domain
 * \tparam refDim	dimensionality of the reference element
 * \tparam WDim		dimensionality of the domain (Do not specify it yourself! It is only for the specializations.)
 */
template <typename TDomain, int refDim, int WDim = TDomain::dim>
class NedelecInterpolation
{
public:

/// type of the geometric positions (WDim-vectors)
	typedef typename TDomain::position_type position_type;

/// max. number of the edges of the full-dimensional elements in the domain
	static const size_t maxNumEdges = (size_t) element_list_traits<typename domain_traits<WDim>::DimElemList>::maxEdges;

public:
/// computes the values at given points
	static void value
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<refDim> local [], /**< [in] local coordinates of the points where to compute */
		const size_t n_pnt, /**< [in] number of the points where to compute */
		MathVector<WDim> values [] /**< [out] where to store the computed n_pnt values */
	)
	{ // This is a generic version: It prints the error message. Cf. the specializations below.
		UG_THROW ("NedelecInterpolation::value: No implementation of the Nedelec elements for "
					"Reference Object " << elem->reference_object_id () <<
					" of reference dim. " << refDim << " in a " << WDim << "d domain.");
	}
	
/// computes curl of the function
	/**
	 * This function computes the value of the curl operator for the Nedelec
	 * representation. Curl is constant over elements. Note that the result
	 * is represented as a vector, not as a Whitney-2-form. In 2d, where the
	 * result has always the form (0, 0, z), it is represented as (z, 0).
	 */
	static void curl
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		MathVector<WDim> & curl_vec /**< [out] where to store the computed n_pnt values */
	)
	{ // This is a generic version: It prints the error message. Cf. the specializations below.
		UG_THROW ("NedelecInterpolation::curl: No implementation of the Nedelec elements for "
					"Reference Object " << elem->reference_object_id () <<
					" of reference dim. " << refDim << " in a " << WDim << "d domain.");
	}
	
///	computes flux of the curl through a given plane in an element
	/**
	 * This function computes the flux of the curl through a given plane inside
	 * of a given grid element. The plane is identified by a point on it and the
	 * normal. The flux is multiplied by the norm of the given normal vector (i.e.
	 * specify the unit normal to get the standard flux). The function returns 
	 * the area of the intersection if the element is intersected by the
	 * plane. Otherwise the function returns exactly 0.0.
	 */
	static number curl_flux
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<WDim> & normal, /**< [in] normal to the plane */
		const position_type pnt, /**< [in] point on the plane (identifying the plane) */
		number & flux /**< [out] the flux */
	)
	{ // This is a generic version: It prints the error message. Cf. the specializations below.
		UG_THROW ("NedelecInterpolation::curl_flux: No implementation of the Nedelec elements for "
					"Reference Object " << elem->reference_object_id () <<
					" of reference dim. " << refDim << " in a " << WDim << "d domain.");
		return false;
	}
};

/// A specialization of NedelecInterpolation for 2d, \see NedelecInterpolation
template <typename TDomain>
class NedelecInterpolation<TDomain, 2, 2>
{
public:

/// type of the geometric positions (WDim-vectors)
	typedef typename TDomain::position_type position_type;

/// max. number of the edges of the full-dimensional elements in the domain
	static const size_t maxNumEdges = (size_t) element_list_traits<typename domain_traits<2>::DimElemList>::maxEdges;

public:
/// computes the values at given points
	static void value
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<2> local [], /**< [in] local coordinates of the points where to compute */
		const size_t n_pnt, /**< [in] number of the points where to compute */
		MathVector<2> values [] /**< [out] where to store the computed n_pnt values */
	)
	{
		if (elem->reference_object_id () != ROID_TRIANGLE)
			UG_THROW ("No implementation of the Nedelec elements for "
						"Reference Object " << elem->reference_object_id () <<
						" of reference dim. 2 in a 2d domain. (This must be a triangle.)");
		NedelecT1_LDisc<TDomain, Triangle>::interpolate
			(domain, (Triangle *) elem, corners, dofs, local, n_pnt, values);
	}
	
/// computes curl of the function
	/**
	 * This function computes the value of the curl operator for the Nedelec
	 * representation. Curl is constant over elements. Note that the result
	 * is represented as a vector, not as a Whitney-2-form. In 2d, where the
	 * result has always the form (0, 0, z), it is represented as (z, 0).
	 */
	static void curl
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		MathVector<2> & curl_vec /**< [out] where to store the computed n_pnt values */
	)
	{
		if (elem->reference_object_id () != ROID_TRIANGLE)
			UG_THROW ("No implementation of the Nedelec elements for "
						"Reference Object " << elem->reference_object_id () <<
						" of reference dim. 2 in a 2d domain. (This must be a triangle.)");
		NedelecT1_LDisc<TDomain, Triangle>::curl
			(domain, (Triangle *) elem, corners, dofs, curl_vec);
	}
	
///	computes flux of the curl through a given plane in an element
	/**
	 * This function computes the flux of the curl through a given plane inside
	 * of a given grid element. The plane is identified by a point on it and the
	 * normal. The flux is multiplied by the norm of the given normal vector (i.e.
	 * specify the unit normal to get the standard flux). The function returns 
	 * the area of the intersection if the element is intersected by the
	 * plane. Otherwise the function returns exactly 0.0.
	 */
	static number curl_flux
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<2> & normal, /**< [in] normal to the plane */
		const position_type pnt, /**< [in] point on the plane (identifying the plane) */
		number & flux /**< [out] the flux */
	)
	{ // This is a generic version: It prints the error message. Cf. the specializations below.
		UG_THROW ("NedelecInterpolation::curl_flux: No implementation of the Nedelec elements for "
					"Reference Object " << elem->reference_object_id () <<
					" of reference dim. 3 in a 3d domain.");
		return false;
	}
};

/// A specialization of NedelecInterpolation for 3d, \see NedelecInterpolation
template <typename TDomain>
class NedelecInterpolation<TDomain, 3, 3>
{
public:

/// type of the geometric positions (WDim-vectors)
	typedef typename TDomain::position_type position_type;

/// max. number of the edges of the full-dimensional elements in the domain
	static const size_t maxNumEdges = (size_t) element_list_traits<typename domain_traits<3>::DimElemList>::maxEdges;

public:
/// computes the values at given points
	static void value
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<3> local [], /**< [in] local coordinates of the points where to compute */
		const size_t n_pnt, /**< [in] number of the points where to compute */
		MathVector<3> values [] /**< [out] where to store the computed n_pnt values */
	)
	{
		if (elem->reference_object_id () != ROID_TETRAHEDRON)
			UG_THROW ("No implementation of the Nedelec elements for "
						"Reference Object " << elem->reference_object_id () <<
						" of reference dim. 3 in a 3d domain. (This must be a tetrahedron.)");
		NedelecT1_LDisc<TDomain, Tetrahedron>::interpolate
			(domain, (Tetrahedron *) elem, corners, dofs, local, n_pnt, values);
	}
	
/// computes curl of the function
	/**
	 * This function computes the value of the curl operator for the Nedelec
	 * representation. Curl is constant over elements. Note that the result
	 * is represented as a vector, not as a Whitney-2-form.
	 */
	static void curl
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		MathVector<3> & curl_vec /**< [out] where to store the computed n_pnt values */
	)
	{
		if (elem->reference_object_id () != ROID_TETRAHEDRON)
			UG_THROW ("No implementation of the Nedelec elements for "
						"Reference Object " << elem->reference_object_id () <<
						" of reference dim. 3 in a 3d domain. (This must be a tetrahedron.)");
		NedelecT1_LDisc<TDomain, Tetrahedron>::curl
			(domain, (Tetrahedron *) elem, corners, dofs, curl_vec);
	}
	
///	computes flux of the curl through a given plane in an element
	/**
	 * This function computes the flux of the curl through a given plane inside
	 * of a given grid element. The plane is identified by a point on it and the
	 * normal. The flux is multiplied by the norm of the given normal vector (i.e.
	 * specify the unit normal to get the standard flux). The function returns 
	 * the area of the intersection if the element is intersected by the
	 * plane. Otherwise the function returns exactly 0.0.
	 */
	static number curl_flux
	(
		const TDomain * domain, /**< [in] the domain */
		GridObject * elem, /**< [in] element */
		const position_type corners [], /**< [in] array of the global corner coordinates */
		const number dofs [], /**< [in] arrays of values of the Nedelec degrees of freedom */
		const MathVector<3> & normal, /**< [in] normal to the plane */
		const position_type pnt, /**< [in] point on the plane (identifying the plane) */
		number & flux /**< [out] the flux */
	);
};

/**
 * Computes the arithm. average of given MathVectors.
 */
template <int dim>
inline void get_ave_vector
(
	const size_t n_vec, /**< [in] number of the vectors in the array */
	const MathVector<dim> * vec, /**< [in] array of the vectors */
	MathVector<dim> & ave /**< [out] the average */
)
{
	ave = 0;
	for (size_t i = 0; i < n_vec; i++) ave += vec [i];
	ave /= n_vec;
}

/** computes the "generalized vector product" of two vectors
 *
 * In 3d, this is a usual cross-product. In 2d, the first component
 * of the result is the determinant of the 2x2-matrix build of the operands,
 * and the second component is always 0.
 *
 * \tparam dim	dimensionality of the vectors
 */

template <size_t dim>
inline void GenVecCross
(
	MathVector<dim> & result,
	const MathVector<dim> & v_1, const MathVector<dim> & v_2
)
{
	/* Cf. the specializations below! */
	UG_THROW ("The generalized vector product is defined only in 2 and 3 dimensions");
};

template <>
inline void GenVecCross<2>
(
	MathVector<2> & result,
	const MathVector<2> & v_1, const MathVector<2> & v_2
)
{
	result[0] = v_1[0] * v_2[1] - v_1[1] * v_2[0];
	result[1] = 0;
};

template <>
inline void GenVecCross<3>
(
	MathVector<3> & result,
	const MathVector<3> & v_1, const MathVector<3> & v_2
)
{
	VecCross (result, v_1, v_2);
};

} // end namespace Electromagnetism
} // end namespace ug

#include "nedelec_local_ass_impl.h"

#endif /* __H__UG__PLUGINS__ELECTROMAGNETISM__ROT_ROT_ASS__ */

/* End of File */
