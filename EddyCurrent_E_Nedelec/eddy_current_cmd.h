/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
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
 * Misc. functions and commands for the eddy current simulations.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_CMD__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_CMD__

#include <vector>

namespace ug{
namespace Electromagnetism{

/// Computes the power of the electromagnetic field (up to the contribution of the boundary)
template <typename TGridFunc>
void calc_power
(
	TGridFunc * pJGGF, ///< grid function of the Nedelec-DoFs of the generator current \f$\mathbf{J}_G\f$
	size_t JG_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{J}_G\f$
	SubsetGroup & JG_ssg, ///< (full-dim.) subsets where \f$\mathbf{J}_G\f$ is defined (non-zero and in the kernel)
	TGridFunc * pEGF, ///< grid function of the Nedelec-DoFs of the electric field \f$\mathbf{E}\f$
	size_t E_fct[], ///< indices of the Re and Im parts in the grid function for \f$\mathbf{E}\f$
	number pow[] ///< to add the integral
);

/// Prints the (complex-valued) power of the electromagnetic field
template <typename TGridFunc>
void CalcPower
(
	SmartPtr<TGridFunc> spJGGF, ///< [in] grid function with the generator current
    const char* JG_cmps, ///< [in] names of the components of the grid function (for Re and Im)
	const char* JG_ss, ///< (full-dim.) subsets where \f$\mathbf{J}_G\f$ is defined (non-zero and in the kernel)
	SmartPtr<TGridFunc> spEGF, ///< [in] grid function with the electric field
    const char* E_cmps ///< [in] names of the components of the grid function (for Re and Im)
);

/// Computation of the magnetic flux through windings of a coil
template <typename TGridFunc>
number calc_magnetic_flux
(
	const TGridFunc * gfE, ///< [in] grid function with the electric field
	const size_t fct [2], ///< [in] function components (Re and Im)
	const SubsetGroup & ss_grp, ///< [in] subsets of the kerner of the coil
	const MathVector<TGridFunc::dim> & Normal, ///< [in] normal to the planes of the windings
	const typename TGridFunc::domain_type::position_type & base_pnt, ///< [in] point on the plane of the 1st winding
	const size_t n_pnt, ///< [in] number of the windings
	const typename TGridFunc::domain_type::position_type & d_pnt, ///< [in] thickness of the windings
	number flux [2] ///< [out] the computed flux (Re and Im)
);

/// Prints of the magnetic flux through windings of a coil
template <typename TGridFunc>
void CalcMagneticFlux
(
	SmartPtr<TGridFunc> spGF, ///< [in] grid function with the electric field
    const char* cmps, ///< [in] names of the components of the grid function (for Re and Im)
	const char* subsets, ///< [in] subsets where to compute the flux
	const std::vector<number>& Normal, ///< [in] normal to the planes of the windings
	const std::vector<number>& base_pnt, ///< [in] point on the plane of the 1st winding
	const size_t n_pnt, ///< [in] number of the windings
	const std::vector<number>& d_pnt ///< [in] thickness of the windings
);

} // end namespace Electromagnetism
} // end namespace ug

#include "eddy_current_cmd_impl.h"

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_CMD__

/* End of File */
