/*
 * Misc. functions and commands for the eddy current simulations.
 *
 * Created on: 26.05.2014
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_CMD__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_CMD__

#include <vector>

namespace ug{
namespace Electromagnetism{

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
