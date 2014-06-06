/*
 * Misc. functions and commands for the eddy current simulations.
 *
 * Created on: 04.06.2014
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_AUX_CMD__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_AUX_CMD__

namespace ug{
namespace Electromagnetism{

/**
 * Computation of the flux of a vector field (given by the Whitney-1 elements)
 * through a low-dimensional subset (surface). The function returns the value
 * of the flux.
 */
template <typename TGridFunc>
number ComputeFlux
(
	TGridFunc * pGF, ///< grid function of the Nedelec-DoFs of the vector field
	size_t fct, ///< index of the function in the grid function
	SubsetGroup & volSSG, ///< full-dim. subsets (adjacent to the surface) to indicate the negative direction
	SubsetGroup & faceSSG ///< the surface (the low-dim. subsets)
);

/**
 * Computation of the flux of a vector field (given by the Whitney-1 elements)
 * through a low-dimensional subset (surface). The function prints the value
 * of the flux in the shell.
 */
template <typename TGridFunc>
void ComputeFlux
(
	SmartPtr<TGridFunc> spGF, ///< grid function of the Nedelec-DoFs of the vector field
	const char * fct_name, ///< name of the function in the grid function
	const char * vol_subsets, ///< full-dim. subsets (adjacent to the surface) to indicate the negative direction
	const char * face_subsets ///< the surface (the low-dim. subsets)
);

} // end namespace Electromagnetism
} // end namespace ug

#include "nedelec_aux_cmd_impl.h"

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__NEDELEC_AUX_CMD__

/* End of File */
