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
