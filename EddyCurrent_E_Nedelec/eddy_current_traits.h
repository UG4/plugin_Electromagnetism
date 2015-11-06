/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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
 * Auxiliary definitions for the eddy currend model.
 */
#ifndef __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_TRAITS__
#define __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_TRAITS__

namespace ug{
namespace Electromagnetism{

/// Auxiliary class defining some important constants
class EddyCurrentTraits
{
public:
	
/// index of the real part in the grid functions
	static const size_t _Re_ = 0;
/// index of the imaginary part in the grid functions
	static const size_t _Im_ = 1;
};
	
} // end namespace Electromagnetism
} // end namespace ug

#endif // __H__UG__PLUGINS__ELECTROMAGNETISM__EDDY_CURRENT_TRAITS__

/* End of File */
