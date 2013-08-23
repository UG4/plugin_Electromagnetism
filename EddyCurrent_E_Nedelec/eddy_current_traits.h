/*
 * Auxiliary definitions for the eddy currend model.
 * Created on: 15.07.2013
 * Author: D. Logashenko
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
