/*
 * functions.h
 *
 *  Created on: 17.02.2016
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__UTIL__FUNCTIONS_H__
#define __UG__PLUGINS__CABLE_NEURON__UTIL__FUNCTIONS_H__


#include "common/util/smart_pointer.h"


namespace ug {
namespace cable_neuron {


/**
 * @brief Checks whether the given domain contains a cycle.
 * If it does, the coordinates of a vertex contained in the cycle are UG_LOGged
 * to help identify the cycle.
 * For this function to work, the domain should ideally not be distributed.
 * However, the function will still yield correct results if the domain is distributed,
 * as long as no neuron is cut by the distribution. (This condition is not checked!)
 *
 * @param dom domain to be checked
 * @return whether domain is acyclic
 */
template <typename TDomain>
bool is_acyclic(SmartPtr<TDomain> dom);



} // namespace cable_neruon
} // namespace ug



#endif // __UG__PLUGINS__CABLE_NEURON__UTIL__FUNCTIONS_H__
