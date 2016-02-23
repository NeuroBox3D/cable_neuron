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
 * @param verbosity output verbosity: choose a value >1 for any output
 * @return whether domain is acyclic
 */
template <typename TDomain>
bool is_acyclic(SmartPtr<TDomain> dom, int verbosity);
template <typename TDomain>
bool is_acyclic(SmartPtr<TDomain> dom) {return is_acyclic(dom, 1);}


/**
 * @brief Checks the presynapse indices of the geometry.
 * The function checks whether the presynaptic indices are consecutive,
 * unique and start at 0. The return value is an error code specifying
 * either of the three:
 * 0: no error,
 * 2: at least one index was not unique,
 * 4: at least one index is not present (i.e., indices not consecutive).
 *
 * @param dom domain to be checked
 * @param verbosity output verbosity: choose a value >1 for any output
 * @return error code (0 for no errors)
 */
template <typename TDomain>
int check_presyn_indices(SmartPtr<TDomain> dom, int verbosity);
template <typename TDomain>
int check_presyn_indices(SmartPtr<TDomain> dom) {return check_presyn_indices(dom, 1);}


/**
 * @brief Performs several checks on the passed domain.
 * The function performs the is_acyclic() and check_presyn_indices() functions
 * and returns an error code containing the possible errors.
 *
 * 0: no error,
 * 1: domain contains a cycle,
 * 2: at least one index was not unique,
 * 4: at least one index is not present (i.e., indices not consecutive).
 *
 * @param dom domain to be checked
 * @param verbosity output verbosity: choose 0 for no output
 *                                           1 for short results
 *                                           >1 for more detailed output from each check
 * @return error code (0 for no errors)
 */
template <typename TDomain>
int check_domain(SmartPtr<TDomain> dom, int verbosity);
template <typename TDomain>
int check_domain(SmartPtr<TDomain> dom) {return check_domain(dom, 1);}


} // namespace cable_neruon
} // namespace ug



#endif // __UG__PLUGINS__CABLE_NEURON__UTIL__FUNCTIONS_H__
