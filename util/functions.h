/*
 * functions.h
 *
 *  Created on: 17.02.2016
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__UTIL__FUNCTIONS_H__
#define __UG__PLUGINS__CABLE_NEURON__UTIL__FUNCTIONS_H__


#include "common/util/smart_pointer.h"
#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_multi_grid.h"


namespace ug {
namespace cable_neuron {


///@{
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
///@}


///@{
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
///@}


///@{
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
///@}


// This is a debugging tool for a very specific problem.
// Do not compile but for the specific debuging purpose.
#if 0
template <typename TDomain>
void test_vertices(SmartPtr<TDomain> dom);
#endif

/**
 * @brief Attach unique IDs to neurons in a network.
 * In a depth-first forest traversal, all vertices of a network are assigned
 * a unique neuron identifier corresponding to the neuron they belong to.
 * Numbering starts at 0 and is without ID gaps.
 *
 * @param g   grid to perform identification on
 */
void neuron_identification(Grid& g);


/**
 * @brief Save one neuron from a network to SWC file.
 * The function will assign each neuron in the network a unique identifier
 * using the neuron_identification() function. It then exports the neuron
 * with the specified identifier to the specified file name.
 * Optionally, a scaling factor can be provided that will be applied to
 * the coordinates and radius prior to exporting.
 *
 * @param ugxFileName   network source file (.ugx)
 * @param neuronIndex   neuron ID to be exported
 * @param swcFileName   output file name
 * @param scale         scaling factor (optional)
 */
void save_neuron_to_swc
(
    std::string ugxFileName,
    size_t neuronIndex,
    std::string swcFileName,
    number scale = 1.0
);



/// length of neurite subset
number subset_length(int si, ConstSmartPtr<MGSubsetHandler> sh);

/// length of neurite subset
number subset_length(const char* subset, ConstSmartPtr<MGSubsetHandler> sh);


} // namespace cable_neruon
} // namespace ug



#endif // __UG__PLUGINS__CABLE_NEURON__UTIL__FUNCTIONS_H__
