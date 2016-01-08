/*
 * order.h
 *
 *  Created on: 12.08.2015
 *      Author: mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__UTIL__ORDER_H__
#define __UG__PLUGINS__CABLE_NEURON__UTIL__ORDER_H__

#include <vector>

#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/dof_manager/dof_distribution.h"


namespace ug {
namespace cable_neuron {


void compute_cuthillmckee_order
(
	std::vector<size_t>& vNewIndex,
	std::vector<std::vector<size_t> >& vvConnection
);

void order_cuthillmckee(DoFDistribution& dofDistr);

template <typename TDomain>
void order_cuthillmckee(ApproximationSpace<TDomain>& approxSpace);


} // namespace cable_neuron
} // namespace ug


#endif // __UG__PLUGINS__CABLE_NEURON__UTIL__ORDER_H__
