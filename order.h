/*
 * order.h
 *
 *  Created on: 12.08.2015
 *      Author: mbreit
 */

#ifndef PLUGINS_EXPERIMENTAL_HH_KABELNEW_ORDER_H_
#define PLUGINS_EXPERIMENTAL_HH_KABELNEW_ORDER_H_

#include <vector>

#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/dof_manager/dof_distribution.h"


namespace ug {
namespace cable {


void compute_cuthillmckee_order
(
	std::vector<size_t>& vNewIndex,
	std::vector<std::vector<size_t> >& vvConnection
);

void order_cuthillmckee(DoFDistribution& dofDistr);

template <typename TDomain>
void order_cuthillmckee(ApproximationSpace<TDomain>& approxSpace);


} // namespace cable
} // namespace ug


#endif // PLUGINS_EXPERIMENTAL_HH_KABELNEW_ORDER_H_
