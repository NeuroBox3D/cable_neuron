/*
 * synapse_distributor_impl.h
 *
 *  Created on: 08.04.2015
 *      Author: mbreit
 */

#include <common/error.h>

namespace ug {
namespace cable_neuron {

////////////////////////////////////////////////////////////////////////////////////////////
//	SynapseDistributor::has_active_synapses(MathVector<dim>& c, number t)
/**
 * Returns true/false whether the given position has an active synapse
 */
template <size_t dim>
bool SynapseDistributor::has_active_synapses(const MathVector<dim>& c, const number time, number& current) const
{
//	vector<SynapseInfo*> vSyn;
//	switch (dim)
//	{
//		case 1:
//			vSyn = find_synapses(c[0],0,0); break;
//		case 2:
//			vSyn = find_synapses(c[0],c[1],0); break;
//		case 3:
//			vSyn = find_synapses(c[0],c[1],c[2]); break;
//		default:
//			UG_THROW("Invalid dimension " << dim << "! Only dims 1, 2, 3 are allowed.");
//	}

///	TODO: Implement me!

	return false;
}


} // namespace cable_neuron
} // namespace ug
