/*
 * synapse_iterator.h
 *
 *  Created on: Sep 22, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_SYNAPSE_ITERATOR_H_
#define SPLIT_SYNAPSE_HANDLER_SYNAPSE_ITERATOR_H_

#include <vector>
#include <cstddef>
#include "base_synapse.h"

#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "split_synapse_handler.h"

#include "alpha_post_synapse.h"
#include "alpha_pre_synapse.h"
#include "exp2_post_synapse.h"
#include "exp2_pre_synapse.h"

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

template <typename T, typename TDomain>
class SynapseIterator {
	typename std::vector<T*>::iterator m_intIt;
	typename std::vector<T*>::iterator m_end;
public:

	SynapseIterator(){}

	/*SynapseIterator(const std::vector<IBaseSynapse*>& v) {
		for(size_t i=0; i<v.size(); ++i){
			if(dynamic_cast<T*>( v[i])) {
				m_vSynapses.push_back((T*)v[i]);
			}
		}

		m_intIt = m_vSynapses.begin();
		m_end = m_vSynapses.end();
	}*/
	//void set_synapse_handler(SmartPtr<synapse_handler::SplitSynapseHandler<TDomain> > sh)
	virtual ~SynapseIterator(){}

	void set_synapse_handler(SmartPtr<synapse_handler::SplitSynapseHandler<TDomain> > sh) {
		std::vector<IBaseSynapse*> v = sh->get_synapses();
		std::vector<T*> vsyn;

		for(size_t i=0; i<v.size(); ++i){
			if(dynamic_cast<T*>( v[i])) {
				vsyn.push_back((T*)v[i]);
			}
		}

		m_intIt = vsyn.begin();
		m_end = vsyn.end();
	}

	T* operator*() {return *m_intIt;}
	T* operator++() {return *(++m_intIt);}
	bool operator==(const SynapseIterator<T, TDomain>& rhs) {return (rhs.m_intIt==m_intIt);}
	bool operator!=(const SynapseIterator<T, TDomain>& rhs) {return (rhs.m_intIt!=m_intIt);}
	bool is_referenceable() {return m_intIt!=m_end;}
};


// ////////////////////////////////////
//	explicit template instantiations //
// ////////////////////////////////////

#ifdef UG_DIM_1
	template class SynapseIterator<AlphaPreSynapse, Domain1d>;
	template class SynapseIterator<AlphaPostSynapse, Domain1d>;
	template class SynapseIterator<Exp2PreSynapse, Domain1d>;
	template class SynapseIterator<Exp2PostSynapse, Domain1d>;
#endif

#ifdef UG_DIM_2
	template class SynapseIterator<AlphaPreSynapse, Domain2d>;
	template class SynapseIterator<AlphaPostSynapse, Domain2d>;
	template class SynapseIterator<Exp2PreSynapse, Domain2d>;
	template class SynapseIterator<Exp2PostSynapse, Domain2d>;
#endif

#ifdef UG_DIM_3
	template class SynapseIterator<AlphaPreSynapse, Domain3d>;
	template class SynapseIterator<AlphaPostSynapse, Domain3d>;
	template class SynapseIterator<Exp2PreSynapse, Domain3d>;
	template class SynapseIterator<Exp2PostSynapse, Domain3d>;
#endif


} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_SYNAPSE_ITERATOR_H_ */
