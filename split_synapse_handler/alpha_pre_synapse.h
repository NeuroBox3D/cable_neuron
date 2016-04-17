/*
 * PreAlphaSynapse.h
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#ifndef PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_
#define PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_

#include "../synapse_handler/function/types.h" 						//SynapseType
#include <string>													//std::string
#include <common/types.h> 											//number
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 	//VectorProxyBase
#include "pre_synapse.h"					//IPreSynapse
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "synapse_container.h"


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class AlphaPreSynapse : public IPreSynapse
{
private:
	number m_onset;

public:
	//ctor & dtor
	AlphaPreSynapse();					//needed for template generator
	AlphaPreSynapse(
			const number& location,
			const number& onset);

	AlphaPreSynapse(
			const unsigned long long id,
			const unsigned long long postsynapse_id,
			const number& location,
			const number& onset);

	virtual ~AlphaPreSynapse();

	//setter & getter
	SynapseType type() const {return ALPHA_PRE_SYNAPSE;}
	std::string name() const {return "ALPHA_PRE_SYNAPSE";}

	//pre synapses are true
	bool split_type() const {return true;}

	void set_onset(const number& onset) {m_onset = onset;}
	number onset() const {return m_onset;}

	//Interface dummy's
	void update(const number& t, VectorProxyBase* up=NULL);
	bool is_active(const number& t, VectorProxyBase* up=NULL);

	//serialization from IBaseSynapse interface
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>
};

/**
 * class for parametrization of a set of AlphaPreSynapses, for usage in lua-scripts
 */
class AlphaPreSynapses : ISynapseContainer
{
private:
	number m_mean_onset;
	number m_dev_onset;

	std::vector<IBaseSynapse*> m_vSynapses;

	bool m_parametrizised;

public:
	AlphaPreSynapses(const size_t num_synapses)
:m_mean_onset(0),
 m_dev_onset(0),
 m_vSynapses(std::vector<IBaseSynapse*>(num_synapses)),
 m_parametrizised(true)
{}

	AlphaPreSynapses(
			const size_t num_synapses,
			number mean_onset,
			number dev_onset)
:m_mean_onset(mean_onset),
 m_dev_onset(dev_onset),
 m_vSynapses(std::vector<IBaseSynapse*>(num_synapses)),
 m_parametrizised(true)
{}

	void set_mean_onset(const number val) {m_mean_onset=val;}
	void set_dev_onset(const number val) {m_dev_onset=val;}
	size_t size() {return m_vSynapses.size();}

	virtual ~AlphaPreSynapses() {}
	/**
	 * Call this in your lua script to assemble the synapses vector
	 */
	std::vector<IBaseSynapse*> get_synapses() {

		//return m_vSynapses;

		if(!m_parametrizised)
			UG_THROW("No parameters for AlphaPreSynapses");

		boost::mt19937 onset_rng;
		boost::normal_distribution<number> onset_dist(m_mean_onset, m_dev_onset);
		boost::variate_generator<boost::mt19937, boost::normal_distribution<number> > onset_var(onset_rng, onset_dist);

		for(size_t i = 0; i < m_vSynapses.size(); ++i) {
			IBaseSynapse *s = new AlphaPreSynapse(0.0, onset_var());
			m_vSynapses[i] = s;
		}

		return m_vSynapses;
	}
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* PLUGINS_CABLE_NEURON_SPLITTEDSYNAPSE_HANDLER_PREALPHASYNAPSE_H_ */
