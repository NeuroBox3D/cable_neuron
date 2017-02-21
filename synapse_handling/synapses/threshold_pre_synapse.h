/*
 * PreExp2Synapse.h
 *
 *  Created on: Mar 16, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__THRESHOLD_PRE_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__THRESHOLD_PRE_SYNAPSE_H

#include <string>           // std::string
#include "pre_synapse.h"    // IPreSynapse

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class ThresholdPreSynapse : public IPreSynapse
{
public:
    static const std::string name_string;

private:
	number m_onset;
	number m_duration;
	number m_threshold;

public:
	//ctor & dtor
	ThresholdPreSynapse(); // needed for template generator
	ThresholdPreSynapse(
			const number location,
			const number duration,
			const number threshold
			);

	ThresholdPreSynapse(
			const synapse_id id,
			const number location,
			const number duration,
			const number threshold
			);

	virtual ~ThresholdPreSynapse();

	// setter & getter
	void set_onset(number onset) {m_onset = onset;}
	void set_duration(number duration) {m_duration = duration;}
	void set_threshold(number val) {m_threshold = val;}

	number onset() const {return m_onset;}
	number duration() const {return m_duration;}
	number threshold() const {return m_threshold;}

	SynapseType type() const {return THRESHOLD_PRE_SYNAPSE;}
	const std::string& name() const {return name_string;}

	//interface methods
	void update(const number& t, const std::vector<number>& u);
	bool is_active(const number& t);

	//serialization interface methods
	void put_to(std::ostream& os) const;        //'put_to' == operator<<
	void get_from(std::istream& is);            //'get_from' == operator>>
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__THRESHOLD_PRE_SYNAPSE_H
