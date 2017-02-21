/*
 * AlphaPreSynapse.h
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__ONSET_PRE_SYNAPSE_H
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__ONSET_PRE_SYNAPSE_H

#include <string>           // std::string
#include "pre_synapse.h"    // IPreSynapse


namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class OnsetPreSynapse : public IPreSynapse
{
public:
    static const std::string name_string;

private:
	number m_onset;
	number m_duration;

public:
	//ctor & dtor
	OnsetPreSynapse();					//needed for template generator
	OnsetPreSynapse(
			const number location,
			const number onset,
			const number duration);

	OnsetPreSynapse(
			const synapse_id id,
			const number location,
			const number onset,
			const number duration);

	virtual ~OnsetPreSynapse();

	//setter & getter
	SynapseType type() const {return ONSET_PRE_SYNAPSE;}
	const std::string&  name() const {return name_string;}

	void set_onset(number onset) {m_onset = onset;}
	void set_duration(number dur) {m_duration = dur;}

	number onset() const {return m_onset;}
	number duration() const {return m_duration;}

	//Interface dummy's
	void update(const number& t, const std::vector<number>& u);
	bool is_active(const number& t);

	//serialization from IBaseSynapse interface
	void put_to(std::ostream& os) const;			//'put_to' == operator<<
	void get_from(std::istream& is);				//'get_from' == operator>>
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__ONSET_PRE_SYNAPSE_H
