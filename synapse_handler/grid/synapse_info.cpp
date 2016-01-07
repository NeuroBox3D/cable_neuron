/*!
 * \file synapse_info.cpp
 * \brief synapse info struct
 *
 *  Created on: May 11, 2015
 *      Author: stephan
 */

#include "synapse_info.h"

namespace ug {
namespace synapse_handler {


void SynapseInfo::print_to_log() const
{
	switch (this->m_type)
	{
		case EMPTY_SYNAPSE:
			UG_LOG("type: EMPTY_SYNAPSE" << std::endl);
			return;

		case ALPHA_SYNAPSE:
			typedef synapse_traits<AlphaSynapse> STA;
			UG_LOG("type        : " << STA::name() << std::endl);
			UG_LOG("local coords: " << STA::loc_coord(*this) << std::endl);
			UG_LOG("onset       : " << STA::onset(*this) << std::endl);
			UG_LOG("max conduct.: " << STA::g_max(*this) << std::endl);
			UG_LOG("revers. pot.: " << STA::v_rev(*this) << std::endl);
			UG_LOG("tau         : " << STA::tau(*this) << std::endl);
			UG_LOG("frequency   : " << STA::freq(*this) << std::endl);
			UG_LOG("#spikes     : " << STA::nSpikes(*this) << std::endl);
			return;

		case EXP2_SYNAPSE:
			typedef synapse_traits<Exp2Syn> STB;
			UG_LOG("type        : " << STB::name() << std::endl);
			UG_LOG("local coords: " << STB::loc_coord(*this) << std::endl);
			UG_LOG("onset       : " << STB::onset(*this) << std::endl);
			UG_LOG("max conduct.: " << STB::g_max(*this) << std::endl);
			UG_LOG("revers. pot.: " << STB::v_rev(*this) << std::endl);
			UG_LOG("tau1        : " << STB::tau1(*this) << std::endl);
			UG_LOG("tau2        : " << STB::tau2(*this) << std::endl);
			UG_LOG("delay       : " << STB::delay(*this) << std::endl);
			UG_LOG("threshold   : " << STB::threshold(*this) << std::endl);
			UG_LOG("presyn. ind.: " << STB::presyn_ind(*this) << std::endl);
			UG_LOG("activated?  : " << (STB::activated(*this) ? "true" : "false") << std::endl);
			UG_LOG("duration    : " << STB::activity_time(*this) << std::endl);
			return;

		case CUSTOM_SYNAPSE:
			UG_LOG("type: CUSTOM_SYNAPSE" << std::endl);
			return;

		default:
			UG_THROW("Unknown synapse type.");
	}
}

std::string SynapseInfo::format() const {
	std::ostringstream os;
	switch (this->m_type)
	{
		case EMPTY_SYNAPSE:
			os << "type: EMPTY_SYNAPSE" << std::endl;
			break;

		case ALPHA_SYNAPSE:
			typedef synapse_traits<AlphaSynapse> STA;
			os << "type        : " << STA::name() << std::endl;
			os << "local coords: " << STA::loc_coord(*this) << std::endl;
			os << "onset       : " << STA::onset(*this) << std::endl;
			os << "max conduct.: " << STA::g_max(*this) << std::endl;
			os << "revers. pot.: " << STA::v_rev(*this) << std::endl;
			os << "tau         : " << STA::tau(*this) << std::endl;
			os << "frequency   : " << STA::freq(*this) << std::endl;
			os << "#spikes     : " << STA::nSpikes(*this) << std::endl;
			break;

		case EXP2_SYNAPSE:
			typedef synapse_traits<Exp2Syn> STB;
			os << "type        : " << STB::name() << std::endl;
			os << "local coords: " << STB::loc_coord(*this) << std::endl;
			os << "onset       : " << STB::onset(*this) << std::endl;
			os << "max conduct.: " << STB::g_max(*this) << std::endl;
			os << "revers. pot.: " << STB::v_rev(*this) << std::endl;
			os << "tau1        : " << STB::tau1(*this) << std::endl;
			os << "tau2        : " << STB::tau2(*this) << std::endl;
			os << "delay       : " << STB::delay(*this) << std::endl;
			os << "threshold   : " << STB::threshold(*this) << std::endl;
			os << "presyn. ind.: " << STB::presyn_ind(*this) << std::endl;
			os << "activated?  : " << (STB::activated(*this) ? "true" : "false") << std::endl;
			os << "duration    : " << STB::activity_time(*this) << std::endl;
			break;

		case CUSTOM_SYNAPSE:
			os << "type: CUSTOM_SYNAPSE" << std::endl;
			break;

		default:
			os << "Unknown synapse type." << std::endl;
			break;

	}
	return os.str();

}

std::ostream& operator<<(std::ostream& os, const std::vector<SynapseInfo>& synapses) {
	std::vector<SynapseInfo>::const_iterator it = synapses.begin();
	for (; it != synapses.end(); ++it)
		os << it->format();
	return os;
}


// DO NOT CHANGE! Needed for serialization! //
std::ostream& operator<<(std::ostream &os, const SynapseInfo& synapse)
{
	using std::ostringstream;
	ostringstream strs;
	strs << synapse.m_locCoords << " ";
	strs << static_cast<int>(synapse.m_type) << " ";
	strs << synapse.m_onset << " ";
	strs << synapse.m_tau << " ";
	strs << synapse.m_gMax << " ";
	strs << synapse.m_vRev << " ";
	strs << synapse.m_param1 << " ";
	strs << synapse.m_param2 << " ";
	strs << synapse.m_param3 << " ";
	strs << synapse.m_param4;
	os << strs.str();
	return os;
}


// DO NOT CHANGE! Needed for serialization! //
std::istream& operator>>(std::istream& in, SynapseInfo& synapse)
{
	std::string temp;
	using boost::lexical_cast;
	// loc coordinates
	in >> temp;
	number locCoords = lexical_cast<number>(temp);
	temp.clear();

	// synapse type
	in >> temp;
	// cast to int first, since a direct lexical_cast to char would yield
	// '1' (aka 49) from "1" instead of 1 etc.
	unsigned char type = static_cast<unsigned char>(lexical_cast<unsigned int>(temp));
	temp.clear();

	// onset
	in >> temp;
	number onset = lexical_cast<number>(temp);
	temp.clear();

	// tau
	in >> temp;
	number tau = lexical_cast<number>(temp);
	temp.clear();

	// gmax
	in >> temp;
	number gmax = lexical_cast<number>(temp);
	temp.clear();

	// vrev
	in >> temp;
	number vrev = lexical_cast<number>(temp);
	temp.clear();

	// param 1
	in >> temp;
	number param1 = lexical_cast<number>(temp);
	temp.clear();

	// param 2
	in >> temp;
	number param2 = lexical_cast<number>(temp);
	temp.clear();

	// param 3
	in >> temp;
	number param3 = lexical_cast<number>(temp);
	temp.clear();

	// param 4
	in >> temp;
	unsigned int param4 = lexical_cast<unsigned int>(temp);
	temp.clear();

	// populate synapse
	synapse = SynapseInfo(locCoords, type, onset, tau, gmax, vrev, param1, param2, param3, param4);
	return in;
}


} // namespace sh
} // namespace ug

