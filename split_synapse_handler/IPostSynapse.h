/*
 * IPostSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_IPOSTSYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_IPOSTSYNAPSE_H_

#include <common/types.h> 											//number
#include "../synapse_handler/function/types.h" 						//SynapseType
#include <string>													//std::string
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 	//VectorProxyBase
#include "IBaseSynapse.h"											//IBaseSynapse

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPostSynapse : public IBaseSynapse
{
private:
	unsigned long long m_id;					//own postsynapse id / alternative: shared id if there is a 1:1 relation?
	unsigned long long m_presynapse_id;			//presynapse id (if needed/existing)
	number m_location;							//location on edge

public:
	//ctor & dtor
	IPostSynapse() {} 							//needed for SynapseDealer
	IPostSynapse(
			const unsigned long long id,
			const unsigned long long presynapseid,
			const number& location);

	virtual ~IPostSynapse();

	//setter & getter
	void set_presynapse_id(const unsigned long long id) {m_presynapse_id = id;}
	void set_id(const unsigned long long id) {m_id = id;}
	void set_location(const number& loc) {m_location = loc;}

	unsigned long long presynapse_id() const {return m_presynapse_id;}
	unsigned long long id() const {return m_id;}
	number location() const {return m_location;}

	//virtual interface methods
	virtual SynapseType type() const = 0;
	virtual std::string name() const = 0;

	virtual number current(const number& t) = 0;

	//from serialization interface IBaseSynapse
	virtual void put_to(std::ostream& os) const = 0;			//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;				//'get_from' == operator>>
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_IPOSTSYNAPSE_H_ */
