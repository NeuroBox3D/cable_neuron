/*
 * IPostSynapse.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lreinhardt
 */

#ifndef SPLIT_SYNAPSE_HANDLER_POST_SYNAPSE_H_
#define SPLIT_SYNAPSE_HANDLER_POST_SYNAPSE_H_

#include <common/types.h> 											//number
#include "../synapse_handler/function/types.h" 						//SynapseType
#include <string>													//std::string

#include "base_synapse.h"											//IBaseSynapse
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 	//VectorProxyBase

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

class IPostSynapse : public IBaseSynapse
{
public:
	//ctor & dtor
	IPostSynapse(); 							//needed for SynapseDealer
	IPostSynapse(
			const SYNAPSE_ID id,
			const number location
	);

	virtual ~IPostSynapse();

	//setter & getter

	virtual bool is_presynapse() const {return false;}
	virtual bool is_postsynapse() const {return true;}

	virtual number current(const number& t, const number& vm) = 0;

	//from serialization interface IBaseSynapse
	virtual void put_to(std::ostream& os) const = 0;			//'put_to' == operator<<
	virtual void get_from(std::istream& is) = 0;				//'get_from' == operator>>

	virtual void activate(number time) = 0;
	virtual void deactivate(number time) = 0;
	virtual bool is_active(number time) = 0;
};

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */

#endif /* SPLIT_SYNAPSE_HANDLER_POST_SYNAPSE_H_ */
