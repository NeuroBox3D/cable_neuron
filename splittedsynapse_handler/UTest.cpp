/*
 * UTest.cpp
 *
 *  Created on: Mar 1, 2016
 *      Author: lreinhardt
 */

#include "UTest.h"
#include "IPreSynapse.h"
#include "PreAlphaSynapse.h"
#include <iostream>

using namespace std;

namespace ug {
namespace cable_neuron {
namespace synapse_handler {

UTest::UTest() {
	IPreSynapse* s1 = new PreAlphaSynapse(0.1);
	cout<<s1->type()<<endl;
	cout<<s1->name()<<endl;
	delete s1;

}

UTest::~UTest() {
	// TODO Auto-generated destructor stub
}

} /* namespace synapse_handler */
} /* namespace cable_neuron */
} /* namespace ug */
