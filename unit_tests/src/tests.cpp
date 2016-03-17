/*!
 * \file plugins/experimental/synapse_handler/unit_tests/src/tests.cpp
 * \brief boost unit test for the synapse_handler plugin
 *
 * \author Stephan Grein
 */

#define BOOST_TEST_MODULE __CPP__UNIT_TESTS__UG__SYNAPSE_HANDLER__

// boost includes
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

// ug includes
#include "../../synapse_handler/synapse_handler.h"
#include "../../synapse_handler/synapse.h"
#include "../../splittedsynapse_handler/IPreSynapse.h"
#include "../../splittedsynapse_handler/PreAlphaSynapse.h"
#include "../../splittedsynapse_handler/PreExp2Synapse.h"
#include "../../splittedsynapse_handler/IPostSynapse.h"
#include "../../splittedsynapse_handler/PostAlphaSynapse.h"
#include "../../splittedsynapse_handler/PostExp2Synapse.h"
#include "../../synapse_handler/function/types.h" 				//SynapseType


// using directives
using namespace boost::unit_test;
using namespace ug::cable_neuron;
using namespace synapse_handler;
using namespace std;

/// begin of test suite presynapses
BOOST_AUTO_TEST_SUITE(PRESYNAPSES);

BOOST_AUTO_TEST_CASE(ALPHASYNAPSE) {
	number onset = 1e-6;
	number location = 0.3;
	number t1 = 1e-7;
	number t2 = 3e-6;

	IPreSynapse *s1 = new PreAlphaSynapse(0.0, 0.0);
	IPreSynapse *s2 = new PreAlphaSynapse(0.0, 0.0);
	IPreSynapse *s3 = new PreAlphaSynapse(3, 0.0, 0.0);


	BOOST_REQUIRE_MESSAGE(s1->name() == "PRE_ALPHA_SYNAPSE","s1 is a PRE_ALPHA_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s1->type() == PRE_ALPHA_SYNAPSE,"s1's type is PRE_ALPHA_SYNAPSE");

	s1->set_id(0);
	s2->set_id(1);
	BOOST_REQUIRE_MESSAGE(s1->id() == 0,"s1's id is 0");
	BOOST_REQUIRE_MESSAGE(s2->id() == 1,"s2's id is 1");
	BOOST_REQUIRE_MESSAGE(s3->id() == 3,"s3's id is 3");

	s1->set_location(location);
	BOOST_REQUIRE_MESSAGE( fabs(s1->location() - location) < 1e-16,"location setter and getter test"); //"perils of floating point comparisons"

	static_cast<PreAlphaSynapse*>(s1)->set_onset(onset);
	BOOST_REQUIRE_MESSAGE( fabs(static_cast<PreAlphaSynapse*>(s1)->onset() -onset) < 1e-16 ,"onset setter and getter test");

	BOOST_REQUIRE_MESSAGE(s1->active(t1) == false ,"synapse should not be active");
	BOOST_REQUIRE_MESSAGE(s1->active(t2) == true ,"synapse should be active");

	//todo:
	//test update() somehow


	delete s1;
	delete s2;
	delete s3;

	BOOST_MESSAGE("PreAlphaSynapse working.");
	cout<<"PreAlphaSynapse working."<<endl;
}

BOOST_AUTO_TEST_CASE(EXP2SYNAPSE) {
	number location = 0.3;
	number onset = 1e-6;
	number t1 = 1e-7;
	number t2 = 3e-6;

	IPreSynapse *s1 = new PreExp2Synapse(0,0);
	IPreSynapse *s2 = new PreExp2Synapse(1,0,0);
	IPreSynapse *s3 = new PreExp2Synapse(0,0);

	s3->set_id(2);
	BOOST_REQUIRE_MESSAGE(s1->id() == 0,"id check");
	BOOST_REQUIRE_MESSAGE(s2->id() == 1,"id check");
	BOOST_REQUIRE_MESSAGE(s3->id() == 2,"id check");

	s1->set_location(location);
	BOOST_REQUIRE_MESSAGE( fabs(s1->location() - location) < 1e-16,"location setter and getter test"); //"perils of floating point comparisons"

	BOOST_REQUIRE_MESSAGE(s1->name() == "PRE_EXP2_SYNAPSE","name check");
	BOOST_REQUIRE_MESSAGE(s2->type() == PRE_EXP2_SYNAPSE,"type check");

	static_cast<PreExp2Synapse*>(s1)->set_onset(onset);
	BOOST_REQUIRE_MESSAGE( fabs(static_cast<PreExp2Synapse*>(s1)->onset() -onset) < 1e-16 ,"onset setter and getter test");

	BOOST_REQUIRE_MESSAGE(s1->active(t1) == false ,"synapse should not be active");
	BOOST_REQUIRE_MESSAGE(s1->active(t2) == true ,"synapse should be active");

	//todo:
	//test update() somehow

	delete s1;
	delete s2;
	delete s3;

	BOOST_MESSAGE("PreExp2Synapse working.");
	cout<<"PreExp2Synapse working."<<endl;
}

/// end of test suite presynapses
BOOST_AUTO_TEST_SUITE_END();




/// begin of test suite postsynapses
BOOST_AUTO_TEST_SUITE(POSTSYNAPSES);

BOOST_AUTO_TEST_CASE(ALPHASYNAPSE) {
	number gmax = 1;
	number onset = 2;
	number tau = 3;
	number vm = 4;
	number e = 5;
	number location = 6;

	unsigned long id1 = 0;
	unsigned long id2 = 1;

	IPostSynapse *s1 = new PostAlphaSynapse(0,0,0,0,0,0,0);
	IPostSynapse *s2 = new PostAlphaSynapse(0,0,0,0,0,0,0);

	BOOST_REQUIRE_MESSAGE(s1->name() == "POST_ALPHA_SYNAPSE","s1 is a POST_ALPHA_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s1->type() == POST_ALPHA_SYNAPSE,"s1's type is POST_ALPHA_SYNAPSE");

	s1->set_presynapse_id(id1);
	s2->set_presynapse_id(id2);
	s1->set_location(location);
	static_cast<PostAlphaSynapse*>(s1)->set_gMax(gmax);
	static_cast<PostAlphaSynapse*>(s1)->set_onset(onset);
	static_cast<PostAlphaSynapse*>(s1)->set_tau(tau);
	static_cast<PostAlphaSynapse*>(s1)->set_vm(vm);
	static_cast<PostAlphaSynapse*>(s1)->set_e(e);

	BOOST_REQUIRE_MESSAGE(s1->presynapse_id() == id1, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s2->presynapse_id() == id2, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(fabs(s1->location() - location) < 1e-16, "location setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostAlphaSynapse*>(s1)->gMax() - gmax) < 1e-16 ,"gmax setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostAlphaSynapse*>(s1)->onset() - onset) < 1e-16 ,"onset setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostAlphaSynapse*>(s1)->tau() - tau) < 1e-16 ,"tau setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostAlphaSynapse*>(s1)->vm() - vm) < 1e-16 ,"vm etter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostAlphaSynapse*>(s1)->e() - e) < 1e-16 ,"e setter and getter test");

	//todo:
	//test current?

	delete s1;
	delete s2;

	BOOST_MESSAGE("PostAlphaSynapse working.");
	cout<<"PostAlphaSynapse working."<<endl;
}

BOOST_AUTO_TEST_CASE(EXP2SYNAPSE) {
	number tau1 = 2;
	number tau2 = 3;
	number w = 1;
	number vm = 4;
	number e = 5;
	number location = 6;

	unsigned long id1 = 0;
	unsigned long id2 = 1;
	unsigned long id3 = 2;

	IPostSynapse *s1 = new PostExp2Synapse(0,0,0,0,0,0);
	IPostSynapse *s2 = new PostExp2Synapse(id2,0,0,0,0,0,0);
	IPostSynapse *s3 = new PostExp2Synapse(0,0,0,0,0,0);

	s3->set_presynapse_id(id3);
	BOOST_REQUIRE_MESSAGE(s1->presynapse_id() == id1, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s2->presynapse_id() == id2, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s3->presynapse_id() == id3, "presynapse id getter and setter test");

	BOOST_REQUIRE_MESSAGE(s1->name() == "POST_EXP2_SYNAPSE","s1 is a POST_ALPHA_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s1->type() == POST_EXP2_SYNAPSE,"s1's type is POST_ALPHA_SYNAPSE");

	static_cast<PostExp2Synapse*>(s1)->set_tau1(tau1);
	static_cast<PostExp2Synapse*>(s1)->set_w(w);
	static_cast<PostExp2Synapse*>(s1)->set_tau2(tau2);
	static_cast<PostExp2Synapse*>(s1)->set_vm(vm);
	static_cast<PostExp2Synapse*>(s1)->set_e(e);
	s1->set_location(location);
	BOOST_REQUIRE_MESSAGE(fabs(s1->location() - location) < 1e-16, "location setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostExp2Synapse*>(s1)->tau1() - tau1) < 1e-16 ,"gmax setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostExp2Synapse*>(s1)->tau2() - tau2) < 1e-16 ,"onset setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostExp2Synapse*>(s1)->w() - w) < 1e-16 ,"tau setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostExp2Synapse*>(s1)->vm() - vm) < 1e-16 ,"vm etter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<PostExp2Synapse*>(s1)->e() - e) < 1e-16 ,"e setter and getter test");

	delete s1;
	delete s2;
	delete s3;

	BOOST_MESSAGE("PostExp2Synapse working.");
	cout<<"PostExp2Synapse working."<<endl;
}

/// end of test suite postsynapses
BOOST_AUTO_TEST_SUITE_END();
