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
#include "../../split_synapse_handler/IPreSynapse.h"
#include "../../split_synapse_handler/AlphaPreSynapse.h"
#include "../../split_synapse_handler/Exp2PreSynapse.h"
#include "../../split_synapse_handler/IPostSynapse.h"
#include "../../split_synapse_handler/AlphaPostSynapse.h"
#include "../../split_synapse_handler/Exp2PostSynapse.h"
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

	IPreSynapse *s1 = new AlphaPreSynapse(0.0, 0.0);
	IPreSynapse *s2 = new AlphaPreSynapse(1, 1, 0.0, 0.0);
	IPreSynapse *s3 = new AlphaPreSynapse(3, 2, 0.0, 0.0);

	BOOST_REQUIRE_MESSAGE(s1->type() == ALPHA_PRE_SYNAPSE,"s1's type is PRE_ALPHA_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s1->name() == "ALPHA_PRE_SYNAPSE","s1 is a PRE_ALPHA_SYNAPSE");

	s1->set_id(0);
	s1->set_postsynapse_id(3);
	BOOST_REQUIRE_MESSAGE(s1->id() == 0,"id check");
	BOOST_REQUIRE_MESSAGE(s2->id() == 1,"id check");
	BOOST_REQUIRE_MESSAGE(s3->id() == 3,"id check");
	BOOST_REQUIRE_MESSAGE(s3->postsynapse_id() == 2,"postsynapse_id check");

	s1->set_location(location);
	BOOST_REQUIRE_MESSAGE( fabs(s1->location() - location) < 1e-16,"location setter and getter test"); //"perils of floating point comparisons"

	static_cast<AlphaPreSynapse*>(s1)->set_onset(onset);
	BOOST_REQUIRE_MESSAGE( fabs(static_cast<AlphaPreSynapse*>(s1)->onset() -onset) < 1e-16 ,"onset setter and getter test");

	BOOST_REQUIRE_MESSAGE(s1->is_active(t1) == false ,"synapse should not be active");
	BOOST_REQUIRE_MESSAGE(s1->is_active(t2) == true ,"synapse should be active");

	//todo:
	//test update() somehow


	delete s1;
	delete s2;
	delete s3;

	BOOST_MESSAGE("AlphaPreSynapse working.");
	cout<<"AlphaPreSynapse working."<<endl;
}

BOOST_AUTO_TEST_CASE(EXP2SYNAPSE) {
	number location = 0.3;
	number onset = 1e-6;
	number t1 = 1e-7;
	number t2 = 3e-6;

	IPreSynapse *s1 = new Exp2PreSynapse(0.0, 0.0);
	IPreSynapse *s2 = new Exp2PreSynapse(0.0, 0.0);
	IPreSynapse *s3 = new Exp2PreSynapse(3, 5, 0.0, 0.0);

	s2->set_id(2);
	s2->set_postsynapse_id(42);
	BOOST_REQUIRE_MESSAGE(s1->id() == 0,"id check");
	BOOST_REQUIRE_MESSAGE(s2->id() == 2,"id check");
	BOOST_REQUIRE_MESSAGE(s2->postsynapse_id() == 42,"postsynapse_id check");
	BOOST_REQUIRE_MESSAGE(s3->id() == 3,"id check");
	BOOST_REQUIRE_MESSAGE(s3->postsynapse_id() == 5,"id check");

	s1->set_location(location);
	BOOST_REQUIRE_MESSAGE( fabs(s1->location() - location) < 1e-16,"location setter and getter test"); //"perils of floating point comparisons"

	BOOST_REQUIRE_MESSAGE(s2->type() == EXP2_PRE_SYNAPSE,"type check");
	BOOST_REQUIRE_MESSAGE(s1->name() == "EXP2_PRE_SYNAPSE","name check");

	static_cast<Exp2PreSynapse*>(s1)->set_onset(onset);
	BOOST_REQUIRE_MESSAGE( fabs(static_cast<Exp2PreSynapse*>(s1)->onset() -onset) < 1e-16 ,"onset setter and getter test");

	BOOST_REQUIRE_MESSAGE(s1->is_active(t1) == false ,"synapse should not be active");
	BOOST_REQUIRE_MESSAGE(s1->is_active(t2) == true ,"synapse should be active");

	//todo:
	//test update() somehow

	delete s1;
	delete s2;
	delete s3;

	BOOST_MESSAGE("Exp2PreSynapse working.");
	cout<<"Exp2PreSynapse working."<<endl;
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

	IPostSynapse *s1 = new AlphaPostSynapse(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	IPostSynapse *s2 = new AlphaPostSynapse(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	IPostSynapse *s3 = new AlphaPostSynapse(1, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

	BOOST_REQUIRE_MESSAGE(s1->name() == "POST_ALPHA_SYNAPSE","s1 is a POST_ALPHA_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s1->type() == ALPHA_POST_SYNAPSE,"s1's type is POST_ALPHA_SYNAPSE");

	s1->set_id(3);
	s1->set_presynapse_id(4);
	s1->set_location(location);
	static_cast<AlphaPostSynapse*>(s1)->set_gMax(gmax);
	static_cast<AlphaPostSynapse*>(s1)->set_onset(onset);
	static_cast<AlphaPostSynapse*>(s1)->set_tau(tau);
	static_cast<AlphaPostSynapse*>(s1)->set_vm(vm);
	static_cast<AlphaPostSynapse*>(s1)->set_e(e);

	BOOST_REQUIRE_MESSAGE(s1->id() == 3, "id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s2->id() == 0, "id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s3->id() == 1, "id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s1->presynapse_id() == 4, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s2->presynapse_id() == 0, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s3->presynapse_id() == 2, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(fabs(s1->location() - location) < 1e-16, "location setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<AlphaPostSynapse*>(s1)->gMax() - gmax) < 1e-16 ,"gmax setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<AlphaPostSynapse*>(s1)->onset() - onset) < 1e-16 ,"onset setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<AlphaPostSynapse*>(s1)->tau() - tau) < 1e-16 ,"tau setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<AlphaPostSynapse*>(s1)->vm() - vm) < 1e-16 ,"vm etter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<AlphaPostSynapse*>(s1)->e() - e) < 1e-16 ,"e setter and getter test");

	//todo:
	//test current?

	delete s1;
	delete s2;
	delete s3;

	BOOST_MESSAGE("AlphaPostSynapse working.");
	cout<<"AlphaPostSynapse working."<<endl;
}

BOOST_AUTO_TEST_CASE(EXP2SYNAPSE) {
	number tau1 = 2;
	number tau2 = 3;
	number w = 1;
	number vm = 4;
	number e = 5;
	number location = 6;

	IPostSynapse *s1 = new Exp2PostSynapse(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	IPostSynapse *s2 = new Exp2PostSynapse(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	IPostSynapse *s3 = new Exp2PostSynapse(5, 6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

	s1->set_id(32);
	s1->set_presynapse_id(42);
	BOOST_REQUIRE_MESSAGE(s1->id() == 32, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s1->presynapse_id() == 42, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s2->id() == 0, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s2->presynapse_id() == 0, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s3->id() == 5, "presynapse id getter and setter test");
	BOOST_REQUIRE_MESSAGE(s3->presynapse_id() == 6, "presynapse id getter and setter test");

	BOOST_REQUIRE_MESSAGE(s1->name() == "POST_EXP2_SYNAPSE","s1 is a POST_ALPHA_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s1->type() == EXP2_POST_SYNAPSE,"s1's type is POST_ALPHA_SYNAPSE");

	static_cast<Exp2PostSynapse*>(s1)->set_tau1(tau1);
	static_cast<Exp2PostSynapse*>(s1)->set_w(w);
	static_cast<Exp2PostSynapse*>(s1)->set_tau2(tau2);
	static_cast<Exp2PostSynapse*>(s1)->set_vm(vm);
	static_cast<Exp2PostSynapse*>(s1)->set_e(e);
	s1->set_location(location);
	BOOST_REQUIRE_MESSAGE(fabs(s1->location() - location) < 1e-16, "location setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<Exp2PostSynapse*>(s1)->tau1() - tau1) < 1e-16 ,"gmax setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<Exp2PostSynapse*>(s1)->tau2() - tau2) < 1e-16 ,"onset setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<Exp2PostSynapse*>(s1)->w() - w) < 1e-16 ,"tau setter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<Exp2PostSynapse*>(s1)->vm() - vm) < 1e-16 ,"vm etter and getter test");
	BOOST_REQUIRE_MESSAGE(fabs(static_cast<Exp2PostSynapse*>(s1)->e() - e) < 1e-16 ,"e setter and getter test");

	delete s1;
	delete s2;
	delete s3;

	BOOST_MESSAGE("Exp2PostSynapse working.");
	cout<<"Exp2PostSynapse working."<<endl;
}

/// end of test suite postsynapses
BOOST_AUTO_TEST_SUITE_END();
