/*!
 * \file plugins/experimental/synapse_handler/unit_tests/src/tests.cpp
 * \brief boost unit test for the synapse_handler plugin
 *
 * \author Stephan Grein
 */

#define BOOST_TEST_MODULE __CPP__UNIT_TESTS__UG__SYNAPSE_HANDLER__

//#include <sstream>

// boost includes
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

// ug includes
#include "../../split_synapse_handler/pre_synapse.h"
#include "../../split_synapse_handler/alpha_pre_synapse.h"
#include "../../split_synapse_handler/exp2_pre_synapse.h"
#include "../../split_synapse_handler/post_synapse.h"
#include "../../split_synapse_handler/alpha_post_synapse.h"
#include "../../split_synapse_handler/exp2_post_synapse.h"
#include "../../synapse_handler/function/types.h"
#include "../../split_synapse_handler/synapse_dealer.h"
#include "../../split_synapse_distributor/split_synapse_distributor.h"



// using directives
using namespace boost::unit_test;
using namespace ug::cable_neuron;
using namespace synapse_handler;
using namespace std;

/// begin of test suite presynapses
BOOST_AUTO_TEST_SUITE(PRESYNAPSES);

BOOST_AUTO_TEST_CASE(ALPHAPRESYNAPSE) {
	SYNAPSE_ID id0 = 1;
	SYNAPSE_ID id1 = 10;
	SYNAPSE_ID id2 = 23;

	number location0 = 0;
	number location1 = 1;
	number location2 = 0.4;

	number onset0 = 1e-5;
	number onset1 = 3e-5;
	number onset2 = 6e-5;

	number duration0 = 8e-5;
	number duration1 = 8e-5;
	number duration2 = 8e-5;

	number t0 = 0.0;
	number t1 = 2e-5;
	number t2 = 5e-5;
	number t3 = 8e-5;
	number t4 = 12e-5;

	AlphaPreSynapse* s = new AlphaPreSynapse();
	AlphaPreSynapse* u = new AlphaPreSynapse(location1, onset1, duration1);

	AlphaPreSynapse* s0 = new AlphaPreSynapse(id0, id0, location0, onset0, duration0);
	AlphaPreSynapse* s1 = new AlphaPreSynapse(id1, id1, location1, onset1, duration1);
	AlphaPreSynapse* s2 = new AlphaPreSynapse(id2, id2, location2, onset2, duration2);
	AlphaPreSynapse* s3 = new AlphaPreSynapse();

	/**
	 * Getter and setter tests
	 */
	s->set_id(id1);
	s->set_location(location2);
	s->set_onset(onset1);
	s->set_duration(duration0);
	BOOST_REQUIRE_MESSAGE(s->id() == id1,"id check0");
	BOOST_REQUIRE_MESSAGE(s->location() == location2,"location check0");
	BOOST_REQUIRE_MESSAGE(s->onset() == onset1,"onset check0");
	BOOST_REQUIRE_MESSAGE(s->duration() == duration0,"duration check0");
	BOOST_REQUIRE_MESSAGE(s->type() == ALPHA_PRE_SYNAPSE,"type check");

	u->set_id(id2);
	BOOST_REQUIRE_MESSAGE(u->id() == id2,"id check1");
	BOOST_REQUIRE_MESSAGE(u->name() == "ALPHA_PRE_SYNAPSE", "name check");


	/**
	 * Functionality tests
	 */
	BOOST_REQUIRE_MESSAGE( !s0->is_active(t0) &&
						   !s1->is_active(t0) &&
						   !s2->is_active(t0), "active check (t0)");

	BOOST_REQUIRE_MESSAGE(  s0->is_active(t1) &&
						   !s1->is_active(t1) &&
						   !s2->is_active(t1), "active check (t1)");

	BOOST_REQUIRE_MESSAGE(  s0->is_active(t2) &&
						    s1->is_active(t2) &&
						   !s2->is_active(t2), "active check (t2)");

	BOOST_REQUIRE_MESSAGE(  s0->is_active(t3) &&
						    s1->is_active(t3) &&
						    s2->is_active(t3), "active check (t3)");

	BOOST_REQUIRE_MESSAGE( !s0->is_active(t4) &&
						   !s1->is_active(t4) &&
						    s2->is_active(t4), "active check (t4)");

	//serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("ALPHA_PRE_SYNAPSE 111 10 1.4 3e-08 4");
	oss1 << s1;

	oss2 << "ALPHA_PRE_SYNAPSE" << " ";
	oss2 << id1 << " ";
	oss2 << id1 << " ";
	oss2 << location1 << " ";
	oss2 << onset1 << " ";
	oss2 << duration1;

	BOOST_REQUIRE_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; //pop identifier away
	iss >> s3;
	BOOST_REQUIRE_MESSAGE(s3->type() == ALPHA_PRE_SYNAPSE, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->name() == "ALPHA_PRE_SYNAPSE", "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->id() == 111, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->postsynapse_id() == 10, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->location() == 1.4, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->onset() == 3e-08, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->duration() == 4, "deserialization check");

//	synapse dealer
	SynapseDealer::instance()->register_synapse_type<AlphaPreSynapse>("ALPHA_PRE_SYNAPSE");
	IBaseSynapse* s4 = SynapseDealer::instance()->deal("ALPHA_PRE_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s4->name() == "ALPHA_PRE_SYNAPSE", "SynapseDealer check");

	delete s;
	delete u;
	delete s0;
	delete s1;
	delete s2;
	delete s3;
	delete s4;
}

BOOST_AUTO_TEST_CASE(EXP2PRESYNAPSE) {

	SYNAPSE_ID id0 = 0;
	SYNAPSE_ID id1 = 1;
	SYNAPSE_ID id2 = 4;

	number location0 = 0;
	number location1 = 1;
	number location2 = 0.3;

	number onset0 = 0;
	number onset1 = 1e-5;
	number onset2 = 5e-5;

	number duration0 = 8e-5;
	number duration1 = 8e-5;
	number duration2 = 8e-5;

	number threshold0 = -0.01;
	number threshold1 = -0.01;
	number threshold2 = -0.01;

	number t0 = 1e-5;
	number t1 = 3e-5;
	number t2 = 1e-4;

	number v0 = -0.065; //-65mV
	number v1 = -0.020; //-20mV
	number v2 = -0.008; //-8mV
	std::vector<number> x0; x0.push_back(v0); //unknown vectors
	std::vector<number> x1; x1.push_back(v1);
	std::vector<number> x2; x2.push_back(v2);

	/**
	 * Test getter and setter with these
	 */
	Exp2PreSynapse* s = new Exp2PreSynapse();
	Exp2PreSynapse* u = new Exp2PreSynapse(location2, onset2, duration2, threshold2);


	/**
	 * Test current/activation/deactivation functionality with these
	 */
	Exp2PreSynapse* s0 = new Exp2PreSynapse(id0, id0, location0, onset0, duration0, threshold0);
	Exp2PreSynapse* s1 = new Exp2PreSynapse(id1, id1, location1, onset1, duration1, threshold1);
	Exp2PreSynapse* s2 = new Exp2PreSynapse(id2, id2, location2, onset2, duration2, threshold2);
	Exp2PreSynapse* s3 = new Exp2PreSynapse();


	/**
	 * Getter and setter tests
	 */
	s->set_id(id1);
	s->set_location(location1);
	s->set_onset(onset1);
	s->set_duration(duration1);
	s->set_threshold(threshold1);

	u->set_id(id2);

	BOOST_REQUIRE_MESSAGE(s->id() == id1,"id getter and setter check1");
	BOOST_REQUIRE_MESSAGE(s->location() == location1,"location getter and setter check1");
	BOOST_REQUIRE_MESSAGE(s->onset() == onset1,"onset getter and setter check1");
	BOOST_REQUIRE_MESSAGE(s->duration() == duration1,"duration getter and setter check1");
	BOOST_REQUIRE_MESSAGE(s->threshold() == threshold1,"threshold getter and setter check1");
	BOOST_REQUIRE_MESSAGE(s->type() == EXP2_PRE_SYNAPSE,"type check");

	BOOST_REQUIRE_MESSAGE(u->id() == id2,"id getter and setter check2");
	BOOST_REQUIRE_MESSAGE(u->name() == "EXP2_PRE_SYNAPSE","name check");

	/**
	 * Functionality tests
	 */

	//rising potential until -13mV
	{
		number vm = v0;
		number t = 0.0;
		number dt = 1e-5;
		number dv = 1e-3;

		while(vm <= -1.3e-2) {
			std::vector<number> x; x.push_back(vm);
			s0->update(t, x);

			BOOST_REQUIRE_MESSAGE((s0->is_active(t)) == false, "active check0(t="<<t<<")");

			t += dt;
			vm += dv;
		}
	}


	//rising potential until -7mV
	{
		number vm = v0;
		number t = 0.0;
		number dt = 1e-5;
		number dv = 1e-3;

		while(vm <= -7e-3) {
			std::vector<number> x; x.push_back(vm);
			s0->update(t, x);

			BOOST_REQUIRE_MESSAGE((s0->is_active(t)) == (vm >= threshold0), "active check1");

			t += dt;
			vm += dv;
		}
	}



	//serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("EXP2_PRE_SYNAPSE 111 10 1.4 3e-08 3 -0.010");
	oss1 << s2;

	oss2 << "EXP2_PRE_SYNAPSE" << " ";
	oss2 << id2 << " ";
	oss2 << id2 << " ";
	oss2 << location2 << " ";
	oss2 << nan("") << " ";
	oss2 << duration2 << " ";
	oss2 << threshold2;
	BOOST_REQUIRE_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; //pop identifier away
	iss >> s3;
	BOOST_REQUIRE_MESSAGE(s3->type() == EXP2_PRE_SYNAPSE, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->name() == "EXP2_PRE_SYNAPSE", "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->id() == 111, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->postsynapse_id() == 10, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->location() == 1.4, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->onset() != s3->onset(), "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->duration() == 3, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->threshold() == -0.010, "deserialization check");

//	synapse dealer
	SynapseDealer::instance()->register_synapse_type<Exp2PreSynapse>("EXP2_PRE_SYNAPSE");
	IBaseSynapse* s4 = SynapseDealer::instance()->deal("EXP2_PRE_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s4->name() == "EXP2_PRE_SYNAPSE", "SynapseDealer check");

	delete s;
	delete u;
	delete s0;
	delete s1;
	delete s2;
	delete s3;
	delete s4;
}

/// end of test suite presynapses
BOOST_AUTO_TEST_SUITE_END();



/// begin of test suite postsynapses
BOOST_AUTO_TEST_SUITE(POSTSYNAPSES);

BOOST_AUTO_TEST_CASE(ALPHAPOSTSYNAPSE) {
	SYNAPSE_ID id0 = 1;
	SYNAPSE_ID id1 = 3;
	SYNAPSE_ID id2 = 195;

	number location0 = 0;
	number location1 = 1;
	number location2 = 0.6;

	number onset0 = 2;
	number onset1 = 2;
	number onset2 = 2;

	number gmax0 = 0.000511591;
	number gmax1 = 0.000511591;
	number gmax2 = 0.000511591;

	number tau0 = 1.7;
	number tau1 = 1.7;
	number tau2 = 1.7;

	number rev0 = 0.0;
	number rev1 = 0.0;
	number rev2 = 0.0;

	number vm0 = -0.065;
	number vm1 = -0.065;
	number vm2 = -0.065;

	number t0 = 1e-5;
	number t1 = 5e-5;
	number t2 = 7e-5;
	number t3 = 9e-5;


	AlphaPostSynapse* s = new AlphaPostSynapse();
	s->set_id(id2);
	s->set_location(location2);
	s->set_onset(onset2);
	s->set_gMax(gmax2);
	s->set_tau(tau2);
	s->set_rev(rev2);
	BOOST_REQUIRE_MESSAGE(s->name() == "ALPHA_POST_SYNAPSE","name check");
	BOOST_REQUIRE_MESSAGE(s->id() == id2,"getter and setter check0");
	BOOST_REQUIRE_MESSAGE(s->location() == location2,"getter and setter check0");
	BOOST_REQUIRE_MESSAGE(s->onset() == onset2,"getter and setter check0");
	BOOST_REQUIRE_MESSAGE(s->gMax() == gmax2,"getter and setter check0");
	BOOST_REQUIRE_MESSAGE(s->tau() == tau2,"getter and setter check0");
	BOOST_REQUIRE_MESSAGE(s->rev() == rev2,"getter and setter check0");
	BOOST_REQUIRE_MESSAGE(s->type() == ALPHA_POST_SYNAPSE,"type check");

	AlphaPostSynapse* s1 = new AlphaPostSynapse(id0, id0, location0, onset0, gmax0, tau0, rev0);
	AlphaPostSynapse* s2 = new AlphaPostSynapse(id1, id1, location1, onset1, gmax1, tau1, rev1);
	AlphaPostSynapse* s3 = new AlphaPostSynapse(id2, id2, location2, onset2, gmax2, tau2, rev2);
	AlphaPostSynapse* s4 = new AlphaPostSynapse();

	BOOST_REQUIRE_MESSAGE(s1->is_active(t0) == false,"is_active check");

	s2->activate(t1);
	BOOST_REQUIRE_MESSAGE(s2->is_active(t1) == true,"is_active check");

	s3->activate(t2);
	s3->deactivate(t3);
	BOOST_REQUIRE_MESSAGE(s3->is_active(t3) == false,"is_active check");


	//todo:
	//test current?

	//serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("ALPHA_POST_SYNAPSE 1 2 4 3 1 2 1");
	oss1 << s3;
	oss2 << "ALPHA_POST_SYNAPSE" << " ";
	oss2 << id2 << " ";
	oss2 << id2 << " ";
	oss2 << location2 << " ";
	oss2 << nan("") << " ";
	oss2 << gmax2 << " ";
	oss2 << tau2 << " ";
	oss2 << rev2;
	BOOST_REQUIRE_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; //pop identifier away

	iss >> s4;
	BOOST_REQUIRE_MESSAGE(s4->id()== 1, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s4->presynapse_id() == 2, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s4->location() == 4, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s4->onset() != s4->onset(), "deserialization check");
	BOOST_REQUIRE_MESSAGE(s4->gMax() == 1, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s4->tau() == 2, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s4->rev() == 1, "deserialization check");

	//synapse dealer
	SynapseDealer::instance()->register_synapse_type<AlphaPostSynapse>("ALPHA_POST_SYNAPSE");
	IBaseSynapse* s0 = SynapseDealer::instance()->deal("ALPHA_POST_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s4->name() == "ALPHA_POST_SYNAPSE", "SynapseDealer check");

	delete s;
	delete s0;
	delete s1;
	delete s2;
	delete s3;
	delete s4;
}

BOOST_AUTO_TEST_CASE(EXP2POSTSYNAPSE) {
	SYNAPSE_ID id0 = 1;
	SYNAPSE_ID id1 = 3;
	SYNAPSE_ID id2 = 195;

	number location0 = 0;
	number location1 = 1;
	number location2 = 0.6;

	number onset0 = 2;
	number onset1 = 2;
	number onset2 = 2;

	number gmax0 = 0.000511591;
	number gmax1 = 0.000511591;
	number gmax2 = 0.000511591;

	number tau10 = 1.7;
	number tau11 = 1.7;
	number tau12 = 1.7;

	number tau20 = 0.2;
	number tau21 = 0.2;
	number tau22 = 0.2;

	number rev0 = 0.0;
	number rev1 = 0.0;
	number rev2 = 0.0;

	number vm0 = -0.065;
	number vm1 = -0.065;
	number vm2 = -0.065;

	number t0 = 1e-5;
	number t1 = 5e-5;
	number t2 = 7e-5;
	number t3 = 9e-5;

	Exp2PostSynapse* s = new Exp2PostSynapse();
	s->set_id(id2);
	s->set_location(location2);
	s->set_onset(onset2);
	s->set_gMax(gmax2);
	s->set_tau1(tau12);
	s->set_tau2(tau22);
	s->set_rev(rev2);
	BOOST_REQUIRE_MESSAGE(s->name() == "EXP2_POST_SYNAPSE", "name test");
	BOOST_REQUIRE_MESSAGE(s->type() == EXP2_POST_SYNAPSE, "type test");
	BOOST_REQUIRE_MESSAGE(s->id() == id2, "getter and setter test");
	BOOST_REQUIRE_MESSAGE(s->location() == location2, "getter and setter test");
	BOOST_REQUIRE_MESSAGE(s->onset() == onset2, "getter and setter test");
	BOOST_REQUIRE_MESSAGE(s->gMax() == gmax2, "getter and setter test");
	BOOST_REQUIRE_MESSAGE(s->tau1() == tau12, "getter and setter test");
	BOOST_REQUIRE_MESSAGE(s->tau2() == tau22, "getter and setter test");
	BOOST_REQUIRE_MESSAGE(s->rev() == rev2, "getter and setter test");

	Exp2PostSynapse* s0 = new Exp2PostSynapse(id0, id0, location0, onset0, gmax0, tau10, tau20, rev0);
	Exp2PostSynapse* s1 = new Exp2PostSynapse(id1, id1, location1, onset1, gmax1, tau11, tau21, rev1);
	Exp2PostSynapse* s2 = new Exp2PostSynapse(id2, id2, location2, onset2, gmax2, tau12, tau22, rev2);
	Exp2PostSynapse* s3 = new Exp2PostSynapse();


	BOOST_REQUIRE_MESSAGE(s1->is_active(t0) == false,"is_active check");

	s2->activate(t1);
	BOOST_REQUIRE_MESSAGE(s2->is_active(t1) == true,"is_active check");

	s3->activate(t2);
	s3->deactivate(t3);
	BOOST_REQUIRE_MESSAGE(s3->is_active(t3) == false,"is_active check");

	//serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("EXP2_POST_SYNAPSE 1 22 3 1 1 2 4 1");
	oss1 << s1;
	oss2 << "EXP2_POST_SYNAPSE" << " ";
	oss2 << id1 << " ";
	oss2 << id1 << " ";
	oss2 << location1 << " ";
	oss2 << nan("") << " ";
	oss2 << gmax1 << " ";
	oss2 << tau11 << " ";
	oss2 << tau21 << " ";
	oss2 << rev1;
	BOOST_REQUIRE_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; //pop identifier away
	iss >> s3;
	BOOST_REQUIRE_MESSAGE(s3->id() == 1, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->presynapse_id() == 22, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->location() == 3, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->onset() != s3->onset(), "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->gMax() == 1, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->tau1() == 2, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->tau2() == 4, "deserialization check");
	BOOST_REQUIRE_MESSAGE(s3->rev() == 1, "deserialization check");

	//synapse dealer
	SynapseDealer::instance()->register_synapse_type<Exp2PostSynapse>("EXP2_POST_SYNAPSE");
	IBaseSynapse* s4 = SynapseDealer::instance()->deal("EXP2_POST_SYNAPSE");
	BOOST_REQUIRE_MESSAGE(s4->name() == "EXP2_POST_SYNAPSE", "SynapseDealer check");

	delete s;
	delete s0;
	delete s1;
	delete s2;
	delete s3;
}

/// end of test suite postsynapses
BOOST_AUTO_TEST_SUITE_END();
