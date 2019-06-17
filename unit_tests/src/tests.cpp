/*
 * Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 *
 * Authors: Markus Breit, Stephan Grein, Lukas Reinhardt
 * Creation date: 2019-06-12
 *
 * This file is part of NeuroBox, which is based on UG4.
 *
 * NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3
 * (as published by the Free Software Foundation) with the following additional
 * attribution requirements (according to LGPL/GPL v3 §7):
 *
 * (1) The following notice must be displayed in the appropriate legal notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 *
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 *
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating PDE based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 *   NeuroBox: computational mathematics in multiscale neuroscience.
 *   Computing and visualization in science (2019).
 * "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 *   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*!
 * \brief boost unit test for the synapse_handler plugin
 *
 * \author Stephan Grein
 */

#define BOOST_TEST_MODULE __CPP__UNIT_TESTS__UG__SYNAPSE_HANDLING__

//#include <sstream>

// boost includes
#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

// ug includes
#include "../../synapse_handling/synapses/pre_synapse.h"
#include "../../synapse_handling/synapses/onset_pre_synapse.h"
#include "../../synapse_handling/synapses/threshold_pre_synapse.h"
#include "../../synapse_handling/synapses/post_synapse.h"
#include "../../synapse_handling/synapses/alpha_post_synapse.h"
#include "../../synapse_handling/synapses/exp2_post_synapse.h"
#include "../../synapse_handling/types.h"
#include "../../synapse_handling/synapse_dealer.h"
#include "../../synapse_handling/synapse_distributor.h"
#include "../../util/neuronal_topology_importer.h"
#include "geometry_information.h"

#include "ug.h"  // UGInit
#include "bridge/bridge.h"  // GetUGRegistry
#include "lib_disc/domain.h"  // Domain3d



// using directives
using namespace boost::unit_test;
using namespace ug;
using namespace ug::cable_neuron;
using namespace ug::cable_neuron::synapse_handler;
using namespace ug::cable_neuron::neuronal_topology_importer;
using namespace std;

/// begin of test suite presynapses
BOOST_AUTO_TEST_SUITE(PRESYNAPSES);

BOOST_AUTO_TEST_CASE(ALPHAPRESYNAPSE) {
	// init ug (only needed once in very first test)
	try {
		int argc = 1;
		char arg0[8] = {'u', 'g', 's', 'h', 'e', 'l', 'l', '\0'};
		char* arg0p = &arg0[0];
		char** argv = &arg0p;
		UGInit(&argc, &argv);
		script::RegisterDefaultLuaBridge(&bridge::GetUGRegistry());
	} catch (const UGError& e) {
		throw framework::setup_error(e.get_stacktrace());
	}

	synapse_id id0 = 1;
	synapse_id id1 = 10;
	synapse_id id2 = 23;

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

	OnsetPreSynapse* s = new OnsetPreSynapse();
	OnsetPreSynapse* u = new OnsetPreSynapse(location1, onset1, duration1);

	OnsetPreSynapse* s0 = new OnsetPreSynapse(id0, location0, onset0, duration0);
	OnsetPreSynapse* s1 = new OnsetPreSynapse(id1, location1, onset1, duration1);
	OnsetPreSynapse* s2 = new OnsetPreSynapse(id2, location2, onset2, duration2);
	OnsetPreSynapse* s3 = new OnsetPreSynapse();

	// getter and setter tests
	s->set_id(id1);
	s->set_location(location2);
	s->set_onset(onset1);
	s->set_duration(duration0);
	BOOST_CHECK_MESSAGE(s->id() == id1,"id check0");
	BOOST_CHECK_MESSAGE(s->location() == location2,"location check0");
	BOOST_CHECK_MESSAGE(s->onset() == onset1,"onset check0");
	BOOST_CHECK_MESSAGE(s->duration() == duration0,"duration check0");
	BOOST_CHECK_MESSAGE(s->type() == ONSET_PRE_SYNAPSE,"type check");

	u->set_id(id2);
	BOOST_CHECK_MESSAGE(u->id() == id2,"id check1");
	BOOST_CHECK_MESSAGE(u->name() == OnsetPreSynapse::name_string, "name check");


	// functionality tests
	BOOST_CHECK_MESSAGE( !s0->is_active(t0) &&
						   !s1->is_active(t0) &&
						   !s2->is_active(t0), "active check (t0)");

	BOOST_CHECK_MESSAGE(  s0->is_active(t1) &&
						   !s1->is_active(t1) &&
						   !s2->is_active(t1), "active check (t1)");

	BOOST_CHECK_MESSAGE(  s0->is_active(t2) &&
						    s1->is_active(t2) &&
						   !s2->is_active(t2), "active check (t2)");

	BOOST_CHECK_MESSAGE(  s0->is_active(t3) &&
						    s1->is_active(t3) &&
						    s2->is_active(t3), "active check (t3)");

	BOOST_CHECK_MESSAGE( !s0->is_active(t4) &&
						   !s1->is_active(t4) &&
						    s2->is_active(t4), "active check (t4)");

	//serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("OnsetPreSynapse 111 1.4 3e-08 4");
	oss1 << s1;

	oss2 << "OnsetPreSynapse" << " ";
	oss2 << id1 << " ";
	oss2 << location1 << " ";
	oss2 << onset1 << " ";
	oss2 << duration1;

	BOOST_CHECK_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; //pop identifier away
	iss >> s3;
	BOOST_CHECK_MESSAGE(s3->type() == ONSET_PRE_SYNAPSE, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->name() == "OnsetPreSynapse", "deserialization check");
	BOOST_CHECK_MESSAGE(s3->id() == 111, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->location() == 1.4, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->onset() == 3e-08, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->duration() == 4, "deserialization check");

	// synapse dealer
	SynapseDealer::instance()->register_synapse_type<OnsetPreSynapse>();
	IBaseSynapse* s4 = SynapseDealer::instance()->deal(OnsetPreSynapse::name_string);
	BOOST_CHECK_MESSAGE(s4->name() == OnsetPreSynapse::name_string, "SynapseDealer check");

	delete s;
	delete u;
	delete s0;
	delete s1;
	delete s2;
	delete s3;
	delete s4;
}

BOOST_AUTO_TEST_CASE(EXP2PRESYNAPSE) {

	synapse_id id0 = 0;
	synapse_id id1 = 1;
	synapse_id id2 = 4;

	number location0 = 0;
	number location1 = 1;
	number location2 = 0.3;

	number onset1 = 1e-5;
	number onset2 = 5e-5;

	number duration0 = 8e-5;
	number duration1 = 8e-5;
	number duration2 = 8e-5;

	number threshold0 = -0.01;
	number threshold1 = -0.01;
	number threshold2 = -0.01;

	number v0 = -0.065; //-65mV
	number v1 = -0.020; //-20mV
	number v2 = -0.008; //-8mV
	std::vector<number> x0; x0.push_back(v0); //unknown vectors
	std::vector<number> x1; x1.push_back(v1);
	std::vector<number> x2; x2.push_back(v2);

	// test getter and setter with these
	ThresholdPreSynapse* s = new ThresholdPreSynapse();
	ThresholdPreSynapse* u = new ThresholdPreSynapse(location2, onset2, duration2, threshold2);

	// test current/activation/deactivation functionality with these
	ThresholdPreSynapse* s0 = new ThresholdPreSynapse(id0, location0, duration0, threshold0);
	ThresholdPreSynapse* s1 = new ThresholdPreSynapse(id1, location1, duration1, threshold1);
	ThresholdPreSynapse* s2 = new ThresholdPreSynapse(id2, location2, duration2, threshold2);
	ThresholdPreSynapse* s3 = new ThresholdPreSynapse();


	// getter and setter tests
	s->set_id(id1);
	s->set_location(location1);
	s->set_onset(onset1);
	s->set_duration(duration1);
	s->set_threshold(threshold1);

	u->set_id(id2);

	BOOST_CHECK_MESSAGE(s->id() == id1,"id getter and setter check1");
	BOOST_CHECK_MESSAGE(s->location() == location1,"location getter and setter check1");
	BOOST_CHECK_MESSAGE(s->onset() == onset1,"onset getter and setter check1");
	BOOST_CHECK_MESSAGE(s->duration() == duration1,"duration getter and setter check1");
	BOOST_CHECK_MESSAGE(s->threshold() == threshold1,"threshold getter and setter check1");
	BOOST_CHECK_MESSAGE(s->type() == THRESHOLD_PRE_SYNAPSE,"type check");

	BOOST_CHECK_MESSAGE(u->id() == id2,"id getter and setter check2");
	BOOST_CHECK_MESSAGE(u->name() == ThresholdPreSynapse::name_string,"name check");


	// functionality tests

	// rising potential until -13mV
	{
		number vm = v0;
		number t = 0.0;
		number dt = 1e-5;
		number dv = 1e-3;

		while(vm <= -1.3e-2) {
			std::vector<number> x; x.push_back(vm);
			s0->update(t, x);

			BOOST_CHECK_MESSAGE((s0->is_active(t)) == false, "active check0(t="<<t<<")");

			t += dt;
			vm += dv;
		}
	}


	// rising potential until -7mV
	{
		number vm = v0;
		number t = 0.0;
		number dt = 1e-5;
		number dv = 1e-3;

		while(vm <= -7e-3) {
			std::vector<number> x; x.push_back(vm);
			s0->update(t, x);

			BOOST_CHECK_MESSAGE((s0->is_active(t)) == (vm >= threshold0), "active check1");

			t += dt;
			vm += dv;
		}
	}



	// serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("ThresholdPreSynapse 111 1.4 3 -0.010");
	oss1 << s2;

	oss2 << "ThresholdPreSynapse" << " ";
	oss2 << id2 << " ";
	oss2 << location2 << " ";
	oss2 << duration2 << " ";
	oss2 << threshold2;
	BOOST_CHECK_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; // pop identifier away
	iss >> s3;
	BOOST_CHECK_MESSAGE(s3->type() == THRESHOLD_PRE_SYNAPSE, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->name() == "ThresholdPreSynapse", "deserialization check");
	BOOST_CHECK_MESSAGE(s3->id() == 111, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->location() == 1.4, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->duration() == 3, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->threshold() == -0.010, "deserialization check");

	// synapse dealer
	SynapseDealer::instance()->register_synapse_type<ThresholdPreSynapse>();
	IBaseSynapse* s4 = SynapseDealer::instance()->deal(ThresholdPreSynapse::name_string);
	BOOST_CHECK_MESSAGE(s4->name() == ThresholdPreSynapse::name_string, "SynapseDealer check");

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
	synapse_id id2 = 195;

	number location0 = 0;
	number location1 = 1;
	number location2 = 0.6;

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
	BOOST_CHECK_MESSAGE(s->name() == AlphaPostSynapse::name_string,"name check");
	BOOST_CHECK_MESSAGE(s->id() == id2,"getter and setter check0");
	BOOST_CHECK_MESSAGE(s->location() == location2,"getter and setter check0");
	BOOST_CHECK_MESSAGE(s->onset() == onset2,"getter and setter check0");
	BOOST_CHECK_MESSAGE(s->gMax() == gmax2,"getter and setter check0");
	BOOST_CHECK_MESSAGE(s->tau() == tau2,"getter and setter check0");
	BOOST_CHECK_MESSAGE(s->rev() == rev2,"getter and setter check0");
	BOOST_CHECK_MESSAGE(s->type() == ALPHA_POST_SYNAPSE,"type check");

	AlphaPostSynapse* s1 = new AlphaPostSynapse(location0, gmax0, tau0, rev0);
	AlphaPostSynapse* s2 = new AlphaPostSynapse(location1, gmax1, tau1, rev1);
	AlphaPostSynapse* s3 = new AlphaPostSynapse(location2, gmax2, tau2, rev2);
	AlphaPostSynapse* s4 = new AlphaPostSynapse();

	BOOST_CHECK_MESSAGE(s1->is_active(t0) == false,"is_active check");

	s2->activate(t1);
	BOOST_CHECK_MESSAGE(s2->is_active(t1) == true,"is_active check");

	s3->activate(t2);
	s3->deactivate();
	BOOST_CHECK_MESSAGE(s3->is_active(t3) == false,"is_active check");


	// todo: test current?

	//serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("AlphaPostSynapse 1 4 1 2 1");
	oss1 << s;
	oss2 << "AlphaPostSynapse" << " ";
	oss2 << id2 << " ";
	oss2 << location2 << " ";
	oss2 << gmax2 << " ";
	oss2 << tau2 << " ";
	oss2 << rev2;
	BOOST_CHECK_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; // pop identifier away

	iss >> s4;
	BOOST_CHECK_MESSAGE(s4->id()== 1, "deserialization check");
	BOOST_CHECK_MESSAGE(s4->name() == "AlphaPostSynapse", "deserialization check");
	BOOST_CHECK_MESSAGE(s4->location() == 4, "deserialization check");
	BOOST_CHECK_MESSAGE(s4->gMax() == 1, "deserialization check");
	BOOST_CHECK_MESSAGE(s4->tau() == 2, "deserialization check");
	BOOST_CHECK_MESSAGE(s4->rev() == 1, "deserialization check");

	// synapse dealer
	SynapseDealer::instance()->register_synapse_type<AlphaPostSynapse>();
	IBaseSynapse* s0 = SynapseDealer::instance()->deal(AlphaPostSynapse::name_string);
	BOOST_CHECK_MESSAGE(s4->name() == AlphaPostSynapse::name_string, "SynapseDealer check");

	delete s;
	delete s0;
	delete s1;
	delete s2;
	delete s3;
	delete s4;
}

BOOST_AUTO_TEST_CASE(EXP2POSTSYNAPSE) {
	synapse_id id0 = 1;
	synapse_id id1 = 3;
	synapse_id id2 = 195;

	number location0 = 0;
	number location1 = 1;
	number location2 = 0.6;

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
	BOOST_CHECK_MESSAGE(s->name() == Exp2PostSynapse::name_string, "name test");
	BOOST_CHECK_MESSAGE(s->type() == EXP2_POST_SYNAPSE, "type test");
	BOOST_CHECK_MESSAGE(s->id() == id2, "getter and setter test");
	BOOST_CHECK_MESSAGE(s->location() == location2, "getter and setter test");
	BOOST_CHECK_MESSAGE(s->onset() == onset2, "getter and setter test");
	BOOST_CHECK_MESSAGE(s->gMax() == gmax2, "getter and setter test");
	BOOST_CHECK_MESSAGE(s->tau1() == tau12, "getter and setter test");
	BOOST_CHECK_MESSAGE(s->tau2() == tau22, "getter and setter test");
	BOOST_CHECK_MESSAGE(s->rev() == rev2, "getter and setter test");

	Exp2PostSynapse* s0 = new Exp2PostSynapse(id0, location0, gmax0, tau10, tau20, rev0);
	Exp2PostSynapse* s1 = new Exp2PostSynapse(id1, location1, gmax1, tau11, tau21, rev1);
	Exp2PostSynapse* s2 = new Exp2PostSynapse(id2, location2, gmax2, tau12, tau22, rev2);
	Exp2PostSynapse* s3 = new Exp2PostSynapse();


	BOOST_CHECK_MESSAGE(s1->is_active(t0) == false,"is_active check");

	s2->activate(t1);
	BOOST_CHECK_MESSAGE(s2->is_active(t1) == true,"is_active check");

	s3->activate(t2);
	s3->deactivate();
	BOOST_CHECK_MESSAGE(s3->is_active(t3) == false,"is_active check");

	// serialization
	ostringstream oss1;
	ostringstream oss2;
	istringstream iss("Exp2PostSynapse 1 3 1 2 4 1");
	oss1 << s1;
	oss2 << "Exp2PostSynapse" << " ";
	oss2 << id1 << " ";
	oss2 << location1 << " ";
	oss2 << gmax1 << " ";
	oss2 << tau11 << " ";
	oss2 << tau21 << " ";
	oss2 << rev1;
	BOOST_CHECK_MESSAGE(oss1.str() == oss2.str(), "serialization check");

	std::string ident;
	iss >> ident; // pop identifier away
	iss >> s3;
	BOOST_CHECK_MESSAGE(s3->id() == 1, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->name() == "Exp2PostSynapse", "deserialization check");
	BOOST_CHECK_MESSAGE(s3->location() == 3, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->gMax() == 1, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->tau1() == 2, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->tau2() == 4, "deserialization check");
	BOOST_CHECK_MESSAGE(s3->rev() == 1, "deserialization check");

	// synapse dealer
	SynapseDealer::instance()->register_synapse_type<Exp2PostSynapse>();
	IBaseSynapse* s4 = SynapseDealer::instance()->deal(Exp2PostSynapse::name_string);
	BOOST_CHECK_MESSAGE(s4->name() == Exp2PostSynapse::name_string, "SynapseDealer check");

	delete s;
	delete s0;
	delete s1;
	delete s2;
	delete s3;
}

/// end of test suite postsynapses
BOOST_AUTO_TEST_SUITE_END();



////////////////////////////////////////////////////////////////////////
/// testsuite to test the swc import
////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_SUITE(SWC_IMPORT);

BOOST_AUTO_TEST_CASE(CONVERT_FROM_SWC)
{
	// convert SWC geometry
	try
	{
		NeuronalTopologyImporter neti;
		BOOST_REQUIRE_MESSAGE(neti.import_and_generate_grid("02a_pyramidal2aFI.swc"),
			"Failed to import '02a_pyramidal2aFI.swc'.");
	}
	catch (const UGError& e)
	{
		throw std::runtime_error(e.get_stacktrace());
	}

	// load converted geometry
	Domain3d dom;
	try {LoadDomain(dom, "02a_pyramidal2aFI.ugx");}
	catch (const UGError& e)
	{
		throw std::runtime_error(e.get_stacktrace());
	}

	// asserts
	BOOST_CHECK_MESSAGE(dom.grid()->num_edges() == 1513,
		"Exactly 1513 edges should be available in the generated UGX grid");
	BOOST_CHECK_MESSAGE(dom.grid()->num_vertices() == 1514,
		"Exactly 1514 vertices should be available in the generated UGX grid");
}


BOOST_AUTO_TEST_CASE(COMPARE_GENERATED_WITH_REFERENCE_GRID)
{
	std::ifstream ifs1("02a_pyramidal2aFI.ugx");
	std::ifstream ifs2("02a_pyramidal2aFI_ref.ugx");

    std::istream_iterator<char> b1(ifs1), e1;
    std::istream_iterator<char> b2(ifs2), e2;

    BOOST_CHECK_EQUAL_COLLECTIONS(b1, e1, b2, e2);
}

BOOST_AUTO_TEST_SUITE_END();


/// testsuite to test the hoc import
BOOST_AUTO_TEST_SUITE(HOC_IMPORT);

BOOST_AUTO_TEST_CASE(CONVERT_FROM_HOC)
{
	// todo implement
	BOOST_CHECK_MESSAGE(true, "Dummy check.");
}

BOOST_AUTO_TEST_SUITE_END();


/// testsuite to test the txt import
BOOST_AUTO_TEST_SUITE(TXT_IMPORT);

BOOST_AUTO_TEST_CASE(CONVERT_FROM_TXT)
{
	// test geometries
	std::vector<std::pair<std::string, GeometryInformation> > testGeometries;
	testGeometries.push_back
	(
		std::pair<std::string, GeometryInformation>(
		"small_network6", GeometryInformation(42340, 42328, 8, 441, 441, 8))
	);

	typedef std::vector<std::pair<std::string, GeometryInformation> >::const_iterator VGeometryIter;
	for (VGeometryIter it = testGeometries.begin(); it != testGeometries.end(); ++it)
	{
		// import geometry
		NeuronalTopologyImporter neti;
		std::ostringstream ossFail;
		ossFail << "Failed to import '" << it->first << "'.";
		try
		{
			BOOST_REQUIRE_MESSAGE(neti.import_and_generate_grid(it->first, "txt"),
				ossFail.str());
		}
		catch (const UGError& e)
		{
			throw std::runtime_error(e.get_stacktrace());
		}

		// load converted geometry
		Domain3d dom;
		try {LoadDomain(dom, (it->first + ".ugx").c_str());}
		catch (const UGError& e)
		{
			throw std::runtime_error(e.get_stacktrace());
		}
		SmartPtr<MultiGrid> grid = dom.grid();

		// # vertices
		BOOST_CHECK_MESSAGE(dom.grid()->num_vertices() == it->second.num_vertices(),
		   "# vertices don't match (Reference #" << it->second.num_vertices()
		   << " vs. generated #" << dom.grid()->num_vertices());

		// # edges
		BOOST_CHECK_MESSAGE(dom.grid()->num_edges() == it->second.num_edges(),
		   "# edges don't match (Reference #" << it->second.num_edges()
		   << " vs. generated #" << dom.grid()->num_edges());


		// synapse handling
		typedef Attachment<std::vector<IBaseSynapse*> > AVSynapse;
		AVSynapse aSyn = GlobalAttachments::attachment<AVSynapse>("synapses");
		Grid::EdgeAttachmentAccessor<AVSynapse> aaSyn(*grid, aSyn);

		// collect all synapses from grid
		size_t nOnsetPresyn = 0;
		size_t nThresholdPresyn = 0;
		size_t nAlphaPostsyn = 0;
		size_t nExp2Postsyn = 0;
		typename geometry_traits<Edge>::const_iterator eit = grid->begin<Edge>();
		typename geometry_traits<Edge>::const_iterator eit_end = grid->end<Edge>();
		for (; eit != eit_end; ++eit)
		{
			const std::vector<IBaseSynapse*>& syns = aaSyn[*eit];
			const size_t sz = syns.size();
			for (size_t i = 0; i < sz; ++i)
			{
				IBaseSynapse* syn = syns[i];
				if (dynamic_cast<OnsetPreSynapse*>(syn))
				{
					++nOnsetPresyn;
					continue;
				}
				if (dynamic_cast<ThresholdPreSynapse*>(syn))
				{
					++nThresholdPresyn;
					continue;
				}
				if (dynamic_cast<AlphaPostSynapse*>(syn))
				{
					++nAlphaPostsyn;
					continue;
				}
				if (dynamic_cast<Exp2PostSynapse*>(syn))
				{
					++nExp2Postsyn;
					continue;
				}
			}
		}

		// # onset presynapses
		BOOST_CHECK_MESSAGE(nOnsetPresyn == it->second.num_onset_presyn(),
		   "# OnsetPreSynapses don't match (Reference #" << it->second.num_onset_presyn()
		   << " vs. generated #" << nOnsetPresyn);

		// # threshold presynapses
		BOOST_CHECK_MESSAGE(nThresholdPresyn == it->second.num_threshold_presyn(),
		   "# ThresholdPreSynapses don't match (Reference #" << it->second.num_threshold_presyn()
		   << " vs. generated #" << nThresholdPresyn);

		// # alpha synapses
		BOOST_CHECK_MESSAGE(nAlphaPostsyn == it->second.num_alphasyn(),
		   "# AlphaPostSynapses don't match (Reference #" << it->second.num_alphasyn()
		   << " vs. generated #" << nAlphaPostsyn);

		// # exp2 synapses
		BOOST_CHECK_MESSAGE(nExp2Postsyn == it->second.num_exp2syn(),
		   "# Exp2PostSynapses don't match (Reference #" << it->second.num_exp2syn()
		   << " vs. generated #" << nExp2Postsyn);
	}
}

BOOST_AUTO_TEST_CASE(COMPARE_GENERATED_WITH_REFERNCE_GRID) {
	std::ifstream ifs1("small_network6.ugx");
	std::ifstream ifs2("small_network6_ref.ugx");

    std::istream_iterator<char> b1(ifs1), e1;
    std::istream_iterator<char> b2(ifs2), e2;

    BOOST_CHECK_EQUAL_COLLECTIONS(b1, e1, b2, e2);
}

BOOST_AUTO_TEST_SUITE_END();

