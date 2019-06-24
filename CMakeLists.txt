##
 # Copyright (c) 2009-2019: G-CSC, Goethe University Frankfurt
 #
 # Authors: Markus Breit, Pascal Gottmann
 # Creation date: 2014-09-02
 #
 # This file is part of NeuroBox, which is based on UG4.
 #
 # NeuroBox and UG4 are free software: You can redistribute it and/or modify it
 # under the terms of the GNU Lesser General Public License version 3
 # (as published by the Free Software Foundation) with the following additional
 # attribution requirements (according to LGPL/GPL v3 §7):
 #
 # (1) The following notice must be displayed in the appropriate legal notices
 # of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 #
 # (2) The following notice must be displayed at a prominent place in the
 # terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 #
 # (3) The following bibliography is recommended for citation and must be
 # preserved in all covered files:
 # "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 #   parallel geometric multigrid solver on hierarchically distributed grids.
 #   Computing and visualization in science 16, 4 (2013), 151-164"
 # "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 #   flexible software system for simulating PDE based models on high performance
 #   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 # "Stepniewski, M., Breit, M., Hoffer, M. and Queisser, G.
 #   NeuroBox: computational mathematics in multiscale neuroscience.
 #   Computing and visualization in science (2019).
 # "Breit, M. et al. Anatomically detailed and large-scale simulations studying
 #   synapse loss and synchrony using NeuroBox. Front. Neuroanat. 10 (2016), 8"
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 # GNU Lesser General Public License for more details.
 ##
 
################################################################################
# cable_neuron
#
# Plugin for cable equation and channel functionality
################################################################################

cmake_minimum_required(VERSION 2.6)

set(pluginName cable_neuron)

project(cable_neuron)

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


message(STATUS "")
message(STATUS "Additional options (-DCN) for the cable_neuron plugin")

# whether to build the test suite
option(CNTestsuite "Build Testsuite" OFF)
message(STATUS "      Testsuite:   " ${CNTestsuite} " (options are: ON, OFF)")



# cable equation sources
set (SOURCES
	cable_disc/cable_equation.cpp
	cable_disc/cable_equation_withOuterPot.cpp
	cable_disc/implicit_active_cable_disc_base.cpp
	cable_disc/implicit_active_cable_disc.cpp
	cable_disc/implicit_active_cable_disc_nernst.cpp
	membrane_transport/cable_membrane_transport_interface.cpp
	membrane_transport/channel_hh.cpp
	membrane_transport/leakage.cpp
	membrane_transport/ion_leakage.cpp
	membrane_transport/na_k_pump.cpp
	membrane_transport/pmca.cpp
	membrane_transport/ncx.cpp
	membrane_transport/vdcc_bg.cpp
	membrane_transport/ka_golding01.cpp
	membrane_transport/kdr_golding01.cpp
	membrane_transport/nax_golding01.cpp
	util/diam_attachment_handler.cpp
	util/order.cpp
	util/functions.cpp
	util/neuronal_topology_importer.cpp
)

if (PARALLEL)
	set (SOURCES ${SOURCES} util/simple_cable_neuron_partitioner.cpp)
endif (PARALLEL)

if (Parmetis)
	set (NC_WITH_PARMETIS 1)
	set (SOURCES ${SOURCES} util/cable_neuron_unificator.cpp)
endif (Parmetis)

# check whether compiler supports C++11 standard (required for HOC import)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if (COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c++11")
	set(NETI_ALL_FEATURES 1)
else (COMPILER_SUPPORTS_CXX11)
	message(STATUS "      Compiler does not support C++11: HOC import will not be available.")
endif (COMPILER_SUPPORTS_CXX11)

# synapse_handling sources
include(synapse_handling/sources.cmake)

# plugin source
set(SOURCES	${SOURCES} plugin_main.cpp)


# decide if boost testsuite executable (CNTestsuite) is to be built
if (CNTestsuite)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
	set(SOURCES_TEST unit_tests/src/tests.cpp)
	
	# add test suite as a target
	add_executable(CNTestsuite ${SOURCES_TEST})
endif (CNTestsuite)



# This will copy the configurable file config_cmake.h to
# the final configuration header cn_config.h,
# setting the #defines as given by the current CMake variables.
configure_file(config_cmake.h ${CMAKE_BINARY_DIR}/cn_config.h)


if(buildEmbeddedPlugins)
	# add the sources to ug4 sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
	
	if (CNTestsuite)
		target_link_libraries (CNTestsuite ug4)
	endif (CNTestsuite)
	
else(buildEmbeddedPlugins)
	add_library(${pluginName} SHARED ${SOURCES})
		
	if (Parmetis)
		# make plugin link against Parmetis
		set(linkLibraries Parmetis ${linkLibraries})
	endif (Parmetis)
	
	target_link_libraries(${pluginName} ug4)
	
	if (CNTestsuite)
		target_link_libraries (CNTestsuite ${pluginName} ug4)
	endif (CNTestsuite)
endif(buildEmbeddedPlugins)

