################################################################################
# cable
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

# compile and register nmodl-converted Transport only if # ENABLE_CONVERTED_Transport is set
option(CNConvertedTransport "ConvertedTransport")
message(STATUS "      ConvertedTransport:    " ${CNConvertedTransport} " (options are: ON, OFF)")

# compile and register nmodl-converted bap only if CNConvertedTransportbAP is set
option(CNConvertedTransportbAP "ConvertedTransportbAP")
message(STATUS "      ConvertedTransportbAP: " ${CNConvertedTransportbAP} " (options are: ON, OFF)")

# compile Converter
option(CNBuildConverter "BuildConverter" OFF)
message(STATUS "      BuildConverter:        " ${CNBuildConverter} " (options are: ON, OFF)")

# whether to build the test suite
option(CNTestsuite "Build Testsuite" OFF)
message(STATUS "      Testsuite:             " ${CNTestsuite} " (options are: ON, OFF)")



# cable equation sources
set (SOURCES
	cable_disc/cable_equation.cpp
	cable_disc/cable_equation_withOuterPot.cpp
	cable_disc/ElemDiscHH_base.cpp
	cable_disc/ElemDiscHH_fv1.cpp
	cable_disc/ElemDiscHH_Nernst_fv1.cpp
	cable_disc/ElemDiscHH_Nernst_neuron_fv1.cpp
	membrane_transport/cable_membrane_transport_interface.cpp
	membrane_transport/channel_hh.cpp
	membrane_transport/leakage.cpp
	membrane_transport/ion_leakage.cpp
	membrane_transport/na_k_pump.cpp
	membrane_transport/pmca.cpp
	membrane_transport/ncx.cpp
	membrane_transport/vdcc_bg.cpp
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

if (CNConvertedTransport)
	set (CONVERTED_TRANSPORT_ENABLED 1)
    include(membrane_transport/nmodl_converter/converted/sources.cmake)
endif (CNConvertedTransport)

if (CNConvertedTransportbAP)
	set (CONVERTED_TRANSPORT_BAP_ENABLED 1)
    include(membrane_transport/nmodl_converter/converted/sources_bAP.cmake)
endif (CNConvertedTransportbAP)

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

# synapse distributor sources
#set(SOURCES ${SOURCES} synapse_distributor/synapse_distributor.cpp)

# plugin source
set(SOURCES	${SOURCES} plugin_main.cpp)


# decide if Converter binary is to be built
if (CNBuildConverter)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
	set(SOURCES_CONVERTER
		membrane_transport/nmodl_converter/Converter.cpp
		membrane_transport/nmodl_converter/main.cpp)
	
	# add converter as a target
	add_executable(nmodl-converter ${SOURCES_CONVERTER})   
endif (CNBuildConverter)


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

