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
message(STATUS "Info: Additional options (-DCN) for the cable_neuron plugin")

# compile and register nmodl-converted Transport only if # ENABLE_CONVERTED_Transport is set
option(CNConvertedTransport "ConvertedTransport")
message(STATUS "      ConvertedTransport: " ${CNConvertedTransport} " (options are: ON, OFF)")

# compile Converter
option(CNBuildConverter "BuildConverter" OFF)
message(STATUS "      BuildConverter:     " ${CNBuildConverter} " (options are: ON, OFF)")

# whether to build the test suite
option(CNTestsuite "Build Testsuite" OFF)
message(STATUS "      Testsuite:          " ${CNTestsuite} " (options are: ON, OFF)")



# cable equation sources
set (SOURCES
	cable_disc/cable_equation.cpp
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
)

if (CNConvertedTransport)
    include(membrane_transport/nmodl_converter/converted/sources.cmake)
endif (CNConvertedTransport)

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



if(buildEmbeddedPlugins)
	# add the sources to ug4 sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
	
	if (CNConvertedTransport)
		# add definition to plugin_main.cpp to include all converted transport mechanisms
		ExportSingleFileDefinitions(${CMAKE_CURRENT_SOURCE_DIR}/plugin_main.cpp CONVERTED_TRANSPORT_ENABLED)
	endif (CNConvertedTransport)
	
	if (CNTestsuite)
		target_link_libraries (CNTestsuite ug4)
	endif (CNTestsuite)
	
else(buildEmbeddedPlugins)
	add_library(${pluginName} SHARED ${SOURCES})
	
	if (CNConvertedTransport)
		# add definition to plugin_main.cpp to include all converted transport mechanisms
		set_source_files_properties(plugin_main.cpp PROPERTIES COMPILE_DEFINITIONS CONVERTED_TRANSPORT_ENABLED)
	endif (CNConvertedTransport)
	
	target_link_libraries(${pluginName} ug4)
	
	if (CNTestsuite)
		target_link_libraries (CNTestsuite ${pluginName} ug4)
	endif (CNTestsuite)
endif(buildEmbeddedPlugins)

