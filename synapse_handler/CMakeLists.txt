################################################################################
# CMakeLists for the the synapse_handler ug plugin
# Please note: to activate use -DSynapseHandler=ON during CMake configuration.
#
# created by: Stephan Grein
################################################################################

# minimum required CMake version for project below
cmake_minimum_required(VERSION 2.6)

# specify project name
project(synapse_handler)

# default values
set(SHC++0x OFF)
set(SHTestsuite ON)

# do not build testsuite in the embedded case
if(buildEmbeddedPlugins)
   set(SHTestsuite OFF)
endif(buildEmbeddedPlugins)

# set plugins to link
unset(pluginsToLink)

# options for building synapse_handler
message(STATUS "")
message(STATUS "Info: Additional options (-DSH) for the SynapseHandler (SH) plugin:")

option(SHTestsuite "Build Testsuite" ${SHTestsuite})
message(STATUS "Info: Testsuite:                " ${SHTestsuite} " (options are: ON, OFF)")
option(CGC++0x "Build C++0x " ${CGC++0x})
message(STATUS "Info: C++0x:                    " ${CGC++0x} " (options are: ON, OFF)")


# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

# set the sources and unit test sources
set(SOURCES 
    plugin_main.cpp
    util/utility.cpp
    grid/synapse_info.cpp
    grid/synapse_info_attachment_handler.cpp
    synapse_handler.cpp)
    
set(SOURCES_TEST 
    unit_tests/src/tests.cpp)
	
# decide if you want to build the boost testsuite executable (SHTestsuite)
if(${SHTestsuite} STREQUAL "ON")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
  add_executable(SHTestsuite ${SOURCES_TEST})
endif(${SHTestsuite} STREQUAL "ON")

# build project above with C++0x extensions (only .cpp files are affected)
IF(${CGC++0x} STREQUAL "ON")                                                          
  SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++0x")
ENDIF(${CGC++0x} STREQUAL "ON")


# create a shared library from the sources and link it against ug4.
if(buildEmbeddedPlugins)
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
	
	if(${SHTestsuite} STREQUAL "ON")
		target_link_libraries (SHTestsuite SynapseHandler ${pluginsToLink} ug4)
	endif(${SHTestsuite} STREQUAL "ON")
	
	if (SynapseDistributor)
		exportdefinitions(-DSH_SYNAPSE_DISTRIBUTOR_ENABLED)
	endif (SynapseDistributor)
	
else(buildEmbeddedPlugins)
    add_library(SynapseHandler SHARED ${SOURCES})
    
    if (SynapseDistributor)
		target_compile_definitions(SynapseHandler PUBLIC -DSH_SYNAPSE_DISTRIBUTOR_ENABLED)
		set(pluginsToLink ${pluginsToLink} SynapseDistributor)
	endif (SynapseDistributor)
	
    if (HH_Kabelnew)
		set(pluginsToLink ${pluginsToLink} HH_Kabelnew)
	endif (HH_Kabelnew)
	
    target_link_libraries(SynapseHandler ${pluginsToLink} ug4)
	if(${SHTestsuite} STREQUAL "ON")
		target_link_libraries (SHTestsuite SynapseHandler ${pluginsToLink} ug4)
	endif(${SHTestsuite} STREQUAL "ON")
endif(buildEmbeddedPlugins)
