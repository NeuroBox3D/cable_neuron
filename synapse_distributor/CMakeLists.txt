################################################################################
#
################################################################################

set(pluginName	SynapseDistributor)
set(SOURCES		plugin_main.cpp
				synapse_distributor.cpp)
			
if (NOT SynapseHandler)
    message(FATAL_ERROR "ERROR: SynapseDistributor needs SynapseHandler plugin. Please activate by setting -DSynapseHandler=ON.")
endif (NOT SynapseHandler)			


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.6)

project(UG_PLUGIN_${pluginName})

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	# create a shared library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} SynapseHandler ug4)
endif(buildEmbeddedPlugins)
