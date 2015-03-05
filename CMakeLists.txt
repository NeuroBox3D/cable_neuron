################################################################################
# HH_Kabel
#
# Plugin for Hodkin und Huxley + Kabelgleichung
################################################################################

cmake_minimum_required(VERSION 2.6)

project(HH_Kabelnew)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(SOURCES	ElemDiscHH_base.cpp
			ElemDiscHH_fv1.cpp
			ElemDiscHH_Nernst_fv1.cpp
			ElemDiscHH_Nernst_neuron_fv1.cpp
			channel_interface.cpp
			#hh_converted_UG.cpp
			VM_Disc.cpp
			plugin_main.cpp

			)




if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/HH_Kabelnew ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(HH_Kabelnew SHARED ${SOURCES})
	target_link_libraries (HH_Kabelnew ug4)
endif(buildEmbeddedPlugins)
