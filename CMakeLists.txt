################################################################################
# HH_Kabel
#
# Plugin for Hodkin und Huxley + Kabelgleichung
################################################################################

cmake_minimum_required(VERSION 2.6)

project(HH_Kabelnew)

set(HHSynapseProvider ON)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

# sources
set(SOURCES	ElemDiscHH_base.cpp   
##Convert/Debug/release_exp_converted_standard_UG.cpp 
Convert/Debug/cad_converted_standard_UG.cpp 
Convert/Debug/ar_converted_standard_UG.cpp 
Convert/Debug/ca_converted_standard_UG.cpp 
Convert/Debug/kca_converted_standard_UG.cpp 
Convert/Debug/kir_converted_standard_UG.cpp 
Convert/Debug/km_converted_standard_UG.cpp 
Convert/Debug/inwardrect_converted_standard_UG.cpp 
Convert/Debug/CaT_converted_standard_UG.cpp 
Convert/Debug/kv_converted_standard_UG.cpp 
Convert/Debug/na_converted_standard_UG.cpp 
Convert/Debug/HH2_converted_standard_UG.cpp 
Convert/Debug/h_converted_standard_UG.cpp 
Convert/Debug/hh_converted_standard_UG.cpp 
Convert/Debug/NMDA_Mg_converted_standard_UG.cpp 
Convert/Debug/NMDA_Mg_T_converted_standard_UG.cpp 
Convert/Debug/passive_converted_standard_UG.cpp 
Convert/Debug/release_BMK_converted_standard_UG.cpp 
			ElemDiscHH_fv1.cpp
			ElemDiscHH_Nernst_fv1.cpp
			ElemDiscHH_Nernst_neuron_fv1.cpp
			channel_interface.cpp
			VM_Disc.cpp
			plugin_main.cpp
			)


# check for SynapseProvider availability and activate
#
message(STATUS "")
message(STATUS "Info: Additional options for the HH_Kabelnew plugin:")
if (NOT HHSynapseProvider OR NOT SynapseProvider)
    message(STATUS "SynapseProvider support deactivated. Enable by -DHHSynapseProvider=ON if -DHH_Kabelnew=ON")
else (NOT HHSynapseProvider OR NOT SynapseProvider)
    message(STATUS "SynapseProvider support activated.")
    add_definitions(-DPLUGIN_SYNAPSE_PROVIDER_ENABLED)
    ExportDefinitions(-DPLUGIN_SYNAPSE_PROVIDER_ENABLED)
    set(pluginsToLink ${pluginsToLink} SynapseProvider)
endif (NOT HHSynapseProvider OR NOT SynapseProvider)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/HH_Kabelnew ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(HH_Kabelnew SHARED ${SOURCES})
	target_link_libraries (HH_Kabelnew ${pluginsToLink} ug4)
endif(buildEmbeddedPlugins)


















































































































































































































































































