################################################################################
# cable
#
# Plugin for cable equation and channel functionality
################################################################################

cmake_minimum_required(VERSION 2.6)

project(HH_Kabelnew)

set(HHSH ${HHSynapseHandler})

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")

# sources
set(SOURCES	ElemDiscHH_base.cpp
#Convert/Debug/ca_converted_allNernst_UG.cpp 
#Convert/Debug/caL3d_converted_standard_UG.cpp 
#Convert/Debug/pump_converted_standard_UG.cpp 
#Convert/Debug/cad_converted_standard_UG.cpp 
#Convert/Debug/ar_converted_standard_UG.cpp 
#Convert/Debug/ca_converted_allNernst_UG.cpp 
#Convert/Debug/kca_converted_standard_UG.cpp 
#Convert/Debug/kir_converted_standard_UG.cpp 
#Convert/Debug/km_converted_standard_UG.cpp 
#Convert/Debug/inwardrect_converted_standard_UG.cpp 
#Convert/Debug/CaT_converted_standard_UG.cpp 
#Convert/Debug/kv_converted_standard_UG.cpp 
#Convert/Debug/na_converted_standard_UG.cpp 
#Convert/Debug/HH2_converted_standard_UG.cpp 
#Convert/Debug/h_converted_standard_UG.cpp 
#Convert/Debug/hh_converted_standard_UG.cpp 
#Convert/Debug/NMDA_Mg_converted_standard_UG.cpp 
#Convert/Debug/NMDA_Mg_T_converted_standard_UG.cpp 
#Convert/Debug/passive_converted_standard_UG.cpp 
#Convert/Debug/release_BMK_converted_standard_UG.cpp 
#Convert/Debug/release_exp_converted_standard_UG.cpp
#Convert/Debug/Kv4_csi_converted_standard_UG.cpp
#Convert/Debug/Kv4_csiosi_converted_standard_UG.cpp
#Convert/Debug/myseclamp_converted_standard_UG.cpp
			ElemDiscHH_fv1.cpp
			ElemDiscHH_Nernst_fv1.cpp
			ElemDiscHH_Nernst_neuron_fv1.cpp
			diam_attachment_handler.cpp
			channel_interface.cpp
			channel_hh.cpp
			leakage.cpp
			ion_leakage.cpp
			vdcc_bg.cpp
			Ca_NCX.cpp
			Ca_PMCA.cpp
			VM_Disc.cpp
			order.cpp
			plugin_main.cpp
			)


# check for SynapseHandler availability and activate
message(STATUS "")
message(STATUS "Info: Additional options (-DHH) for the HH_Kabelnew (HH) plugin:")
option(HHSynapseHandler, "Build SynapseHandler support" ${HHSH})
message(STATUS "Info: SynapseHandler: " ${HHSH} " (options are: ON, OFF)")

if (NOT HHSynapseHandler OR NOT SynapseHandler)
    message(STATUS "SynapseHandler support deactivated. Enable by -DHHSynapseHandler=ON if -DHH_Kabelnew=ON")
else (NOT HHSynapseHandler OR NOT SynapseHandler)
    message(STATUS "SynapseHandler support activated. Disable by -DHHSynapseHandler=OFF")
    add_definitions(-DPLUGIN_SYNAPSE_HANDLER_ENABLED)
    ExportDefinitions(-DPLUGIN_SYNAPSE_HANDLER_ENABLED)
endif (NOT HHSynapseHandler OR NOT SynapseHandler)

if (NOT HHSynapseDistributor OR NOT SynapseDistributor)
    message(STATUS "SynapseDistributor support deactivated. Enable by -DHHSynapseDistributor=ON if -DHH_Kabelnew=ON")
else (NOT HHSynapseDistributor OR NOT SynapseDistributor)
    message(STATUS "SynapseDistributor support activated. Disable by -DHHSynapseDistributor=OFF")
    add_definitions(-DPLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED)
    ExportDefinitions(-DPLUGIN_SYNAPSE_DISTRIBUTOR_ENABLED)
	set(pluginsToLink ${pluginsToLink} SynapseDistributor)
endif (NOT HHSynapseDistributor OR NOT SynapseDistributor)

# indicate hh plugin will be build!
add_definitions(-DPLUGIN_HH_KABEL_NEW_ENABLED)
ExportDefinitions(-DPLUGIN_HH_KABEL_NEW_ENABLED)

if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/HH_Kabelnew ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(HH_Kabelnew SHARED ${SOURCES})
	target_link_libraries (HH_Kabelnew ${pluginsToLink} ug4)
endif(buildEmbeddedPlugins)



