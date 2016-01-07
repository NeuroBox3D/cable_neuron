################################################################################
# cable
#
# Plugin for cable equation and channel functionality
################################################################################

cmake_minimum_required(VERSION 2.6)

project(HH_Kabelnew)

# include the definitions and dependencies for ug-plugins.
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)


message(STATUS "")
message(STATUS "Info: Additional options (-DHH) for the HH_Kabelnew plugin")

# compile and register nmodl-converted channels only if
# ENABLE_CONVERTED_CHANNELS is set
option(HHConvertedChannels "ConvertedChannels")
message(STATUS "      ConvertedChannels:  " ${HHConvertedChannels} " (options are: ON, OFF)")


if (NOT HHConvChan)
	#message(STATUS "Converted channels compilation disabled. Enable by -DHHConvertedChannels=ON")
	set (ConvChanSources)
else (NOT HHConvChan)
	#message(STATUS "Converted channels compilation enabled. Disable by -DHHConvertedChannels=OFF")
    add_definitions(-DCONVERTED_CHANNELS_ENABLED)
    set (ConvChanSources
		Convert/Debug/Na2_converted_standard_UG.cpp 
		#Convert/Debug/Na_converted_standard_UG.cpp 
		Convert/Debug/na_converted_standard_UG.cpp 
		#Convert/Debug/PARAK75_converted_standard_UG.cpp 
		#Convert/Debug/cacum_converted_standard_UG.cpp 
		#Convert/Debug/gpi_converted_standard_UG.cpp 
		#Convert/Debug/NaL_converted_standard_UG.cpp 
		#Convert/Debug/HVA_converted_standard_UG.cpp 
		#Convert/Debug/HCN_converted_standard_UG.cpp 
		#Convert/Debug/KDRs_converted_standard_UG.cpp 
		#Convert/Debug/KDRf_converted_standard_UG.cpp 
		#Convert/Debug/SK_converted_standard_UG.cpp 
		#Convert/Debug/AXNODE75_converted_standard_UG.cpp 
		#Convert/Debug/nadifl_converted_standard_UG.cpp 
		#Convert/Debug/Kv1_converted_standard_UG.cpp 
		#Convert/Debug/hh_Ca_scaled_converted_standard_UG.cpp 
		#Convert/Debug/hh_Cs_scaled_converted_standard_UG.cpp 
		#Convert/Debug/hh_Cp_scaled_converted_standard_UG.cpp 
		#Convert/Debug/nav1p9_converted_standard_UG.cpp 
		#Convert/Debug/nav1p8_converted_standard_UG.cpp 
		#Convert/Debug/nattxs_converted_standard_UG.cpp 
		#Convert/Debug/nap_converted_standard_UG.cpp 
		#Convert/Debug/ks_converted_standard_UG.cpp 
		#Convert/Debug/kf_converted_standard_UG.cpp 
		#Convert/Debug/borgkm_converted_standard_UG.cpp 
		#Convert/Debug/iq_converted_standard_UG.cpp 
		#Convert/Debug/namir_converted_standard_UG.cpp 
		#Convert/Debug/kdr_converted_standard_UG.cpp 
		#Convert/Debug/naxn_converted_standard_UG.cpp 
		#Convert/Debug/na3n_converted_standard_UG.cpp 
		#Convert/Debug/kdrca1_converted_standard_UG.cpp 
		#Convert/Debug/kaprox_converted_standard_UG.cpp 
		#Convert/Debug/kadist_converted_standard_UG.cpp 
		#Convert/Debug/kslow_converted_standard_UG.cpp 
		#Convert/Debug/kfast_converted_standard_UG.cpp 
		#Convert/Debug/myseclamp_converted_standard_UG.cpp 
		#Convert/Debug/Kv4_csiosi_converted_standard_UG.cpp 
		#Convert/Debug/Kv4_csi_converted_standard_UG.cpp 
		#Convert/Debug/CaT_converted_allNernst_UG.cpp 
		#Convert/Debug/pump_converted_allNernst_UG.cpp 
		Convert/Debug/ca_converted_allNernst_UG.cpp 
		#Convert/Debug/pump_converted_standard_UG.cpp 
		#Convert/Debug/hh_converted_allNernst_UG.cpp 
		Convert/Debug/caL3d_converted_standard_UG.cpp 
		Convert/Debug/cad_converted_standard_UG.cpp 
		Convert/Debug/ar_converted_standard_UG.cpp 
		Convert/Debug/ca_converted_standard_UG.cpp 
		Convert/Debug/kca_converted_standard_UG.cpp 
		Convert/Debug/kir_converted_standard_UG.cpp 
		Convert/Debug/km_converted_standard_UG.cpp 
		Convert/Debug/inwardrect_converted_standard_UG.cpp 
		#Convert/Debug/CaT_converted_standard_UG.cpp 
		Convert/Debug/kv_converted_standard_UG.cpp  
		Convert/Debug/HH2_converted_standard_UG.cpp 
		Convert/Debug/h_converted_standard_UG.cpp 
		Convert/Debug/hh_converted_standard_UG.cpp 
		Convert/Debug/NMDA_Mg_converted_standard_UG.cpp 
		Convert/Debug/NMDA_Mg_T_converted_standard_UG.cpp 
		Convert/Debug/passive_converted_standard_UG.cpp 
		Convert/Debug/release_BMK_converted_standard_UG.cpp 
		Convert/Debug/release_exp_converted_standard_UG.cpp
		Convert/Debug/Kv4_csi_converted_standard_UG.cpp
		Convert/Debug/Kv4_csiosi_converted_standard_UG.cpp
		Convert/Debug/myseclamp_converted_standard_UG.cpp
	)
endif (NOT HHConvChan)

# sources
set (SOURCES
	ElemDiscHH_base.cpp
	Ca_PMCA.cpp
	Ca_NCX.cpp
	Na_K_Pump.cpp
	ElemDiscHH_fv1.cpp
	ElemDiscHH_Nernst_fv1.cpp
	ElemDiscHH_Nernst_neuron_fv1.cpp
	diam_attachment_handler.cpp
	channel_interface.cpp
	channel_hh.cpp
	leakage.cpp
	ion_leakage.cpp
	vdcc_bg.cpp
	VM_Disc.cpp
	order.cpp
	plugin_main.cpp
	${ConvChanSources}
)


# indicate to ug and all plugins that HHcable plugin is built
ExportDefinitions(-DPLUGIN_HH_KABEL_NEW_ENABLED)

if(buildEmbeddedPlugins)
	# add the sources to ug4 sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(HH_Kabelnew SHARED ${SOURCES})
	target_link_libraries(HH_Kabelnew ${pluginsToLink} ug4)
endif(buildEmbeddedPlugins)

