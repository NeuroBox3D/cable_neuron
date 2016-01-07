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

# compile and register nmodl-converted channels only if # ENABLE_CONVERTED_CHANNELS is set
option(CNConvertedChannels "ConvertedChannels")
message(STATUS "      ConvertedChannels:  " ${CNConvertedChannels} " (options are: ON, OFF)")

# whether to build the test suite
option(CNTestsuite "Build Testsuite" OFF)
message(STATUS "      Testsuite:          " ${CNTestsuite} " (options are: ON, OFF)")



# cable equation sources
set (SOURCES
	cable_equation/cable_equation.cpp
	cable_equation/ElemDiscHH_base.cpp
	cable_equation/ElemDiscHH_fv1.cpp
	cable_equation/ElemDiscHH_Nernst_fv1.cpp
	cable_equation/ElemDiscHH_Nernst_neuron_fv1.cpp
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
)

if (CNConvertedChannels)
    set (SOURCES ${SOURCES}
		membrane_transport/nmodl_converter/converted/Na2_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/Na_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/na_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/PARAK75_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/cacum_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/gpi_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/NaL_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/HVA_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/HCN_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/KDRs_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/KDRf_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/SK_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/AXNODE75_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/nadifl_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/Kv1_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/hh_Ca_scaled_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/hh_Cs_scaled_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/hh_Cp_scaled_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/nav1p9_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/nav1p8_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/nattxs_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/nap_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/ks_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/kf_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/borgkm_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/iq_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/namir_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/kdr_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/naxn_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/na3n_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/kdrca1_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/kaprox_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/kadist_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/kslow_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/kfast_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/myseclamp_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/Kv4_csiosi_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/Kv4_csi_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/CaT_converted_allNernst_UG.cpp 
		#membrane_transport/nmodl_converter/converted/pump_converted_allNernst_UG.cpp 
		membrane_transport/nmodl_converter/converted/ca_converted_allNernst_UG.cpp 
		#membrane_transport/nmodl_converter/converted/pump_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/hh_converted_allNernst_UG.cpp 
		membrane_transport/nmodl_converter/converted/caL3d_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/cad_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/ar_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/ca_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/kca_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/kir_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/km_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/inwardrect_converted_standard_UG.cpp 
		#membrane_transport/nmodl_converter/converted/CaT_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/kv_converted_standard_UG.cpp  
		membrane_transport/nmodl_converter/converted/HH2_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/h_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/hh_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/NMDA_Mg_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/NMDA_Mg_T_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/passive_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/release_BMK_converted_standard_UG.cpp 
		membrane_transport/nmodl_converter/converted/release_exp_converted_standard_UG.cpp
		membrane_transport/nmodl_converter/converted/Kv4_csi_converted_standard_UG.cpp
		membrane_transport/nmodl_converter/converted/Kv4_csiosi_converted_standard_UG.cpp
		membrane_transport/nmodl_converter/converted/myseclamp_converted_standard_UG.cpp
	)
endif (CNConvertedChannels)


# synapse_handler sources
set(SOURCES ${SOURCES}
    synapse_handler/util/utility.cpp
    synapse_handler/grid/synapse_info.cpp
    synapse_handler/grid/synapse_info_attachment_handler.cpp
    synapse_handler/synapse_handler.cpp
)

# synapse distributor sources
set(SOURCES	${SOURCES} synapse_distributor/synapse_distributor.cpp)

# plugin source
set(SOURCES	${SOURCES} plugin_main.cpp)


# decide if boost testsuite executable (CNTestsuite) is to be built
if (CNTestsuite)
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${UG_ROOT_PATH}/bin/)
	set(SOURCES_TEST unit_tests/src/tests.cpp)
	add_executable(CNTestsuite ${SOURCES_TEST})   
endif (CNTestsuite)



if(buildEmbeddedPlugins)
	# add the sources to ug4 sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
	
	if (CNConvertedChannels)
		exportdefinitions(-DCONVERTED_CHANNELS_ENABLED)
	endif (CNConvertedChannels)
	
	if (CNTestsuite)
		target_link_libraries (SHTestsuite ug4)
	endif (CNTestsuite)
	
else(buildEmbeddedPlugins)
	add_library(${pluginName} SHARED ${SOURCES})
	
	if (CNConvertedChannels)
		target_compile_definitions(${pluginName} PUBLIC -DCONVERTED_CHANNELS_ENABLED)
	endif (CNConvertedChannels)
	
	target_link_libraries(${pluginName} ug4)
	
	if (CNTestsuite)
		target_link_libraries (CNTestsuite ${pluginName} ug4)
	endif (CNTestsuite)
endif(buildEmbeddedPlugins)

