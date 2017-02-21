set (SOURCES ${SOURCES}
		synapse_handling/synapses/base_synapse.cpp
		synapse_handling/synapses/pre_synapse.cpp
		synapse_handling/synapses/post_synapse.cpp
		
		synapse_handling/synapses/onset_pre_synapse.cpp
		synapse_handling/synapses/alpha_post_synapse.cpp
		
		synapse_handling/synapses/threshold_pre_synapse.cpp
		synapse_handling/synapses/exp2_post_synapse.cpp
		
		synapse_handling/synapse_container.cpp
		
		synapse_handling/synapse_handler.cpp
		synapse_handling/synapse_attachment_handler.cpp
		synapse_handling/synapse_attachment_serializer.cpp
		synapse_handling/synapse_dealer.cpp
		
		synapse_handling/synapse_distributor.cpp
	)
