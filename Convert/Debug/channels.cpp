
 
 
#include "hh_converted_standard_UG.h" 
 
 
{ 
	 typedef hh_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("hh_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "hh_converted_standard_UG", tag); 
} 
 
 
 
#include "passive_converted_standard_UG.h" 
 
 
{ 
	 typedef passive_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("passive_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "passive_converted_standard_UG", tag); 
} 
 
 
