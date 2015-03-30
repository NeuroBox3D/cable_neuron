#include "cad_converted_standard_UG.h" 
 
 
{ 
	 typedef cad_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("cad_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "cad_converted_standard_UG", tag); 
} 
 
 
#include "ar_converted_standard_UG.h" 
 
 
{ 
	 typedef ar_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("ar_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "ar_converted_standard_UG", tag); 
} 
 
 
#include "ca_converted_standard_UG.h" 
 
 
{ 
	 typedef ca_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("ca_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "ca_converted_standard_UG", tag); 
} 
 
 
#include "kca_converted_standard_UG.h" 
 
 
{ 
	 typedef kca_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("kca_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kca_converted_standard_UG", tag); 
} 
 
 
#include "kir_converted_standard_UG.h" 
 
 
{ 
	 typedef kir_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("kir_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kir_converted_standard_UG", tag); 
} 
 
 
#include "km_converted_standard_UG.h" 
 
 
{ 
	 typedef km_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("km_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "km_converted_standard_UG", tag); 
} 
 
 
#include "inwardrect_converted_standard_UG.h" 
 
 
{ 
	 typedef inwardrect_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("inwardrect_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "inwardrect_converted_standard_UG", tag); 
} 
 
 
#include "CaT_converted_standard_UG.h" 
 
 
{ 
	 typedef CaT_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("CaT_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "CaT_converted_standard_UG", tag); 
} 
 
 
#include "kv_converted_standard_UG.h" 
 
 
{ 
	 typedef kv_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("kv_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kv_converted_standard_UG", tag); 
} 
 
 
#include "na_converted_standard_UG.h" 
 
 
{ 
	 typedef na_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("na_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "na_converted_standard_UG", tag); 
} 
 
 
#include "HH2_converted_standard_UG.h" 
 
 
{ 
	 typedef HH2_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("HH2_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "HH2_converted_standard_UG", tag); 
} 
 
 
#include "h_converted_standard_UG.h" 
 
 
{ 
	 typedef h_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("h_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "h_converted_standard_UG", tag); 
} 
 
 
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
 
 
#include "NMDA_Mg_converted_standard_UG.h" 
 
 
{ 
	 typedef NMDA_Mg_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("NMDA_Mg_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "NMDA_Mg_converted_standard_UG", tag); 
} 
 
 
#include "NMDA_Mg_T_converted_standard_UG.h" 
 
 
{ 
	 typedef NMDA_Mg_T_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("NMDA_Mg_T_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "NMDA_Mg_T_converted_standard_UG", tag); 
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
 
 
#include "release_BMK_converted_standard_UG.h" 
 
 
{ 
	 typedef release_BMK_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("release_BMK_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "release_BMK_converted_standard_UG", tag); 
} 
 
 
#include "release_exp_converted_standard_UG.h" 
 
 
{ 
	 typedef release_exp_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("release_exp_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "release_exp_converted_standard_UG", tag); 
} 
 
 
