 
 
 
#include "../converted/ar_converted_standard_UG.h"
 
 
{ 
	 typedef ar_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("ar_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setg0" , &T::setg0)
.add_method("sete" , &T::sete)
.add_method("setc" , &T::setc)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "ar_converted_standard_UG", tag); 
} 
 
 
#include "../converted/ca_converted_standard_UG.h" 
 
 
{ 
	 typedef ca_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("ca_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setvshift" , &T::setvshift)
.add_method("setcao" , &T::setcao)
.add_method("setcai" , &T::setcai)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "ca_converted_standard_UG", tag); 
} 
 
 
#include "../converted/cad_converted_standard_UG.h" 
 
 
{ 
	 typedef cad_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("cad_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setdepth" , &T::setdepth)
.add_method("settaur" , &T::settaur)
.add_method("setcainf" , &T::setcainf)
.add_method("setcai" , &T::setcai)
.add_method("set_log_caSGate" , &T::set_log_caSGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "cad_converted_standard_UG", tag); 
} 
 
 
 
#include "../converted/CaT_converted_standard_UG.h" 
 
 
{ 
	 typedef CaT_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("CaT_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setvshift" , &T::setvshift)
.add_method("setcao" , &T::setcao)
.add_method("setcai" , &T::setcai)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("setv12m" , &T::setv12m)
.add_method("setv12h" , &T::setv12h)
.add_method("setvwm" , &T::setvwm)
.add_method("setvwh" , &T::setvwh)
.add_method("setam" , &T::setam)
.add_method("setah" , &T::setah)
.add_method("setvm1" , &T::setvm1)
.add_method("setvm2" , &T::setvm2)
.add_method("setvh1" , &T::setvh1)
.add_method("setvh2" , &T::setvh2)
.add_method("setwm1" , &T::setwm1)
.add_method("setwm2" , &T::setwm2)
.add_method("setwh1" , &T::setwh1)
.add_method("setwh2" , &T::setwh2)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "CaT_converted_standard_UG", tag); 
} 
 
 
#include "../converted/HH2_converted_standard_UG.h" 
 
 
{ 
	 typedef HH2_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("HH2_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setena" , &T::setena)
.add_method("setek" , &T::setek)
.add_method("setvtraub" , &T::setvtraub)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "HH2_converted_standard_UG", tag); 
} 
 
 
#include "../converted/inwardrect_converted_standard_UG.h" 
 
 
{ 
	 typedef inwardrect_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("inwardrect_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("settha" , &T::settha)
.add_method("setqa" , &T::setqa)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "inwardrect_converted_standard_UG", tag); 
} 
 
 
#include "../converted/kca_converted_standard_UG.h" 
 
 
{ 
	 typedef kca_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kca_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setcai" , &T::setcai)
.add_method("setcaix" , &T::setcaix)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kca_converted_standard_UG", tag); 
} 
 
 
#include "../converted/kir_converted_standard_UG.h" 
 
 
{ 
	 typedef kir_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kir_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgkbar" , &T::setgkbar)
.add_method("setmvhalf" , &T::setmvhalf)
.add_method("setmslope" , &T::setmslope)
.add_method("setmshift" , &T::setmshift)
.add_method("setqfact" , &T::setqfact)
.add_method("set_log_mGate" , &T::set_log_mGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kir_converted_standard_UG", tag); 
} 
 
 
#include "../converted/km_converted_standard_UG.h" 
 
 
{ 
	 typedef km_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("km_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("settha" , &T::settha)
.add_method("setqa" , &T::setqa)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "km_converted_standard_UG", tag); 
} 
 
 
#include "../converted/kv_converted_standard_UG.h" 
 
 
{ 
	 typedef kv_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kv_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("settha" , &T::settha)
.add_method("setqa" , &T::setqa)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kv_converted_standard_UG", tag); 
} 
 

#include "../converted/na_converted_standard_UG.h" 
 
 
{ 
	 typedef na_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("na_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setvshift" , &T::setvshift)
.add_method("settha" , &T::settha)
.add_method("setqa" , &T::setqa)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("setthi1" , &T::setthi1)
.add_method("setthi2" , &T::setthi2)
.add_method("setqi" , &T::setqi)
.add_method("setthinf" , &T::setthinf)
.add_method("setqinf" , &T::setqinf)
.add_method("setRg" , &T::setRg)
.add_method("setRd" , &T::setRd)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "na_converted_standard_UG", tag); 
} 

 
#include "../converted/NMDA_Mg_T_converted_standard_UG.h" 
 
 
{ 
	 typedef NMDA_Mg_T_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("NMDA_Mg_T_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setErev" , &T::setErev)
.add_method("setgmax" , &T::setgmax)
.add_method("setmg" , &T::setmg)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("setvalence" , &T::setvalence)
.add_method("setmemb_fraction" , &T::setmemb_fraction)
.add_method("setRb" , &T::setRb)
.add_method("setRu" , &T::setRu)
.add_method("setRo" , &T::setRo)
.add_method("setRc" , &T::setRc)
.add_method("setRd1" , &T::setRd1)
.add_method("setRr1" , &T::setRr1)
.add_method("setRd2" , &T::setRd2)
.add_method("setRr2" , &T::setRr2)
.add_method("setRmb" , &T::setRmb)
.add_method("setRmu" , &T::setRmu)
.add_method("setRmc1b" , &T::setRmc1b)
.add_method("setRmc1u" , &T::setRmc1u)
.add_method("setRmc2b" , &T::setRmc2b)
.add_method("setRmc2u" , &T::setRmc2u)
.add_method("setRmd1b" , &T::setRmd1b)
.add_method("setRmd1u" , &T::setRmd1u)
.add_method("setRmd2b" , &T::setRmd2b)
.add_method("setRmd2u" , &T::setRmd2u)
.add_method("setRbMg" , &T::setRbMg)
.add_method("setRuMg" , &T::setRuMg)
.add_method("setRoMg" , &T::setRoMg)
.add_method("setRcMg" , &T::setRcMg)
.add_method("setRd1Mg" , &T::setRd1Mg)
.add_method("setRr1Mg" , &T::setRr1Mg)
.add_method("setRd2Mg" , &T::setRd2Mg)
.add_method("setRr2Mg" , &T::setRr2Mg)
.add_method("setC" , &T::setC)
.add_method("set_log_UGate" , &T::set_log_UGate)
.add_method("set_log_FGate" , &T::set_log_FGate)
.add_method("set_log_ClGate" , &T::set_log_ClGate)
.add_method("set_log_D1Gate" , &T::set_log_D1Gate)
.add_method("set_log_D2Gate" , &T::set_log_D2Gate)
.add_method("set_log_OGate" , &T::set_log_OGate)
.add_method("set_log_UMgGate" , &T::set_log_UMgGate)
.add_method("set_log_ClMgGate" , &T::set_log_ClMgGate)
.add_method("set_log_D1MgGate" , &T::set_log_D1MgGate)
.add_method("set_log_D2MgGate" , &T::set_log_D2MgGate)
.add_method("set_log_OMgGate" , &T::set_log_OMgGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "NMDA_Mg_T_converted_standard_UG", tag); 
} 
 
 

 
#include "../converted/passive_converted_standard_UG.h" 
 
 
{ 
	 typedef passive_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("passive_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setg" , &T::setg)
.add_method("sete" , &T::sete)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "passive_converted_standard_UG", tag); 
} 
 
 
#include "../converted/release_BMK_converted_standard_UG.h" 
 
 
{ 
	 typedef release_BMK_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("release_BMK_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setdel" , &T::setdel)
.add_method("setdur" , &T::setdur)
.add_method("setamp" , &T::setamp)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "release_BMK_converted_standard_UG", tag); 
} 
 
 
#include "../converted/release_exp_converted_standard_UG.h" 
 
 
{ 
	 typedef release_exp_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("release_exp_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("settau1" , &T::settau1)
.add_method("settau2" , &T::settau2)
.add_method("set_log_AGate" , &T::set_log_AGate)
.add_method("set_log_BGate" , &T::set_log_BGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "release_exp_converted_standard_UG", tag); 
} 
 
 
#include "../converted/caL3d_converted_standard_UG.h" 
 
 
{ 
	 typedef caL3d_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("caL3d_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setp" , &T::setp)
.add_method("setth" , &T::setth)
.add_method("setq" , &T::setq)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("set_log_CGate" , &T::set_log_CGate)
.add_method("set_log_OGate" , &T::set_log_OGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "caL3d_converted_standard_UG", tag); 
} 
 

 
 
#include "../converted/pump_converted_standard_UG.h" 
 
 
{ 
	 typedef pump_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("pump_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setnai" , &T::setnai)
.add_method("setipumpmax" , &T::setipumpmax)
.add_method("setkm" , &T::setkm)
.add_method("setn" , &T::setn)
.add_method("setnainit" , &T::setnainit)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "pump_converted_standard_UG", tag); 
} 
 
 
#include "../converted/h_converted_standard_UG.h" 
 
 
{ 
	 typedef h_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("h_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setehd" , &T::setehd)
.add_method("setghdbar" , &T::setghdbar)
.add_method("setvhalfl" , &T::setvhalfl)
.add_method("setkl" , &T::setkl)
.add_method("setvhalft" , &T::setvhalft)
.add_method("seta0t" , &T::seta0t)
.add_method("setzetat" , &T::setzetat)
.add_method("setgmt" , &T::setgmt)
.add_method("setq10" , &T::setq10)
.add_method("setqtl" , &T::setqtl)
.add_method("set_log_lGate" , &T::set_log_lGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "h_converted_standard_UG", tag); 
} 
 
 
#include "../converted/hh_converted_standard_UG.h" 
 
 
{ 
	 typedef hh_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("hh_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setgl" , &T::setgl)
.add_method("setel" , &T::setel)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "hh_converted_standard_UG", tag); 
} 
 
 
#include "../converted/ca_converted_allNernst_UG.h" 
 
 
{ 
	 typedef ca_converted_allNernst_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("ca_converted_allNernst_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setvshift" , &T::setvshift)
.add_method("setcao" , &T::setcao)
.add_method("setcai" , &T::setcai)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "ca_converted_allNernst_UG", tag); 
} 
 
 
#include "pump_converted_allNernst_UG.h" 
 
 
{ 
	 typedef pump_converted_allNernst_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("pump_converted_allNernst_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setnai" , &T::setnai)
.add_method("setipumpmax" , &T::setipumpmax)
.add_method("setkm" , &T::setkm)
.add_method("setn" , &T::setn)
.add_method("setnainit" , &T::setnainit)
.add_method("setcelsius" , &T::setcelsius)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "pump_converted_allNernst_UG", tag); 
} 
 
 
#include "hh_converted_allNernst_UG.h" 
 
 
{ 
	 typedef hh_converted_allNernst_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("hh_converted_allNernst_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setgl" , &T::setgl)
.add_method("setel" , &T::setel)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "hh_converted_allNernst_UG", tag); 
} 
 
 
#include "CaT_converted_allNernst_UG.h" 
 
 
{ 
	 typedef CaT_converted_allNernst_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("CaT_converted_allNernst_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setvshift" , &T::setvshift)
.add_method("setcao" , &T::setcao)
.add_method("setcai" , &T::setcai)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("setv12m" , &T::setv12m)
.add_method("setv12h" , &T::setv12h)
.add_method("setvwm" , &T::setvwm)
.add_method("setvwh" , &T::setvwh)
.add_method("setam" , &T::setam)
.add_method("setah" , &T::setah)
.add_method("setvm1" , &T::setvm1)
.add_method("setvm2" , &T::setvm2)
.add_method("setvh1" , &T::setvh1)
.add_method("setvh2" , &T::setvh2)
.add_method("setwm1" , &T::setwm1)
.add_method("setwm2" , &T::setwm2)
.add_method("setwh1" , &T::setwh1)
.add_method("setwh2" , &T::setwh2)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "CaT_converted_allNernst_UG", tag); 
} 
 
 
#include "../converted/NMDA_Mg_converted_standard_UG.h" 
 
 
{ 
	 typedef NMDA_Mg_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("NMDA_Mg_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setErev" , &T::setErev)
.add_method("setgmax" , &T::setgmax)
.add_method("setmg" , &T::setmg)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
.add_method("setvalence" , &T::setvalence)
.add_method("setmemb_fraction" , &T::setmemb_fraction)
.add_method("setRb" , &T::setRb)
.add_method("setRu" , &T::setRu)
.add_method("setRo" , &T::setRo)
.add_method("setRc" , &T::setRc)
.add_method("setRd1" , &T::setRd1)
.add_method("setRr1" , &T::setRr1)
.add_method("setRd2" , &T::setRd2)
.add_method("setRr2" , &T::setRr2)
.add_method("setRmb" , &T::setRmb)
.add_method("setRmu" , &T::setRmu)
.add_method("setRmc1b" , &T::setRmc1b)
.add_method("setRmc1u" , &T::setRmc1u)
.add_method("setRmc2b" , &T::setRmc2b)
.add_method("setRmc2u" , &T::setRmc2u)
.add_method("setRmd1b" , &T::setRmd1b)
.add_method("setRmd1u" , &T::setRmd1u)
.add_method("setRmd2b" , &T::setRmd2b)
.add_method("setRmd2u" , &T::setRmd2u)
.add_method("setRbMg" , &T::setRbMg)
.add_method("setRuMg" , &T::setRuMg)
.add_method("setRoMg" , &T::setRoMg)
.add_method("setRcMg" , &T::setRcMg)
.add_method("setRd1Mg" , &T::setRd1Mg)
.add_method("setRr1Mg" , &T::setRr1Mg)
.add_method("setRd2Mg" , &T::setRd2Mg)
.add_method("setRr2Mg" , &T::setRr2Mg)
.add_method("setC" , &T::setC)
.add_method("set_log_UGate" , &T::set_log_UGate)
.add_method("set_log_ClGate" , &T::set_log_ClGate)
.add_method("set_log_D1Gate" , &T::set_log_D1Gate)
.add_method("set_log_D2Gate" , &T::set_log_D2Gate)
.add_method("set_log_OGate" , &T::set_log_OGate)
.add_method("set_log_UMgGate" , &T::set_log_UMgGate)
.add_method("set_log_ClMgGate" , &T::set_log_ClMgGate)
.add_method("set_log_D1MgGate" , &T::set_log_D1MgGate)
.add_method("set_log_D2MgGate" , &T::set_log_D2MgGate)
.add_method("set_log_OMgGate" , &T::set_log_OMgGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "NMDA_Mg_converted_standard_UG", tag); 
} 
 
 
#include "../converted/Kv4_csi_converted_standard_UG.h" 
 
 
{ 
	 typedef Kv4_csi_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("Kv4_csi_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgmax" , &T::setgmax)
.add_method("setek" , &T::setek)
.add_method("seta" , &T::seta)
.add_method("setza" , &T::setza)
.add_method("setb" , &T::setb)
.add_method("setzb" , &T::setzb)
.add_method("setc" , &T::setc)
.add_method("setzc" , &T::setzc)
.add_method("setd" , &T::setd)
.add_method("setzd" , &T::setzd)
.add_method("setk" , &T::setk)
.add_method("setzk" , &T::setzk)
.add_method("setl" , &T::setl)
.add_method("setzl" , &T::setzl)
.add_method("setf" , &T::setf)
.add_method("setq" , &T::setq)
.add_method("setkci" , &T::setkci)
.add_method("setkic" , &T::setkic)
.add_method("set_log_C0Gate" , &T::set_log_C0Gate)
.add_method("set_log_C1Gate" , &T::set_log_C1Gate)
.add_method("set_log_C2Gate" , &T::set_log_C2Gate)
.add_method("set_log_C3Gate" , &T::set_log_C3Gate)
.add_method("set_log_C4Gate" , &T::set_log_C4Gate)
.add_method("set_log_C5Gate" , &T::set_log_C5Gate)
.add_method("set_log_I0Gate" , &T::set_log_I0Gate)
.add_method("set_log_I1Gate" , &T::set_log_I1Gate)
.add_method("set_log_I2Gate" , &T::set_log_I2Gate)
.add_method("set_log_I3Gate" , &T::set_log_I3Gate)
.add_method("set_log_I4Gate" , &T::set_log_I4Gate)
.add_method("set_log_I5Gate" , &T::set_log_I5Gate)
.add_method("set_log_OGate" , &T::set_log_OGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "Kv4_csi_converted_standard_UG", tag); 
} 
 
 
#include "../converted/Kv4_csiosi_converted_standard_UG.h" 
 
 
{ 
	 typedef Kv4_csiosi_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("Kv4_csiosi_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgmax" , &T::setgmax)
.add_method("setek" , &T::setek)
.add_method("seta" , &T::seta)
.add_method("setza" , &T::setza)
.add_method("setb" , &T::setb)
.add_method("setzb" , &T::setzb)
.add_method("setc" , &T::setc)
.add_method("setzc" , &T::setzc)
.add_method("setd" , &T::setd)
.add_method("setzd" , &T::setzd)
.add_method("setk" , &T::setk)
.add_method("setzk" , &T::setzk)
.add_method("setl" , &T::setl)
.add_method("setzl" , &T::setzl)
.add_method("setf" , &T::setf)
.add_method("setq" , &T::setq)
.add_method("setkci" , &T::setkci)
.add_method("setkic" , &T::setkic)
.add_method("setkoi" , &T::setkoi)
.add_method("setkio" , &T::setkio)
.add_method("setkii2" , &T::setkii2)
.add_method("setki2i" , &T::setki2i)
.add_method("set_log_C0Gate" , &T::set_log_C0Gate)
.add_method("set_log_C1Gate" , &T::set_log_C1Gate)
.add_method("set_log_C2Gate" , &T::set_log_C2Gate)
.add_method("set_log_C3Gate" , &T::set_log_C3Gate)
.add_method("set_log_C4Gate" , &T::set_log_C4Gate)
.add_method("set_log_C5Gate" , &T::set_log_C5Gate)
.add_method("set_log_I0Gate" , &T::set_log_I0Gate)
.add_method("set_log_I1Gate" , &T::set_log_I1Gate)
.add_method("set_log_I2Gate" , &T::set_log_I2Gate)
.add_method("set_log_I3Gate" , &T::set_log_I3Gate)
.add_method("set_log_I4Gate" , &T::set_log_I4Gate)
.add_method("set_log_I5Gate" , &T::set_log_I5Gate)
.add_method("set_log_OGate" , &T::set_log_OGate)
.add_method("set_log_I6Gate" , &T::set_log_I6Gate)
.add_method("set_log_I7Gate" , &T::set_log_I7Gate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "Kv4_csiosi_converted_standard_UG", tag); 
} 
 
 
#include "../converted/myseclamp_converted_standard_UG.h" 
 
 
{ 
	 typedef myseclamp_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("myseclamp_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setrs" , &T::setrs)
.add_method("setdur1" , &T::setdur1)
.add_method("setdur2" , &T::setdur2)
.add_method("setdur3" , &T::setdur3)
.add_method("setdur4" , &T::setdur4)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "myseclamp_converted_standard_UG", tag); 
} 
 
 
#include "kfast_converted_standard_UG.h" 
 
 
{ 
	 typedef kfast_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kfast_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgkbar" , &T::setgkbar)
.add_method("setvhalfn" , &T::setvhalfn)
.add_method("setvhalfl" , &T::setvhalfl)
.add_method("setkn" , &T::setkn)
.add_method("setkl" , &T::setkl)
.add_method("setq10" , &T::setq10)
.add_method("setqq" , &T::setqq)
.add_method("settq" , &T::settq)
.add_method("setek" , &T::setek)
.add_method("set_log_nGate" , &T::set_log_nGate)
.add_method("set_log_lGate" , &T::set_log_lGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kfast_converted_standard_UG", tag); 
} 
 
 
#include "kslow_converted_standard_UG.h" 
 
 
{ 
	 typedef kslow_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kslow_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgkbar" , &T::setgkbar)
.add_method("setvhalfn" , &T::setvhalfn)
.add_method("setvhalfl" , &T::setvhalfl)
.add_method("setkn" , &T::setkn)
.add_method("setkl" , &T::setkl)
.add_method("setq10" , &T::setq10)
.add_method("settmp" , &T::settmp)
.add_method("setek" , &T::setek)
.add_method("set_log_nGate" , &T::set_log_nGate)
.add_method("set_log_lGate" , &T::set_log_lGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kslow_converted_standard_UG", tag); 
} 
 
 
 
#include "kadist_converted_standard_UG.h" 
 
 
{ 
	 typedef kadist_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kadist_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgkabar" , &T::setgkabar)
.add_method("setvhalfn" , &T::setvhalfn)
.add_method("setvhalfl" , &T::setvhalfl)
.add_method("seta0l" , &T::seta0l)
.add_method("seta0n" , &T::seta0n)
.add_method("setzetan" , &T::setzetan)
.add_method("setzetal" , &T::setzetal)
.add_method("setgmn" , &T::setgmn)
.add_method("setgml" , &T::setgml)
.add_method("setlmin" , &T::setlmin)
.add_method("setnmin" , &T::setnmin)
.add_method("setpw" , &T::setpw)
.add_method("settq" , &T::settq)
.add_method("setqq" , &T::setqq)
.add_method("setq10" , &T::setq10)
.add_method("setqtl" , &T::setqtl)
.add_method("setek" , &T::setek)
.add_method("set_log_nGate" , &T::set_log_nGate)
.add_method("set_log_lGate" , &T::set_log_lGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kadist_converted_standard_UG", tag); 
} 
 
 
#include "kaprox_converted_standard_UG.h" 
 
 
{ 
	 typedef kaprox_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kaprox_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgkabar" , &T::setgkabar)
.add_method("setvhalfn" , &T::setvhalfn)
.add_method("setvhalfl" , &T::setvhalfl)
.add_method("seta0l" , &T::seta0l)
.add_method("seta0n" , &T::seta0n)
.add_method("setzetan" , &T::setzetan)
.add_method("setzetal" , &T::setzetal)
.add_method("setgmn" , &T::setgmn)
.add_method("setgml" , &T::setgml)
.add_method("setlmin" , &T::setlmin)
.add_method("setnmin" , &T::setnmin)
.add_method("setpw" , &T::setpw)
.add_method("settq" , &T::settq)
.add_method("setqq" , &T::setqq)
.add_method("setq10" , &T::setq10)
.add_method("setqtl" , &T::setqtl)
.add_method("setek" , &T::setek)
.add_method("set_log_nGate" , &T::set_log_nGate)
.add_method("set_log_lGate" , &T::set_log_lGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kaprox_converted_standard_UG", tag); 
} 
 
 
#include "kdrca1_converted_standard_UG.h" 
 
 
{ 
	 typedef kdrca1_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kdrca1_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setek" , &T::setek)
.add_method("setgkdrbar" , &T::setgkdrbar)
.add_method("setvhalfn" , &T::setvhalfn)
.add_method("seta0n" , &T::seta0n)
.add_method("setzetan" , &T::setzetan)
.add_method("setgmn" , &T::setgmn)
.add_method("setnmax" , &T::setnmax)
.add_method("setq10" , &T::setq10)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kdrca1_converted_standard_UG", tag); 
} 
 
 
#include "na3n_converted_standard_UG.h" 
 
 
{ 
	 typedef na3n_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("na3n_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("settha" , &T::settha)
.add_method("setqa" , &T::setqa)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("setthi1" , &T::setthi1)
.add_method("setthi2" , &T::setthi2)
.add_method("setqd" , &T::setqd)
.add_method("setqg" , &T::setqg)
.add_method("setmmin" , &T::setmmin)
.add_method("sethmin" , &T::sethmin)
.add_method("setq10" , &T::setq10)
.add_method("setRg" , &T::setRg)
.add_method("setRd" , &T::setRd)
.add_method("setqq" , &T::setqq)
.add_method("settq" , &T::settq)
.add_method("setthinf" , &T::setthinf)
.add_method("setqinf" , &T::setqinf)
.add_method("setvhalfs" , &T::setvhalfs)
.add_method("seta0s" , &T::seta0s)
.add_method("setzetas" , &T::setzetas)
.add_method("setgms" , &T::setgms)
.add_method("setvvh" , &T::setvvh)
.add_method("setvvs" , &T::setvvs)
.add_method("setar" , &T::setar)
.add_method("setena" , &T::setena)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "na3n_converted_standard_UG", tag); 
} 
 
 
#include "naxn_converted_standard_UG.h" 
 
 
{ 
	 typedef naxn_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("naxn_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("settha" , &T::settha)
.add_method("setqa" , &T::setqa)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("setthi1" , &T::setthi1)
.add_method("setthi2" , &T::setthi2)
.add_method("setqd" , &T::setqd)
.add_method("setqg" , &T::setqg)
.add_method("setmmin" , &T::setmmin)
.add_method("sethmin" , &T::sethmin)
.add_method("setq10" , &T::setq10)
.add_method("setRg" , &T::setRg)
.add_method("setRd" , &T::setRd)
.add_method("setthinf" , &T::setthinf)
.add_method("setqinf" , &T::setqinf)
.add_method("setena" , &T::setena)
.add_method("set_log_mGate" , &T::set_log_mGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "naxn_converted_standard_UG", tag); 
} 
 

 
 
#include "kdr_converted_standard_UG.h" 
 
 
{ 
	 typedef kdr_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kdr_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setek" , &T::setek)
.add_method("setgkdrbar" , &T::setgkdrbar)
.add_method("setvhalfn" , &T::setvhalfn)
.add_method("setvhalfl" , &T::setvhalfl)
.add_method("seta0l" , &T::seta0l)
.add_method("seta0n" , &T::seta0n)
.add_method("setzetan" , &T::setzetan)
.add_method("setzetal" , &T::setzetal)
.add_method("setgmn" , &T::setgmn)
.add_method("setgml" , &T::setgml)
.add_method("setnmax" , &T::setnmax)
.add_method("set_log_nGate" , &T::set_log_nGate)
.add_method("set_log_lGate" , &T::set_log_lGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kdr_converted_standard_UG", tag); 
} 
 
 
#include "namir_converted_standard_UG.h" 
 
 
{ 
	 typedef namir_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("namir_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setena" , &T::setena)
.add_method("setgnabar" , &T::setgnabar)
.add_method("setvhalf" , &T::setvhalf)
.add_method("setzeta" , &T::setzeta)
.add_method("setvhalfn" , &T::setvhalfn)
.add_method("setvhalfr" , &T::setvhalfr)
.add_method("setvhalfl" , &T::setvhalfl)
.add_method("seta0l" , &T::seta0l)
.add_method("seta0n" , &T::seta0n)
.add_method("seta0r" , &T::seta0r)
.add_method("setzetan" , &T::setzetan)
.add_method("setzetar" , &T::setzetar)
.add_method("setzetal" , &T::setzetal)
.add_method("setgmn" , &T::setgmn)
.add_method("setgmr" , &T::setgmr)
.add_method("setgml" , &T::setgml)
.add_method("setlmax" , &T::setlmax)
.add_method("setnmax" , &T::setnmax)
.add_method("setrmin" , &T::setrmin)
.add_method("setvvh" , &T::setvvh)
.add_method("setvvs" , &T::setvvs)
.add_method("setb" , &T::setb)
.add_method("set_log_nGate" , &T::set_log_nGate)
.add_method("set_log_lGate" , &T::set_log_lGate)
.add_method("set_log_rGate" , &T::set_log_rGate)
.add_method("set_log_frGate" , &T::set_log_frGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "namir_converted_standard_UG", tag); 
} 
 
 
#include "iq_converted_standard_UG.h" 
 
 
{ 
	 typedef iq_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("iq_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("seterevq" , &T::seterevq)
.add_method("setgqbar" , &T::setgqbar)
.add_method("setvhalf" , &T::setvhalf)
.add_method("seta0" , &T::seta0)
.add_method("setzeta" , &T::setzeta)
.add_method("setgq" , &T::setgq)
.add_method("set_log_qGate" , &T::set_log_qGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "iq_converted_standard_UG", tag); 
} 
 
 
#include "borgkm_converted_standard_UG.h" 
 
 
{ 
	 typedef borgkm_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("borgkm_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setcai" , &T::setcai)
.add_method("setek" , &T::setek)
.add_method("settm" , &T::settm)
.add_method("setgkmbar" , &T::setgkmbar)
.add_method("setvhalf" , &T::setvhalf)
.add_method("seta0" , &T::seta0)
.add_method("setzeta" , &T::setzeta)
.add_method("setgm" , &T::setgm)
.add_method("set_log_mGate" , &T::set_log_mGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "borgkm_converted_standard_UG", tag); 
} 
 
 
#include "kf_converted_standard_UG.h" 
 
 
{ 
	 typedef kf_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("kf_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setek" , &T::setek)
.add_method("setA_anF" , &T::setA_anF)
.add_method("setB_anF" , &T::setB_anF)
.add_method("setC_anF" , &T::setC_anF)
.add_method("setA_bnF" , &T::setA_bnF)
.add_method("setB_bnF" , &T::setB_bnF)
.add_method("setC_bnF" , &T::setC_bnF)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "kf_converted_standard_UG", tag); 
} 
 
 
#include "ks_converted_standard_UG.h" 
 
 
{ 
	 typedef ks_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("ks_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setek" , &T::setek)
.add_method("setA_anS" , &T::setA_anS)
.add_method("setB_anS" , &T::setB_anS)
.add_method("setC_anS" , &T::setC_anS)
.add_method("setA_bnS" , &T::setA_bnS)
.add_method("setB_bnS" , &T::setB_bnS)
.add_method("setC_bnS" , &T::setC_bnS)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "ks_converted_standard_UG", tag); 
} 
 
 
#include "nap_converted_standard_UG.h" 
 
 
{ 
	 typedef nap_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("nap_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setena" , &T::setena)
.add_method("setA_amp" , &T::setA_amp)
.add_method("setB_amp" , &T::setB_amp)
.add_method("setC_amp" , &T::setC_amp)
.add_method("setA_bmp" , &T::setA_bmp)
.add_method("setB_bmp" , &T::setB_bmp)
.add_method("setC_bmp" , &T::setC_bmp)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "nap_converted_standard_UG", tag); 
} 
 
 
#include "nattxs_converted_standard_UG.h" 
 
 
{ 
	 typedef nattxs_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("nattxs_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setena" , &T::setena)
.add_method("setA_am" , &T::setA_am)
.add_method("setB_am" , &T::setB_am)
.add_method("setC_am" , &T::setC_am)
.add_method("setA_ah" , &T::setA_ah)
.add_method("setB_ah" , &T::setB_ah)
.add_method("setC_ah" , &T::setC_ah)
.add_method("setA_bm" , &T::setA_bm)
.add_method("setB_bm" , &T::setB_bm)
.add_method("setC_bm" , &T::setC_bm)
.add_method("setA_bh" , &T::setA_bh)
.add_method("setB_bh" , &T::setB_bh)
.add_method("setC_bh" , &T::setC_bh)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "nattxs_converted_standard_UG", tag); 
} 
 
 
#include "nav1p8_converted_standard_UG.h" 
 
 
{ 
	 typedef nav1p8_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("nav1p8_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setena" , &T::setena)
.add_method("setA_am8" , &T::setA_am8)
.add_method("setB_am8" , &T::setB_am8)
.add_method("setC_am8" , &T::setC_am8)
.add_method("setA_ah8" , &T::setA_ah8)
.add_method("setB_ah8" , &T::setB_ah8)
.add_method("setC_ah8" , &T::setC_ah8)
.add_method("setA_bm8" , &T::setA_bm8)
.add_method("setB_bm8" , &T::setB_bm8)
.add_method("setC_bm8" , &T::setC_bm8)
.add_method("setA_bh8" , &T::setA_bh8)
.add_method("setB_bh8" , &T::setB_bh8)
.add_method("setC_bh8" , &T::setC_bh8)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "nav1p8_converted_standard_UG", tag); 
} 
 
 
#include "nav1p9_converted_standard_UG.h" 
 
 
{ 
	 typedef nav1p9_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("nav1p9_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("setena" , &T::setena)
.add_method("setA_am9" , &T::setA_am9)
.add_method("setB_am9" , &T::setB_am9)
.add_method("setC_am9" , &T::setC_am9)
.add_method("setA_ah9" , &T::setA_ah9)
.add_method("setB_ah9" , &T::setB_ah9)
.add_method("setC_ah9" , &T::setC_ah9)
.add_method("setA_bm9" , &T::setA_bm9)
.add_method("setB_bm9" , &T::setB_bm9)
.add_method("setC_bm9" , &T::setC_bm9)
.add_method("setA_bh9" , &T::setA_bh9)
.add_method("setB_bh9" , &T::setB_bh9)
.add_method("setC_bh9" , &T::setC_bh9)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "nav1p9_converted_standard_UG", tag); 
} 
 
 
#include "hh_Cp_scaled_converted_standard_UG.h" 
 
 
{ 
	 typedef hh_Cp_scaled_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("hh_Cp_scaled_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setgl" , &T::setgl)
.add_method("setel" , &T::setel)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "hh_Cp_scaled_converted_standard_UG", tag); 
} 
 
 
#include "hh_Cs_scaled_converted_standard_UG.h" 
 
 
{ 
	 typedef hh_Cs_scaled_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("hh_Cs_scaled_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setgl" , &T::setgl)
.add_method("setel" , &T::setel)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "hh_Cs_scaled_converted_standard_UG", tag); 
} 
 
 
#include "hh_Ca_scaled_converted_standard_UG.h" 
 
 
{ 
	 typedef hh_Ca_scaled_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("hh_Ca_scaled_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setgl" , &T::setgl)
.add_method("setel" , &T::setel)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "hh_Ca_scaled_converted_standard_UG", tag); 
} 
 
 
#include "Kv1_converted_standard_UG.h" 
 
 
{ 
	 typedef Kv1_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("Kv1_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgbar" , &T::setgbar)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "Kv1_converted_standard_UG", tag); 
} 
 
 
#include "nadifl_converted_standard_UG.h" 
 
 
{ 
	 typedef nadifl_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("nadifl_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setD" , &T::setD)
.add_method("set_log_naiGate" , &T::set_log_naiGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "nadifl_converted_standard_UG", tag); 
} 
 
 
#include "AXNODE75_converted_standard_UG.h" 
 
 
{ 
	 typedef AXNODE75_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("AXNODE75_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnapbar" , &T::setgnapbar)
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setgl" , &T::setgl)
.add_method("setena" , &T::setena)
.add_method("setek" , &T::setek)
.add_method("setel" , &T::setel)
.add_method("setvshift" , &T::setvshift)
.add_method("setvtraub" , &T::setvtraub)
.add_method("setampA" , &T::setampA)
.add_method("setampB" , &T::setampB)
.add_method("setampC" , &T::setampC)
.add_method("setbmpA" , &T::setbmpA)
.add_method("setbmpB" , &T::setbmpB)
.add_method("setbmpC" , &T::setbmpC)
.add_method("setamA" , &T::setamA)
.add_method("setamB" , &T::setamB)
.add_method("setamC" , &T::setamC)
.add_method("setbmA" , &T::setbmA)
.add_method("setbmB" , &T::setbmB)
.add_method("setbmC" , &T::setbmC)
.add_method("setahA" , &T::setahA)
.add_method("setahB" , &T::setahB)
.add_method("setahC" , &T::setahC)
.add_method("setbhA" , &T::setbhA)
.add_method("setbhB" , &T::setbhB)
.add_method("setbhC" , &T::setbhC)
.add_method("setasA" , &T::setasA)
.add_method("setasB" , &T::setasB)
.add_method("setasC" , &T::setasC)
.add_method("setbsA" , &T::setbsA)
.add_method("setbsB" , &T::setbsB)
.add_method("setbsC" , &T::setbsC)
.add_method("set_log_mpGate" , &T::set_log_mpGate)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
.add_method("set_log_sGate" , &T::set_log_sGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "AXNODE75_converted_standard_UG", tag); 
} 
 
 
#include "SK_converted_standard_UG.h" 
 
 
{ 
	 typedef SK_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("SK_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgk" , &T::setgk)
.add_method("setisKCa" , &T::setisKCa)
.add_method("setek" , &T::setek)
.add_method("setki" , &T::setki)
.add_method("setcai" , &T::setcai)
.add_method("set_log_wGate" , &T::set_log_wGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "SK_converted_standard_UG", tag); 
} 
 
 
#include "KDRf_converted_standard_UG.h" 
 
 
{ 
	 typedef KDRf_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("KDRf_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgk" , &T::setgk)
.add_method("setek" , &T::setek)
.add_method("setki" , &T::setki)
.add_method("set_log_pfastGate" , &T::set_log_pfastGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "KDRf_converted_standard_UG", tag); 
} 
 
 
#include "KDRs_converted_standard_UG.h" 
 
 
{ 
	 typedef KDRs_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("KDRs_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgk" , &T::setgk)
.add_method("setek" , &T::setek)
.add_method("setki" , &T::setki)
.add_method("set_log_pslowdGate" , &T::set_log_pslowdGate)
.add_method("set_log_pslowiGate" , &T::set_log_pslowiGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "KDRs_converted_standard_UG", tag); 
} 
 
 
#include "HCN_converted_standard_UG.h" 
 
 
{ 
	 typedef HCN_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("HCN_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgh" , &T::setgh)
.add_method("seteih" , &T::seteih)
.add_method("set_log_fGate" , &T::set_log_fGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "HCN_converted_standard_UG", tag); 
} 
 
 
#include "HVA_converted_standard_UG.h" 
 
 
{ 
	 typedef HVA_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("HVA_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgcaL" , &T::setgcaL)
.add_method("setgcaN" , &T::setgcaN)
.add_method("setiNCa" , &T::setiNCa)
.add_method("setiLCa" , &T::setiLCa)
.add_method("setinactLtau" , &T::setinactLtau)
.add_method("setinactLmax" , &T::setinactLmax)
.add_method("seteca" , &T::seteca)
.add_method("setcai" , &T::setcai)
.add_method("setcao" , &T::setcao)
.add_method("set_log_qGate" , &T::set_log_qGate)
.add_method("set_log_uGate" , &T::set_log_uGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "HVA_converted_standard_UG", tag); 
} 
 
 
/*#include "Na_converted_standard_UG.h"
 
 
{ 
	 typedef Na_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("Na_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgna" , &T::setgna)
.add_method("setrest" , &T::setrest)
.add_method("setena" , &T::setena)
.add_method("setnai" , &T::setnai)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "Na_converted_standard_UG", tag); 
} */
 
 
#include "NaL_converted_standard_UG.h" 
 
 
{ 
	 typedef NaL_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("NaL_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgna" , &T::setgna)
.add_method("setinaL" , &T::setinaL)
.add_method("setena" , &T::setena)
.add_method("setnai" , &T::setnai)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "NaL_converted_standard_UG", tag); 
} 
 
 
#include "gpi_converted_standard_UG.h" 
 
 
{ 
	 typedef gpi_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("gpi_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgpas" , &T::setgpas)
.add_method("setepas" , &T::setepas)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "gpi_converted_standard_UG", tag); 
} 
 
 
#include "cacum_converted_standard_UG.h" 
 
 
{ 
	 typedef cacum_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("cacum_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setdepth" , &T::setdepth)
.add_method("settau" , &T::settau)
.add_method("setcai0" , &T::setcai0)
.add_method("setcai0_ca_ion" , &T::setcai0_ca_ion)
.add_method("set_log_caiGate" , &T::set_log_caiGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "cacum_converted_standard_UG", tag); 
} 
 
 
#include "PARAK75_converted_standard_UG.h" 
 
 
{ 
	 typedef PARAK75_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("PARAK75_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgkbar" , &T::setgkbar)
.add_method("setek" , &T::setek)
.add_method("setvshift" , &T::setvshift)
.add_method("setanA" , &T::setanA)
.add_method("setanB" , &T::setanB)
.add_method("setanC" , &T::setanC)
.add_method("setbnA" , &T::setbnA)
.add_method("setbnB" , &T::setbnB)
.add_method("setbnC" , &T::setbnC)
.add_method("set_log_nGate" , &T::set_log_nGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "PARAK75_converted_standard_UG", tag); 
} 
 
 

 
 
#include "Na2_converted_standard_UG.h" 
 
 
{ 
	 typedef Na2_converted_standard_UG<TDomain> T; 
	 typedef ICableMembraneTransport<TDomain> TBase; 
	 string name = string("Na2_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgna" , &T::setgna)
.add_method("setrest" , &T::setrest)
.add_method("setena" , &T::setena)
.add_method("setnai" , &T::setnai)
.add_method("set_log_SGate" , &T::set_log_SGate)
.add_method("set_log_mGate" , &T::set_log_mGate)
.add_method("set_log_hGate" , &T::set_log_hGate)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "Na2_converted_standard_UG", tag); 
} 
 
 
