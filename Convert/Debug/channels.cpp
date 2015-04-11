
 
 
#include "ar_converted_standard_UG.h" 
 
 
{ 
	 typedef ar_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
 
 
#include "ca_converted_standard_UG.h" 
 
 
{ 
	 typedef ca_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "ca_converted_standard_UG", tag); 
} 
 
 
#include "cad_converted_standard_UG.h" 
 
 
{ 
	 typedef cad_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("cad_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setdepth" , &T::setdepth)
.add_method("settaur" , &T::settaur)
.add_method("setcainf" , &T::setcainf)
.add_method("setcai" , &T::setcai)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "cad_converted_standard_UG", tag); 
} 
 
 
#include "CaT_converted_standard_UG.h" 
 
 
{ 
	 typedef CaT_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "CaT_converted_standard_UG", tag); 
} 
 
 
#include "h_converted_standard_UG.h" 
 
 
{ 
	 typedef h_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setgl" , &T::setgl)
.add_method("setel" , &T::setel)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "hh_converted_standard_UG", tag); 
} 
 
 
#include "HH2_converted_standard_UG.h" 
 
 
{ 
	 typedef HH2_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
	 string name = string("HH2_converted_standard_UG").append(suffix); 
	 reg.add_class_<T, TBase >(name, grp) 
	 	 .template add_constructor<void (*)(const char*, const char*)>("Function(s)#Subset(s)") 
	 	 .template add_constructor<void (*)(const std::vector<std::string>&, const std::vector<std::string>&)>("Function(s)#Subset(s)") 
.add_method("setgnabar" , &T::setgnabar)
.add_method("setgkbar" , &T::setgkbar)
.add_method("setena" , &T::setena)
.add_method("setek" , &T::setek)
.add_method("setcelsius" , &T::setcelsius)
.add_method("setvtraub" , &T::setvtraub)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "HH2_converted_standard_UG", tag); 
} 
 
 
#include "inwardrect_converted_standard_UG.h" 
 
 
{ 
	 typedef inwardrect_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "inwardrect_converted_standard_UG", tag); 
} 
 
 
#include "kca_converted_standard_UG.h" 
 
 
{ 
	 typedef kca_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
.add_method("setgkbar" , &T::setgkbar)
.add_method("setmvhalf" , &T::setmvhalf)
.add_method("setmslope" , &T::setmslope)
.add_method("setmshift" , &T::setmshift)
.add_method("setqfact" , &T::setqfact)
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
.add_method("setgbar" , &T::setgbar)
.add_method("settha" , &T::settha)
.add_method("setqa" , &T::setqa)
.add_method("setRa" , &T::setRa)
.add_method("setRb" , &T::setRb)
.add_method("settemp" , &T::settemp)
.add_method("setq10" , &T::setq10)
.add_method("setvmin" , &T::setvmin)
.add_method("setvmax" , &T::setvmax)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "km_converted_standard_UG", tag); 
} 
 
 
#include "kv_converted_standard_UG.h" 
 
 
{ 
	 typedef kv_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "na_converted_standard_UG", tag); 
} 
 
 
#include "NMDA_Mg_converted_standard_UG.h" 
 
 
{ 
	 typedef NMDA_Mg_converted_standard_UG<TDomain> T; 
	 typedef IChannel<TDomain> TBase; 
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
.add_method("setg" , &T::setg)
.add_method("sete" , &T::sete)
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
.add_method("setdel" , &T::setdel)
.add_method("setdur" , &T::setdur)
.add_method("setamp" , &T::setamp)
	 	 .set_construct_as_smart_pointer(true); 
	 reg.add_class_to_group(name, "release_BMK_converted_standard_UG", tag); 
} 
 
 
 
 
