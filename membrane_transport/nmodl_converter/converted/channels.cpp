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
 
 
