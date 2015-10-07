#include "NMDA_Mg_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
 	F = m_pVMDisc->F; 
 R = m_pVMDisc->R; 
 K = m_pVMDisc->temperature(); 
 celsius = m_pVMDisc->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getErev() 
{ 
return Erev; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getgmax() 
{ 
return gmax; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getmg() 
{ 
return mg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getvmin() 
{ 
return vmin; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getvmax() 
{ 
return vmax; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getvalence() 
{ 
return valence; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getmemb_fraction() 
{ 
return memb_fraction; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRb() 
{ 
return Rb; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRu() 
{ 
return Ru; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRo() 
{ 
return Ro; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRc() 
{ 
return Rc; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRd1() 
{ 
return Rd1; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRr1() 
{ 
return Rr1; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRd2() 
{ 
return Rd2; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRr2() 
{ 
return Rr2; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmb() 
{ 
return Rmb; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmu() 
{ 
return Rmu; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmc1b() 
{ 
return Rmc1b; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmc1u() 
{ 
return Rmc1u; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmc2b() 
{ 
return Rmc2b; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmc2u() 
{ 
return Rmc2u; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmd1b() 
{ 
return Rmd1b; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmd1u() 
{ 
return Rmd1u; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmd2b() 
{ 
return Rmd2b; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRmd2u() 
{ 
return Rmd2u; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRbMg() 
{ 
return RbMg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRuMg() 
{ 
return RuMg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRoMg() 
{ 
return RoMg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRcMg() 
{ 
return RcMg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRd1Mg() 
{ 
return Rd1Mg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRr1Mg() 
{ 
return Rr1Mg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRd2Mg() 
{ 
return Rd2Mg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getRr2Mg() 
{ 
return Rr2Mg; 
} 
template<typename TDomain> 
double NMDA_Mg_converted_standard_UG<TDomain>::getC() 
{ 
return C; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setErev(double val) 
{ 
Erev = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setgmax(double val) 
{ 
gmax = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setmg(double val) 
{ 
mg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setvmin(double val) 
{ 
vmin = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setvmax(double val) 
{ 
vmax = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setvalence(double val) 
{ 
valence = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setmemb_fraction(double val) 
{ 
memb_fraction = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRb(double val) 
{ 
Rb = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRu(double val) 
{ 
Ru = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRo(double val) 
{ 
Ro = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRc(double val) 
{ 
Rc = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRd1(double val) 
{ 
Rd1 = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRr1(double val) 
{ 
Rr1 = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRd2(double val) 
{ 
Rd2 = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRr2(double val) 
{ 
Rr2 = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmb(double val) 
{ 
Rmb = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmu(double val) 
{ 
Rmu = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmc1b(double val) 
{ 
Rmc1b = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmc1u(double val) 
{ 
Rmc1u = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmc2b(double val) 
{ 
Rmc2b = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmc2u(double val) 
{ 
Rmc2u = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmd1b(double val) 
{ 
Rmd1b = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmd1u(double val) 
{ 
Rmd1u = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmd2b(double val) 
{ 
Rmd2b = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRmd2u(double val) 
{ 
Rmd2u = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRbMg(double val) 
{ 
RbMg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRuMg(double val) 
{ 
RuMg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRoMg(double val) 
{ 
RoMg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRcMg(double val) 
{ 
RcMg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRd1Mg(double val) 
{ 
Rd1Mg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRr1Mg(double val) 
{ 
Rr1Mg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRd2Mg(double val) 
{ 
Rd2Mg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setRr2Mg(double val) 
{ 
Rr2Mg = val; 
} 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::setC(double val) 
{ 
C = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->UGate)) 
UG_THROW("Attachment necessary (UGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->UGate); 
this->aaUGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->UGate); 
 
if (spGrid->has_vertex_attachment(this->ClGate)) 
UG_THROW("Attachment necessary (ClGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->ClGate); 
this->aaClGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->ClGate); 
 
if (spGrid->has_vertex_attachment(this->D1Gate)) 
UG_THROW("Attachment necessary (D1Gate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->D1Gate); 
this->aaD1Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->D1Gate); 
 
if (spGrid->has_vertex_attachment(this->D2Gate)) 
UG_THROW("Attachment necessary (D2Gate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->D2Gate); 
this->aaD2Gate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->D2Gate); 
 
if (spGrid->has_vertex_attachment(this->OGate)) 
UG_THROW("Attachment necessary (OGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->OGate); 
this->aaOGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->OGate); 
 
if (spGrid->has_vertex_attachment(this->UMgGate)) 
UG_THROW("Attachment necessary (UMgGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->UMgGate); 
this->aaUMgGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->UMgGate); 
 
if (spGrid->has_vertex_attachment(this->ClMgGate)) 
UG_THROW("Attachment necessary (ClMgGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->ClMgGate); 
this->aaClMgGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->ClMgGate); 
 
if (spGrid->has_vertex_attachment(this->D1MgGate)) 
UG_THROW("Attachment necessary (D1MgGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->D1MgGate); 
this->aaD1MgGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->D1MgGate); 
 
if (spGrid->has_vertex_attachment(this->D2MgGate)) 
UG_THROW("Attachment necessary (D2MgGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->D2MgGate); 
this->aaD2MgGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->D2MgGate); 
 
if (spGrid->has_vertex_attachment(this->OMgGate)) 
UG_THROW("Attachment necessary (OMgGate) for NMDA_Mg_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->OMgGate); 
this->aaOMgGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->OMgGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> NMDA_Mg_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
{ 
	 //var for output 
	 std::vector<number> GatingAccesors; 
 
	 typedef ug::MathVector<TDomain::dim> position_type; 
 
	 position_type coord; 
 
	 if (coord.size()==1) 
	 	 coord[0]=x; 
	 if (coord.size()==2) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 } 
	 if (coord.size()==3) 
	 { 
	 	 coord[0] = x;
	 	 coord[1] = y;
	 	 coord[2] = z;
	 } 
	 //accesors 
	 typedef Attachment<position_type> position_attachment_type; 
	 typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accesor_type; 
 
	 // Definitions for Iteration over all Elements 
	 typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; 
	 SubsetGroup ssGrp; 
	 try { ssGrp = SubsetGroup(m_pVMDisc->approx_space()->domain()->subset_handler(), this->m_vSubset);} 
	 UG_CATCH_THROW("Subset group creation failed."); 
 
	 itType iter; 
	 number bestDistSq, distSq; 
	 Vertex* bestVrt; 
 
	 // Iterate only if there is one Gtting needed 
	 if (m_log_UGate || m_log_ClGate || m_log_D1Gate || m_log_D2Gate || m_log_OGate || m_log_UMgGate || m_log_ClMgGate || m_log_D1MgGate || m_log_D2MgGate || m_log_OMgGate )
	 { 
	 	 // iterating over all elements 
	 	 for (size_t si=0; si < ssGrp.size(); si++) 
	 	 { 
	 	 	 itType iterBegin = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]); 
	 	 	 itType iterEnd = m_pVMDisc->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
	 	 	 const position_accesor_type& aaPos = m_pVMDisc->approx_space()->domain()->position_accessor(); 
	 	 	 if (si==0) 
	 	 	 { 
	 	 	 	 bestVrt = *iterBegin; 
	 	 	 	 bestDistSq = VecDistanceSq(coord, aaPos[bestVrt]); 
	 	 	 } 
	 	 	 iter = iterBegin; 
	 	 	 iter++; 
	 	 	 while(iter != iterEnd) 
	 	 	 { 
	 	 	 	 distSq = VecDistanceSq(coord, aaPos[*iter]); 
	 	 	 	 { 
	 	 	 	 	 bestDistSq = distSq; 
	 	 	 	 	 bestVrt = *iter; 
	 	 	 	 } 
	 	 	 	 ++iter; 
	 	 	 } 
	 	 } 
	 	 if (m_log_UGate == true) 
	 	 	 GatingAccesors.push_back(this->aaUGate[bestVrt]); 
	 	 if (m_log_ClGate == true) 
	 	 	 GatingAccesors.push_back(this->aaClGate[bestVrt]); 
	 	 if (m_log_D1Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaD1Gate[bestVrt]); 
	 	 if (m_log_D2Gate == true) 
	 	 	 GatingAccesors.push_back(this->aaD2Gate[bestVrt]); 
	 	 if (m_log_OGate == true) 
	 	 	 GatingAccesors.push_back(this->aaOGate[bestVrt]); 
	 	 if (m_log_UMgGate == true) 
	 	 	 GatingAccesors.push_back(this->aaUMgGate[bestVrt]); 
	 	 if (m_log_ClMgGate == true) 
	 	 	 GatingAccesors.push_back(this->aaClMgGate[bestVrt]); 
	 	 if (m_log_D1MgGate == true) 
	 	 	 GatingAccesors.push_back(this->aaD1MgGate[bestVrt]); 
	 	 if (m_log_D2MgGate == true) 
	 	 	 GatingAccesors.push_back(this->aaD2MgGate[bestVrt]); 
	 	 if (m_log_OMgGate == true) 
	 	 	 GatingAccesors.push_back(this->aaOMgGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_UGate(bool bLogUGate) { m_log_UGate = bLogUGate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_ClGate(bool bLogClGate) { m_log_ClGate = bLogClGate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_D1Gate(bool bLogD1Gate) { m_log_D1Gate = bLogD1Gate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_D2Gate(bool bLogD2Gate) { m_log_D2Gate = bLogD2Gate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_OGate(bool bLogOGate) { m_log_OGate = bLogOGate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_UMgGate(bool bLogUMgGate) { m_log_UMgGate = bLogUMgGate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_ClMgGate(bool bLogClMgGate) { m_log_ClMgGate = bLogClMgGate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_D1MgGate(bool bLogD1MgGate) { m_log_D1MgGate = bLogD1MgGate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_D2MgGate(bool bLogD2MgGate) { m_log_D2MgGate = bLogD2MgGate; }
template<typename TDomain> void NMDA_Mg_converted_standard_UG<TDomain>::set_log_OMgGate(bool bLogOMgGate) { m_log_OMgGate = bLogOMgGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
number dt = m_pVMDisc->time(); 
// make preparing vor getting values of every edge 
number v = vrt_values[VMDisc<TDomain>::_v_]; 

 
aaUGate[vrt] = 1; 
}  
 
 
 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number celsius = m_pVMDisc->temperature_celsius(); 
 number FARADAY = m_pVMDisc->F; 
 number dt = newTime - m_pVMDisc->time(); 
number v = vrt_values[VMDisc<TDomain>::_v_]; 

 
double U = aaUGate[vrt]; 
double Cl = aaClGate[vrt]; 
double D1 = aaD1Gate[vrt]; 
double D2 = aaD2Gate[vrt]; 
double O = aaOGate[vrt]; 
double UMg = aaUMgGate[vrt]; 
double ClMg = aaClMgGate[vrt]; 
double D1Mg = aaD1MgGate[vrt]; 
double D2Mg = aaD2MgGate[vrt]; 
double OMg = aaOMgGate[vrt]; 

 
 

 
 
double 	rb 	= Rb 	* (1e3) * C;
double 	rbMg 	= RbMg 	* (1e3) * C;
double 	rmb 	= Rmb 	* mg * (1e3) * exp((v-40) * valence * memb_fraction /25);
double 	rmu 	= Rmu 	* exp((-1)*(v-40) * valence * (1-memb_fraction) /25);
double 	rmc1b 	= Rmc1b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25);
double 	rmc1u 	= Rmc1u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25);
double 	rmc2b 	= Rmc2b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25);
double 	rmc2u 	= Rmc2u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25);
double 	rmd1b 	= Rmd1b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25);
double 	rmd1u 	= Rmd1u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25);
double 	rmd2b 	= Rmd2b * mg * (1e3) * exp((v-40) * valence * memb_fraction /25);
double 	rmd2u 	= Rmd2u * exp((-1)*(v-40) * valence * (1-memb_fraction) /25);
 
 
 
U+=(-U*rb+Cl*Ru)*dt; 
Cl+=(U*rb+-Cl*Ru)*dt; 
Cl+=(-Cl*Ro+O*Rc)*dt; 
O+=(Cl*Ro+-O*Rc)*dt; 
Cl+=(-Cl*Rd1+D1*Rr1)*dt; 
D1+=(Cl*Rd1+-D1*Rr1)*dt; 
D1+=(-D1*Rd2+D2*Rr2)*dt; 
D2+=(D1*Rd2+-D2*Rr2)*dt; 
O+=(-O*rmb+OMg*rmu)*dt; 
OMg+=(O*rmb+-OMg*rmu)*dt; 
UMg+=(-UMg*rbMg+ClMg*RuMg)*dt; 
ClMg+=(UMg*rbMg+-ClMg*RuMg)*dt; 
ClMg+=(-ClMg*RoMg+OMg*RcMg)*dt; 
OMg+=(ClMg*RoMg+-OMg*RcMg)*dt; 
ClMg+=(-ClMg*Rd1Mg+D1Mg*Rr1Mg)*dt; 
D1Mg+=(ClMg*Rd1Mg+-D1Mg*Rr1Mg)*dt; 
D1Mg+=(-D1Mg*Rd2Mg+D2Mg*Rr2Mg)*dt; 
D2Mg+=(D1Mg*Rd2Mg+-D2Mg*Rr2Mg)*dt; 
U+=(-U*rmc1b+UMg*rmc1u)*dt; 
UMg+=(U*rmc1b+-UMg*rmc1u)*dt; 
Cl+=(-Cl*rmc2b+ClMg*rmc2u)*dt; 
ClMg+=(Cl*rmc2b+-ClMg*rmc2u)*dt; 
D1+=(-D1*rmd1b+D1Mg*rmd1u)*dt; 
D1Mg+=(D1*rmd1b+-D1Mg*rmd1u)*dt; 
D2+=(-D2*rmd2b+D2Mg*rmd2u)*dt; 
D2Mg+=(D2*rmd2b+-D2Mg*rmd2u)*dt; 
 
 
 
aaUGate[vrt] = U; 
aaClGate[vrt] = Cl; 
aaD1Gate[vrt] = D1; 
aaD2Gate[vrt] = D2; 
aaOGate[vrt] = O; 
aaUMgGate[vrt] = UMg; 
aaClMgGate[vrt] = ClMg; 
aaD1MgGate[vrt] = D1Mg; 
aaD2MgGate[vrt] = D2Mg; 
aaOMgGate[vrt] = OMg; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number U = aaUGate[ver]; 
number Cl = aaClGate[ver]; 
number D1 = aaD1Gate[ver]; 
number D2 = aaD2Gate[ver]; 
number O = aaOGate[ver]; 
number UMg = aaUMgGate[ver]; 
number ClMg = aaClMgGate[ver]; 
number D1Mg = aaD1MgGate[ver]; 
number D2Mg = aaD2MgGate[ver]; 
number OMg = aaOMgGate[ver]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 

 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
 
 
number g = gmax * O; 

 
 
outCurrentValues.push_back( (1e-6) * g * (v - Erev)); 
} 
 
 
template<typename TDomain> 
void NMDA_Mg_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class NMDA_Mg_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class NMDA_Mg_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class NMDA_Mg_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


