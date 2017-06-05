#include "nax_g01_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable_neuron { 
 
 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::trap0(double v, double th, double a, double q) 
{ 
double trap0; 
	if (fabs((v-th)*1(/mV)) > 1e-6) {
	        return  1(ms/mV)*a * (v - th) / (1 - exp(-(v - th)/q)); 
	} else {
	        return  1(ms/mV)*a * q; 
 	}
}	

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::ce_obj_available()  
{  
	init_attachments();  
 	F = m_pCE->F; 
 R = m_pCE->R; 
 K = m_pCE->temperature(); 
 celsius = m_pCE->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getgbar() 
{ 
return gbar; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::gettha() 
{ 
return tha; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getqa() 
{ 
return qa; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getRa() 
{ 
return Ra; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getRb() 
{ 
return Rb; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getthi1() 
{ 
return thi1; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getthi2() 
{ 
return thi2; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getqd() 
{ 
return qd; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getqg() 
{ 
return qg; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getmmin() 
{ 
return mmin; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::gethmin() 
{ 
return hmin; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getq10() 
{ 
return q10; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getRg() 
{ 
return Rg; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getRd() 
{ 
return Rd; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getthinf() 
{ 
return thinf; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getqinf() 
{ 
return qinf; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getena() 
{ 
return ena; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::getmscale() 
{ 
return mscale; 
} 
template<typename TDomain> 
double nax_g01_converted_standard_UG<TDomain>::gethscale() 
{ 
return hscale; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setgbar(double val) 
{ 
gbar = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::settha(double val) 
{ 
tha = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setqa(double val) 
{ 
qa = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setRa(double val) 
{ 
Ra = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setRb(double val) 
{ 
Rb = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setthi1(double val) 
{ 
thi1 = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setthi2(double val) 
{ 
thi2 = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setqd(double val) 
{ 
qd = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setqg(double val) 
{ 
qg = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setmmin(double val) 
{ 
mmin = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::sethmin(double val) 
{ 
hmin = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setq10(double val) 
{ 
q10 = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setRg(double val) 
{ 
Rg = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setRd(double val) 
{ 
Rd = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setthinf(double val) 
{ 
thinf = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setqinf(double val) 
{ 
qinf = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setena(double val) 
{ 
ena = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::setmscale(double val) 
{ 
mscale = val; 
} 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::sethscale(double val) 
{ 
hscale = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for nax_g01_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for nax_g01_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> nax_g01_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 try { ssGrp = SubsetGroup(m_pCE->approx_space()->domain()->subset_handler(), this->m_vSubset);} 
	 UG_CATCH_THROW("Subset group creation failed."); 
 
	 itType iter; 
	 number bestDistSq, distSq; 
	 Vertex* bestVrt; 
 
	 // Iterate only if there is one Gtting needed 
	 if (m_log_mGate || m_log_hGate )
	 { 
	 	 // iterating over all elements 
	 	 for (size_t si=0; si < ssGrp.size(); si++) 
	 	 { 
	 	 	 itType iterBegin = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]); 
	 	 	 itType iterEnd = m_pCE->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
	 	 	 const position_accesor_type& aaPos = m_pCE->approx_space()->domain()->position_accessor(); 
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
	 	 if (m_log_mGate == true) 
	 	 	 GatingAccesors.push_back(this->aamGate[bestVrt]); 
	 	 if (m_log_hGate == true) 
	 	 	 GatingAccesors.push_back(this->aahGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void nax_g01_converted_standard_UG<TDomain>::set_log_mGate(bool bLogmGate) { m_log_mGate = bLogmGate; }
template<typename TDomain> void nax_g01_converted_standard_UG<TDomain>::set_log_hGate(bool bLoghGate) { m_log_hGate = bLoghGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
//get celsius and time
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number celsius = m_pCE->temperature_celsius(); 
number dt = m_pCE->time(); 
// make preparing vor getting values of every edge 
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number na = vrt_values[CableEquation<TDomain>::_na_]; 

 
number vm = v;
double           a, b, qt; 
qt= pow(q10 , ((celsius-24)/10(degC))); 
	a = trap0(vm,tha,Ra,qa); 
	b = trap0(-vm,-tha,Rb,qa); 
double 	mtau = 1(ms)/(a+b)/qt; 
        if (mtau<mmin) {mtau=mmin;};; 
        mtau = mtau/mscale; 
double 	minf = a/(a+b); 
double         mexp = 1 - exp(-dt/mtau); 
	a = trap0(vm,thi1,Rd,qd); 
	b = trap0(-vm,-thi2,Rg,qg); 
double 	htau =  1(ms)/(a+b)/qt; 
        if (htau<hmin) {htau=hmin;};; 
        htau = htau/hscale; 
double 	hinf = 1/(1+exp((vm-thinf)/qinf)); 
double         hexp = 1 - exp(-dt/htau); 
aamGate[vrt] =minf; 
aahGate[vrt] =hinf; 
double thegna =  gbar*m*m*m*h; 
double ina =  thegna * (v - ena); 
}  
 
 
 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number celsius = m_pCE->temperature_celsius(); 
 number FARADAY = m_pCE->F; 
 number dt = newTime - m_pCE->time(); 
number v = vrt_values[CableEquation<TDomain>::_v_]; 
number na = vrt_values[CableEquation<TDomain>::_na_]; 

 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 

 
 
number vm = v;
double           a, b, qt; 
qt= pow(q10 , ((celsius-24)/10(degC))); 
	a = trap0(vm,tha,Ra,qa); 
	b = trap0(-vm,-tha,Rb,qa); 
double 	mtau = 1(ms)/(a+b)/qt; 
        if (mtau<mmin) {mtau=mmin;};; 
        mtau = mtau/mscale; 
double 	minf = a/(a+b); 
double         mexp = 1 - exp(-dt/mtau); 
	a = trap0(vm,thi1,Rd,qd); 
	b = trap0(-vm,-thi2,Rg,qg); 
double 	htau =  1(ms)/(a+b)/qt; 
        if (htau<hmin) {htau=hmin;};; 
        htau = htau/hscale; 
double 	hinf = 1/(1+exp((vm-thinf)/qinf)); 
double         hexp = 1 - exp(-dt/htau); 
         m  +=  (minf - m)/mtau*dt; 
; 
 
         h  +=  (hinf - h)/htau*dt; 
; 
 

 
 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pCE->temperature(); 
m_R = m_pCE->R; 
m_F = m_pCE->F; 
 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number na = vrt_values[m_pCE->_na_]; 
number v =  vrt_values[m_pCE->_v_]; 
 
 
number t = m_pCE->time(); 
 
 

 
 
const number helpV = 1e3*(m_pCE->R*m_pCE->temperature())/m_pCE->F; 
 
 
number thegna = gbar*m*m*m*h; 

 
 
outCurrentValues.push_back( thegna * (v - ena)); 
} 
 
 
template<typename TDomain> 
void nax_g01_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(CableEquation<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class nax_g01_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class nax_g01_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class nax_g01_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable_neuron
} // namespace ug


