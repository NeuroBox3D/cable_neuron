#include "hh_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
namespace cable { 
 
 
template<typename TDomain> 
double hh_converted_standard_UG<TDomain>::vtrap(double x, double y) 
{ 
double vtrap; 
        if (fabs(x/y) < 1e-6) {
                return  y*(1 - x/y/2); 
        }else{
                return  x/(exp(x/y) - 1); 
        }
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::vm_disc_available()  
{  
	init_attachments();  
 	F = m_pVMDisc->F; 
 R = m_pVMDisc->R; 
 K = m_pVMDisc->temperature(); 
 celsius = m_pVMDisc->temperature_celsius(); 
}  
 
 
 
template<typename TDomain> 
double hh_converted_standard_UG<TDomain>::getgnabar() 
{ 
return gnabar; 
} 
template<typename TDomain> 
double hh_converted_standard_UG<TDomain>::getgkbar() 
{ 
return gkbar; 
} 
template<typename TDomain> 
double hh_converted_standard_UG<TDomain>::getgl() 
{ 
return gl; 
} 
template<typename TDomain> 
double hh_converted_standard_UG<TDomain>::getel() 
{ 
return el; 
} 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::setgnabar(double val) 
{ 
gnabar = val; 
} 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::setgkbar(double val) 
{ 
gkbar = val; 
} 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::setgl(double val) 
{ 
gl = val; 
} 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::setel(double val) 
{ 
el = val; 
} 
 // creating Method for attachments 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::init_attachments() 
{ 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->SGate)) 
UG_THROW("Attachment necessary (SGate) for hh_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->SGate); 
this->aaSGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->SGate); 
 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for hh_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for hh_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
if (spGrid->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for hh_converted_standard_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->nGate); 
 
} 
 
 
 
template<typename TDomain> 
std::vector<number> hh_converted_standard_UG<TDomain>::state_values(number x, number y, number z) 
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
	 if (m_log_SGate || m_log_mGate || m_log_hGate || m_log_nGate )
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
	 	 if (m_log_SGate == true) 
	 	 	 GatingAccesors.push_back(this->aaSGate[bestVrt]); 
	 	 if (m_log_mGate == true) 
	 	 	 GatingAccesors.push_back(this->aamGate[bestVrt]); 
	 	 if (m_log_hGate == true) 
	 	 	 GatingAccesors.push_back(this->aahGate[bestVrt]); 
	 	 if (m_log_nGate == true) 
	 	 	 GatingAccesors.push_back(this->aanGate[bestVrt]); 
	 } 
	 return GatingAccesors; 
} 
 
//Setters for states_outputs 
template<typename TDomain> void hh_converted_standard_UG<TDomain>::set_log_SGate(bool bLogSGate) { m_log_SGate = bLogSGate; }
template<typename TDomain> void hh_converted_standard_UG<TDomain>::set_log_mGate(bool bLogmGate) { m_log_mGate = bLogmGate; }
template<typename TDomain> void hh_converted_standard_UG<TDomain>::set_log_hGate(bool bLoghGate) { m_log_hGate = bLoghGate; }
template<typename TDomain> void hh_converted_standard_UG<TDomain>::set_log_nGate(bool bLognGate) { m_log_nGate = bLognGate; }
 // Init Method for using gatings 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::init(Vertex* vrt, const std::vector<number>& vrt_values) 
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
number na = vrt_values[VMDisc<TDomain>::_na_]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 

 
double           alpha, beta, sum, q10; 
                      //--//Call once from HOC to initialize inf at resting v.
q10= pow(3 , ((celsius-6.3)/10)); 
                //--//"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
                //--//"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
                //--//"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(q10*sum); 
        ninf = alpha/sum; 
aamGate[vrt] = minf; 
aahGate[vrt] = hinf; 
aanGate[vrt] = ninf; 
}  
 
 
 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) 
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
number na = vrt_values[VMDisc<TDomain>::_na_]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 

 
double S = aaSGate[vrt]; 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 
double n = aanGate[vrt]; 

 
 
double           alpha, beta, sum, q10; 
                      //--//Call once from HOC to initialize inf at resting v.; 
q10= pow(3 , ((celsius-6.3)/10)); 
                //--//"m" sodium activation system; 
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
                //--//"h" sodium inactivation system; 
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
                //--//"n" potassium activation system; 
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(q10*sum); 
        ninf = alpha/sum; 
        m  +=   (minf-m)/mtau*dt; 
; 
 
        h  +=  (hinf-h)/htau*dt; 
; 
 
        n  +=  (ninf-n)/ntau*dt; 
; 
 

 
 
aaSGate[vrt] = S; 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
aanGate[vrt] = n; 
 
 
 
} 
 
 
 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
// inits temperatur from kalvin to celsius and some other typical neuron values
number m_T, m_R, m_F; 
m_T = m_pVMDisc->temperature(); 
m_R = m_pVMDisc->R; 
m_F = m_pVMDisc->F; 
 
 
number S = aaSGate[ver]; 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number n = aanGate[ver]; 
number na = vrt_values[m_pVMDisc->_na_]; 
number k = vrt_values[m_pVMDisc->_k_]; 
number v =  vrt_values[m_pVMDisc->_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 

 
 
const number helpV = 1e3*(m_pVMDisc->R*m_pVMDisc->temperature())/m_pVMDisc->F; 
number ena; 
if (m_pVMDisc->ena() == 0) 
{ 
	  ena = helpV*(log(m_pVMDisc->na_out()/na)); 
} 
else 
{ 
	  ena = m_pVMDisc->ena(); 
} 
number ek; 
if (m_pVMDisc->ek() == 0) 
{ 
	  ek = helpV*(log(m_pVMDisc->na_out()/na)); 
} 
else 
{ 
	  ek = m_pVMDisc->ek(); 
} 
 
 
number gna = gnabar*m*m*m*h; 
number gk = gkbar*n*n*n*n; 

 
 
outCurrentValues.push_back( gna*(v - ena) +  gk*(v - ek)       +  gl*(v - el)); 
} 
 
 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::specify_write_function_indices() 
{ 
 
this->m_vWFctInd.push_back(VMDisc<TDomain>::_v_); 
} 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class hh_converted_standard_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class hh_converted_standard_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class hh_converted_standard_UG<Domain3d>; 
#endif 
 
 
} // namespace cable
} // namespace ug


