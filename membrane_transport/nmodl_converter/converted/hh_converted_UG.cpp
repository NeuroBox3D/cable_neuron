#include "../converted/hh_converted_UG.h"
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain> 
double hh_converted_UG<TDomain>::vtrap(double x, double y)
{ 
        if (fabs(x/y) < 1e-6) {
                return  y*(1 - x/y/2); 
        }else{
                return  x/(exp(x/y) - 1); 
	}
}

 
 
// adding function which always inits_attachments 
template<typename TDomain> 
void hh_converted_UG<TDomain>::ce_obj_available()  
{  
	init_attachments();  
}  
 
 
 
 // creating Method for attachments 
template<typename TDomain> 
void hh_converted_UG<TDomain>::init_attachments() 
{ 
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pCE->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pCE->approx_space()->domain()->grid(); 
if (spGrid->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for hh_converted_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->mGate); 
 
if (spGrid->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for hh_converted_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->hGate); 
 
if (spGrid->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for hh_converted_UG channel dynamics "
"could not be made, since it already exists."); 
spGrid->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGrid, this->nGate); 
 
} 
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void hh_converted_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
{ 
//get celsius 
number celsius = m_pCE->celsius; 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pCE->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t l = 0; l< vl.size(); l++) 
{ 
	 Vertex* vrt = vl[l]; 
 
 
number v = u(m_pCE->_v_, l); 
number na = u(m_pCE->_na_, l); 
number k = u(m_pCE->_k_, l); 

 
double           alpha, beta, sum, q10; 
                      //Call once from HOC to initialize inf at resting v.
q10= pow(3 , ((celsius-6.3)/10)); 
                //"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
                //"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
                //"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(q10*sum); 
        ninf = alpha/sum; 
aamGate[vrt] = minf; 
aahGate[vrt] = hinf; 
aanGate[vrt] = ninf; 
}  
}  
 
 
 
template<typename TDomain> 
void hh_converted_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
{ 
number celsius = m_pCE->celsius; 
 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pCE->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t l = 0; l< vl.size(); l++) 
{ 
	 Vertex* vrt = vl[l]; 
 
 
number dt = newTime - m_pCE->m_aaTime[vrt]; 
number v = u(m_pCE->_v_, l); 
number na = u(m_pCE->_na_, l); 
number k = u(m_pCE->_k_, l); 

 
double m = aamGate[vrt]; 
double h = aahGate[vrt]; 
double n = aanGate[vrt]; 

 
 
double           alpha, beta, sum, q10; 
                      //Call once from HOC to initialize inf at resting v.
q10= pow(3 , ((celsius-6.3)/10)); 
                //"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
                //"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
                //"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(q10*sum); 
        ninf = alpha/sum; 
        m  +=   (minf-m)/mtau*dt; 
        h  +=  (hinf-h)/htau*dt; 
        n  +=  (ninf-n)/ntau*dt; 

 
 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
aanGate[vrt] = n; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void hh_converted_UG<TDomain>::current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number n = aanGate[ver]; 
number na = vrt_values[CableEquation<TDomain>::_na_]; 
number k = vrt_values[CableEquation<TDomain>::_k_]; 
number v =  vrt_values[CableEquation<TDomain>::_v_]; 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number ena = helpV*(log(m_pCE->na_out/na)); 
number ek = helpV*(log(m_pCE->k_out/k)); 
 
 
number gna = gnabar*m*m*m*h; 
number gk = gkbar*n*n*n*n; 

 
 
outCurrentValues.push_back( gna*(v - ena) +  gk*(v - ek)       +  gl*(v - el)); 
 } 
 
 
//////////////////////////////////////////////////////////////////////////////// 
//	explicit template instantiations 
//////////////////////////////////////////////////////////////////////////////// 
#ifdef UG_DIM_1 
template class hh_converted_UG<Domain1d>; 
#endif 
 
 
#ifdef UG_DIM_2 
template class hh_converted_UG<Domain2d>; 
#endif 
 
 
#ifdef UG_DIM_3 
template class hh_converted_UG<Domain3d>; 
#endif 
 
 
}  
  
  
