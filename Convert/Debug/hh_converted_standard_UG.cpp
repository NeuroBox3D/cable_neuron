#include "hh_converted_standard_UG.h"	
#include "lib_grid/lg_base.h" 
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include "lib_disc/function_spaces/grid_function.h" 
#include "lib_disc/function_spaces/local_transfer_interface.h" 
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath> 
namespace ug { 
 
 
template<typename TDomain> 
double hh_converted_standard_UG<TDomain>::vtrap(double x, double y) 
{ 
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
// inits temperatur from kalvin to celsius and some other typical neuron values
m_pVMDisc->celsius = m_T - 273; 
 
 
SmartPtr<Grid> spGrid = m_pVMDisc->approx_space()->domain()->grid(); 
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
 
 
 
 // Init Method for using gatings 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::init(const LocalVector& u, Edge* edge) 
{ 
//get celsius and time
number celsius = m_pVMDisc->celsius; 
number dt = m_pVMDisc->time(); 
// make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t size_l = 0; size_l< vl.size(); size_l++) 
{ 
	 Vertex* vrt = vl[size_l]; 
 
 
number v = u(m_pVMDisc->_v_, size_l); 
number na = u(m_pVMDisc->_na_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 

 
double           alpha, beta, sum, q10; 
                      ;//Call once from HOC to initialize inf at resting v.
q10= pow(3 , ((celsius-6.3)/10)); 
                ;//"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
                ;//"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
                ;//"n" potassium activation system
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
void hh_converted_standard_UG<TDomain>::update_gating(number newTime, const LocalVector& u, Edge* edge) 
{ 
number celsius = m_pVMDisc->celsius; 
 number FARADAY = m_F; 
 // make preparing vor getting values of every edge 
typedef typename MultiGrid::traits<Vertex>::secure_container vrt_list; 
vrt_list vl; 
m_pVMDisc->approx_space()->domain()->grid()->associated_elements_sorted(vl, edge); 
 
 
//over all edges 
for (size_t size_l = 0; size_l< vl.size(); size_l++) 
{ 
	 Vertex* vrt = vl[size_l]; 
 
 
number dt = newTime - m_pVMDisc->m_aaTime[vrt]; 
number v = u(m_pVMDisc->_v_, size_l); 
number na = u(m_pVMDisc->_na_, size_l); 
number k = u(m_pVMDisc->_k_, size_l); 

 
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
	mtau = 1/(sum);
        minf = alpha/sum; 

        //m = (m/dt + alpha)/(1/dt + sum);
        //m = (m + alpha*dt)/ (1 + sum*dt);
        //m = (alpha*dt)/ (1 + sum*dt) + m /(1 + sum*dt);

                //"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(sum);
        hinf = alpha/sum; 

        //h  +=  (hinf-h)/htau*dt;

        //h = (h/dt + alpha)/(1/dt + sum);
        //h = (h + alpha*dt)/ (1 + sum*dt);
        //h = (alpha*dt)/ (1 + sum*dt) + h /(1 + sum*dt);
                //"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(sum);
        ninf = alpha/sum; 

        //n = (alpha*dt)/ (1 + sum*dt) + n /(1 + sum*dt);
        //n = (n + alpha*dt)/ (1 + sum*dt);
        //n = (n/dt + alpha)/(1/dt + sum);

        //n = n / (1 + ntau) + alpha + ninf*dt;


       // u(t+∆t) = u(t) - ∆t*(alpha(V)*u(t) + beta(V)*(1-u(t)))

        //Im impliziten Fall würdest du einfach die beiden rechten u(t) durch u(t+∆t) ersetzen und dann nach u(t+∆t) auflösen, was gehen wird, da das eine lineare Gleichung ist. Da kommt dann typischerweise irgendein Bruch mit (1+∆t*…) im Nenner raus.


       /*n = (n/dt + alpha)/(1/dt + sum);

       n/dt / (1/dt +sum)) + (alpha / (1/dt + sum))
    		   	   	   	   + (alpha/1/dt + alpha/sum)
		n/dt / 1/dt + sum + (alpha/1/dt + ninf) // /*dt
        n / (1 + ntau) + alpha + ninf*dt


        n = (n/dt + alpha)/(1/dt + sum); // *dt

        n = n + alpha*dt / 1 + sum*dt
        n = n / 1 + sum*dt +  alpha*dt/1+sum*dt
		n = n / 1 +sum *dt + alpha/1+sum




		n = (n / 1 + sum*dt) + alpha + ninf


        n = (n + alpha*dt) + (n/sum*dt + alpha*dt/sum*dt) */



        h  +=  (hinf-h)/htau*dt;

        m += (minf-m)/mtau*dt;
 
        n  +=  (ninf-n)/ntau*dt;
        //vector3 test = vertex_position(vrt);
//std::cout << "mau" << std::endl;


typedef typename TDomain::position_type pos_type;
typedef typename TDomain::position_accessor_type AAPos_type;

const AAPos_type aaPos = m_pVMDisc->approx_space()->domain()->position_accessor();
const pos_type coords = aaPos[vrt];

//std::cout << coords << std::endl;

number x=coords[0];
number y=coords[1];
number z=coords[2];


string sh_file = "h_file.txt";
string sm_file = "m_file.txt";
string sn_file = "n_file.txt";

const char* h_file = sh_file.c_str();
const char* m_file = sm_file.c_str();
const char* n_file = sn_file.c_str();

//const char* filenameh = fnameh.c_str();

ofstream myhfile;

ofstream myh_file, mym_file, myn_file;



if (x==0.0 && y==0.0 && z==0.0)
{
	myh_file.open (h_file, std::ios::app);
	myh_file <<  h << "\n";
	myh_file.close();

	mym_file.open (m_file, std::ios::app);
	mym_file <<  m << "\n";
	mym_file.close();

	myn_file.open (n_file, std::ios::app);
	myn_file <<  n << "\n";
	myn_file.close();
}


 
aamGate[vrt] = m; 
aahGate[vrt] = h; 
aanGate[vrt] = n; 
 
 
 
} 
} 
 
 
 
template<typename TDomain> 
void hh_converted_standard_UG<TDomain>::ionic_current(Vertex* ver, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) 
{ 
 
number m = aamGate[ver]; 
number h = aahGate[ver]; 
number n = aanGate[ver]; 
number na = vrt_values[VMDisc<TDomain>::_na_]; 
number k = vrt_values[VMDisc<TDomain>::_k_]; 
number v =  vrt_values[VMDisc<TDomain>::_v_]; 
 
 
number t = m_pVMDisc->time(); 
 
 
const number helpV = 1e3*(m_R*m_T)/m_F; 
number ena; 
if (m_pVMDisc->get_ena() == 0) 
{ 
	  ena = helpV*(log(m_pVMDisc->na_out/na)); 
} 
else 
{ 
	  ena = m_pVMDisc->get_ena(); 
} 
number ek; 
if (m_pVMDisc->get_ek() == 0) 
{ 
	  ek = helpV*(log(m_pVMDisc->k_out/k)); 
} 
else 
{ 
	  ek = m_pVMDisc->get_ek(); 
} 
 
 
number gna = gnabar*m*m*m*h; 
number gk = gkbar*n*n*n*n; 

 
 
outCurrentValues.push_back( gna*(v - ena) +  gk*(v - ek)       +  gl*(v - el)); 
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
 
 
}  
  
  
