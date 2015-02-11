#include " haha.h"	
#include " lib_grid/lg_base.h" 
#include " lib_disc/spatial_disc/elem_disc/elem_disc_interface.h" 
#include " lib_disc/function_spaces/grid_function.h" 
#include " lib_disc/function_spaces/local_transfer_interface.h" 
namespace ug { 
 
 
double haha::vtrap(double x, double y, double z, double d) 
{ 
        if (fabs(x/y) < 1e-6) {
                return  y*(1 - x/y/2)
        }else{
                return  x/(exp(x/y) - 1)
        }
}

 
 // Init Method for using gatings 
template<typename TDomain, typename TAlgebra> 
void haha<TDomain, TAlgebra>::init(number time, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) 
{ 
// attach attachments 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->naGate)) 
UG_THROW("Attachment necessary (naGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->naGate); 
this->aanaGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->naGate); 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->kGate)) 
UG_THROW("Attachment necessary (kGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->kGate); 
this->aakGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->kGate); 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->mGate)) 
UG_THROW("Attachment necessary (mGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->mGate); 
this->aamGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->mGate); 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->hGate)) 
UG_THROW("Attachment necessary (hGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->hGate); 
this->aahGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->hGate); 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->nGate)) 
UG_THROW("Attachment necessary (nGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->nGate); 
this->aanGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->nGate); 
 
if (spGridFct->approx_space()->domain()->grid()->has_vertex_attachment(this->vGate)) 
UG_THROW("Attachment necessary (vGate) for Hodgkin and Huxley channel dynamics "
"could not be made, since it already exists."); 
spGridFct->approx_space()->domain()->grid()->attach_to_vertices(this->vGate); 
this->aavGate = Grid::AttachmentAccessor<Vertex, ADouble>(*spGridFct->approx_space()->domain()->grid(), this->vGate); 
 
// creates multiindeces 
std::vector<DoFIndex> multIndna; 
std::vector<DoFIndex> multIndk; 
std::vector<DoFIndex> multIndv; 
 
typedef typename DoFDistribution::traits<Vertex>::const_iterator itType;  
SubsetGroup ssGrp;  
try { ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);}  
UG_CATCH_THROW("Subset group creation failed."); 
  
for (std::size_t si = 0; si < ssGrp.size(); si++)  
{  
itType iterBegin = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);  
itType iterEnd = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
  
for (itType iter = iterBegin; iter != iterEnd; ++iter)  
{  
 
size_t vm = spGridFct->fct_id_by_name("v");  
size_t na = spGridFct->fct_id_by_name("na"); 
size_t k = spGridFct->fct_id_by_name("k"); 

 
spGridFct->dof_indices(*iter, na, multIndna); 
spGridFct->dof_indices(*iter, k, multIndk); 
spGridFct->dof_indices(*iter, vm, multIndv);  
 
 for (size_t i=0; i<multIndv.size(); i++)  
{  
aavGate[*iter] = DoFRef(*spGridFct, multIndv[i]);  
aanaGate[*iter] = DoFRef(*spGridFct, multIndna[i]);  
aakGate[*iter] = DoFRef(*spGridFct, multIndk[i]);  
}  
 
 
double v = aavGate[*iter]; 
 
double           alpha, beta, sum, q10; 
        q10 = 3^((celsius - 6.3)/10); 
        alpha = .1 * vtrap(-(v+40),10); 
        beta =  4 * exp(-(v+65)/18); 
        sum = alpha + beta; 
	mtau = 1/(q10*sum); 
        minf = alpha/sum; 
        alpha = .07 * exp(-(v+65)/20); 
        beta = 1 / (exp(-(v+35)/10) + 1); 
        sum = alpha + beta; 
	htau = 1/(q10*sum); 
        hinf = alpha/sum; 
        alpha = .01*vtrap(-(v+55),10) ; 
        beta = .125*exp(-(v+65)/80); 
	sum = alpha + beta; 
        ntau = 1/(q10*sum); 
        ninf = alpha/sum; 
aamGate[*iter] = minf; 
aahGate[*iter] = hinf; 
aanGate[*iter] = ninf; 
}  
}  
}  
 
 
 
template<typename TDomain, typename TAlgebra> 
void haha<TDomain, TAlgebra>::update_gating(number newTime, SmartPtr<GridFunction<TDomain, TAlgebra> > spGridFct) 
{ 
typedef typename DoFDistribution::traits<Vertex>::const_iterator itType; 
SubsetGroup ssGrp; 
try { ssGrp = SubsetGroup(spGridFct->domain()->subset_handler(), this->m_vSubset);} 
UG_CATCH_THROW("Subset group creation failed."); 
 
for (std::size_t si = 0; si < ssGrp.size(); si++) 
{ 
 
itType iterBegin = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template begin<Vertex>(ssGrp[si]);  
itType iterEnd = spGridFct->approx_space()->dof_distribution(GridLevel::TOP)->template end<Vertex>(ssGrp[si]); 
 
for (itType iter = iterBegin; iter != iterEnd; ++iter) 
{ 
 
Vertex* vrt = *iter; 
// needed konzentration has to be set 
size_t vm = spGridFct->fct_id_by_name("v");  
size_t na = spGridFct->fct_id_by_name("na"); 
size_t k = spGridFct->fct_id_by_name("k"); 

 
spGridFct->dof_indices(*iter, na, multIndna); 
spGridFct->dof_indices(*iter, k, multIndk); 
spGridFct->dof_indices(*iter, vm, multIndv);  
 
 for (size_t i=0; i<multIndv.size(); i++)  
{  
aavGate[*iter] = DoFRef(*spGridFct, multIndv[i]);  
aanaGate[*iter] = DoFRef(*spGridFct, multIndna[i]);  
aakGate[*iter] = DoFRef(*spGridFct, multIndk[i]);  
}  
 
 
double m = aamGate[*iter]; 
double h = aahGate[*iter]; 
double n = aanGate[*iter]; 

 
 
double v = aavGate[*iter]; 
 
        m +=  (minf-m)/mtau*dt; 
        h += (hinf-h)/htau*dt; 
        n += (ninf-n)/ntau*dt; 

 
 
aamGate[*iter] = m; 
aahGate[*iter] = h; 
aanGate[*iter] = n; 
 
 
 
} 
} 
} 
 
 
 
template<typename TDomain, typename TAlgebra> 
void haha<TDomain, TAlgebra>::ionic_current(Vertex* ver, std::vector<number>& outCurrentValues) 
{ 
 
double m = aam[ver]; 
double h = aah[ver]; 
double n = aan[ver]; 
double na = aana[ver];  
double k = aak[ver];  
double v = m_aav[ver]; 
 
 
const number helpV = (m_R*m_T)/m_F; 
number ena = helpV*(log(na_out/na)); 
number ek = helpV*(log(k_out/k)); 
 
 
number gna = gnabar*m*m*m*h; 
number gk = gkbar*n*n*n*n; 

 
 
outCurrentValues.push_back( +  gna*(v - ena) +  gk*(v - ek)       +  gl*(v - el)/m_F ) 
outCurrentValues.push_back( gna*(v - ena)/m_F ) 
outCurrentValues.push_back( gk*(v - ek)      /m_F ) 
  
}  
  
  
