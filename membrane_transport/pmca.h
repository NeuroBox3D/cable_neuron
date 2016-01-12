/*
 * leakage.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann, mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__PMCA_H__
#define __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__PMCA_H__


#include "cable_membrane_transport_interface.h"


namespace ug {
namespace cable_neuron {


template <typename TDomain>
class PMCA_cable
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		PMCA_cable(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		PMCA_cable(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

		/// destructor
		virtual ~PMCA_cable() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		void set_kd(number kd);

		void set_max_flux(number maxFlux);

		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);
		virtual number lin_dep_on_pot(Vertex* vrt, const std::vector<number>& vrt_values);

	private:
		virtual void specify_write_function_indices();

	private:
		number m_kd;					// mM (Elwess et al.)
		//number KD_P = 3.4e-04;		// mM (Graupner)
		number m_maxFlux;				// mol/s

};


} // namespace cable_neuron
} // namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__PMCA_H__
