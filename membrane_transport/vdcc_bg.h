/*
 * vdcc_bg.h
 *
 *  Created on: 24.08.2015
 *      Author: pgottmann, mbreit
 */

#ifndef __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__VDCC_BG_H__
#define __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__VDCC_BG_H__


#include "cable_membrane_transport_interface.h"

namespace ug {
namespace cable_neuron {


template <typename TDomain>
class VDCC_BG_cable
	: public ICableMembraneTransport<TDomain>
{
	public:
		using ICableMembraneTransport<TDomain>::m_pCE;

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const char*)
		VDCC_BG_cable(const char* functions, const char* subsets);

		/// @copydoc ICableMembraneTransport<TDomain>::ICableMembraneTransport(const std::vector<std::string>&)
		VDCC_BG_cable(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);;

		/// destructor
		virtual ~VDCC_BG_cable() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// setter for output behavior of gatings
		void set_log_nGate(bool bLogNGate);
		void set_log_hGate(bool bLogHGate);
		void set_log_mGate(bool bLogMGate);

	public:
		// inherited from ICableMembraneTransport
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void ce_obj_available();
		virtual std::vector<number> state_values(number x, number y, number z);

	private:
		virtual void specify_write_function_indices();

	private:
		// membrane conductivities
		int m_mp;
		number m_z_m;
		number m_v12_m;		// in units of mV
		number m_tau0_m;	// in units of ms

		int m_hp;
		number m_z_h;
		number m_v12_h;		// in units of mV
		number m_tau0_h;	// in units of ms

		// permeability
		number m_perm;		// in units of m/s

		// one attachment per state variable
		ANumber m_MGate;
		ANumber m_HGate;

		Grid::AttachmentAccessor<Vertex, ANumber> m_aaMGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaHGate;

		// attachment log_files
		bool m_log_hGate, m_log_mGate;
};



} // namespace cable_neuron
} // namespace ug

#endif // __UG__PLUGINS__CABLE_NEURON__MEMBRANE_TRANSPORT__VDCC_BG_H__
