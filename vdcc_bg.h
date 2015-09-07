/*
 * vdcc_bg.h
 *
 *  Created on: 24.08.2015
 *      Author: pgottmann, mbreit
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__VDCC_BG_CABLE_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__VDCC_BG_CABLE_H__


#include "channel_interface.h"

namespace ug {
namespace cable {


template <typename TDomain>
class VDCC_BG_Cable
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		VDCC_BG_Cable(const char* functions, const char* subsets);

		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		VDCC_BG_Cable(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);;

		/// destructor
		virtual ~VDCC_BG_Cable() {};

		/// name
		virtual std::string name();

		/// create attachments and accessors
		void init_attachments();

		/// setter for output behavior of gatings
		void set_log_nGate(bool bLogNGate);
		void set_log_hGate(bool bLogHGate);
		void set_log_mGate(bool bLogMGate);

	public:
		// inherited from IChannel
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vm_disc_available();
		virtual std::vector<number> state_values(number x, number y, number z);

	private:
		virtual void specify_write_function_indices();

	private:
		// membrane conductivities
		int m_mp;
		number m_z_m;
		number m_v12_m;
		number m_tau0_m;

		int m_hp;
		number m_z_h;
		number m_v12_h;
		number m_tau0_h;

		// permeability
		number m_perm;	// in units of m/ms

		// one attachment per state variable
		ANumber m_MGate;
		ANumber m_HGate;

		Grid::AttachmentAccessor<Vertex, ANumber> m_aaMGate;
		Grid::AttachmentAccessor<Vertex, ANumber> m_aaHGate;

		// attachment log_files
		bool m_log_hGate, m_log_mGate;
};



} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__VDCC_BG_CABLE_H__
