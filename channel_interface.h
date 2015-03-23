/*
 * channel_interface.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__

#include "lib_grid/lg_base.h"
#include "lib_grid/grid/grid_base_objects.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/interpolate.h"

#include <vector>
#include <stdio.h>

#include "bindings/lua/lua_user_data.h"

#include "common/util/smart_pointer.h"
#include "common/util/vector_util.h"

#include "VM_Disc.h"



namespace ug
{

// forward declaration
template <typename TDomain>
class VMDisc;


template <typename TDomain>
class IChannel
{
	public:
		///	constructor with comma-separated c-string
		IChannel(const char* functions, const char* subsets)
		: m_pVMDisc(NULL)
		{
			m_vSubset = TokenizeString(subsets);
			m_vWFct = TokenizeString(functions);

			// remove white space
			for (size_t i = 0; i < m_vSubset.size(); ++i)
				RemoveWhitespaceFromString(m_vSubset[i]);
			for (size_t i = 0; i < m_vWFct.size(); ++i)
				RemoveWhitespaceFromString(m_vWFct[i]);

			// if no function/subsets passed, clear functions
			if (m_vWFct.size() == 1 && m_vWFct[0].empty()) m_vWFct.clear();
			if (m_vSubset.size() == 1 && m_vSubset[0].empty()) m_vSubset.clear();

			// if functions/subsets passed with separator, but not all tokens filled, throw error
			for (size_t i = 0; i < m_vWFct.size(); ++i)
			{
				if (m_vWFct.empty())
				{
					UG_THROW("Error while setting functions in IChannel: Passed function string lacks\n"
							 "a function specification at position " << i << " (of " << m_vWFct.size()-1 << ")");
				}
			}
			for(size_t i = 0; i < m_vSubset.size(); ++i)
			{
				if (m_vSubset.empty())
				{
					UG_THROW("Error while setting subsets in IChannel: Passed subset string lacks\n"
							 "a subset specification at position " << i << " (of " << m_vSubset.size()-1 << ")");
				}
			}
		};

		/// constructor with vector of string
		IChannel(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
		: m_pVMDisc(NULL)
		{
			m_vSubset = subsets;
			m_vWFct = functions;
		};

		///	destructor
		virtual ~IChannel() {};

		/**
		 * @brief Initializes the defined channel type.
		 *
		 * During the initialization, the necessary attachments are attached to the vertices
		 * and their values calculated by the equilibrium state for the initial membrane potential
		 * given in the grid function passed
		 *
		 * @param time			initial time
		 * @param spGridFct		initial solution (containing membrane potential and ion concentrations)
		 */
		virtual void init(const LocalVector& u, Edge* e) = 0;

		/// updates the gating parameters
		virtual void update_gating(number newTime, const LocalVector& u, Edge* e) = 0;

		/// provides the ionic current (mol*s^-1) at a given vertex
		virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) = 0;

		/// called when access to the root VM disc is possible (i.e. after call to set_vm_disc)
		virtual void vm_disc_available() {};

		/// adding some Jacobian infos at given vertex
		//virtual void Jacobi_sets(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outJFlux) = 0;

		const std::vector<std::string>& write_fcts() {return m_vWFct;}
		const std::vector<std::string>& write_subsets() {return m_vSubset;}
		void set_vm_disc(VMDisc<TDomain>*  vmdisc) {m_pVMDisc = vmdisc; vm_disc_available();}

	protected:
		/// functions whose defect will be written to by this channel
		std::vector<std::string> m_vWFct;

		/// vector of subsets this channel is declared on
		// TODO: somehow transmit subset information to VMDisc (channel must only be considered on this subset!)
		std::vector<std::string> m_vSubset;

		/// joint VMDisc
		VMDisc<TDomain>* m_pVMDisc;
};


template <typename TDomain>
class ChannelHH
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		ChannelHH(const char* functions, const char* subsets)
		try : IChannel<TDomain>(functions, subsets),
		m_g_K(3.6e-4), m_g_Na(1.2e-3), m_g_I(3.0e-6),
		m_rev_pot_K(-74.1266), m_rev_pot_Na(63.5129),
		m_accuracy(1e-12) {}
		UG_CATCH_THROW("Error in ChannelHH initializer list.");

		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		ChannelHH(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
		try : IChannel<TDomain>(functions, subsets),
		m_g_K(3.6e-4), m_g_Na(1.2e-3), m_g_I(3.0e-6),
		m_rev_pot_K(-74.1266), m_rev_pot_Na(63.5129),
		m_accuracy(1e-12) {}
		UG_CATCH_THROW("Error in ChannelHH initializer list.");

		/// destructor
		virtual ~ChannelHH() {};

		/// create attachments and accessors
		void init_attachments();


	/// functions for setting some HH params
		void set_accuracy(double ac);

		void set_conductivities(number Na, number K, number L);

		void set_rev_pot(number R_Na, number R_K);

		number vtrap(number x, number y);

		// inherited from IChannel
		virtual void init(const LocalVector& u, Edge* e);
		virtual void update_gating(number newTime, const LocalVector& u, Edge* e);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vm_disc_available();
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);

	private:
		// membrane conductivities
		number m_g_K;		// C / (m^2 * mV * ms)
		number m_g_Na;		// C / (m^2 * mV * ms)
		number m_g_I;		// C / (m^2 * mV * ms)

		// reversal pot of sodium and potassium
		number m_rev_pot_K;
		number m_rev_pot_Na;

		// params gating
		number m_accuracy;

		// one attachment per state variable
		ADouble m_MGate;
		ADouble m_HGate;
		ADouble m_NGate;

		Grid::AttachmentAccessor<Vertex, ADouble> m_aaMGate;
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaHGate;
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaNGate;
};



template <typename TDomain>
class ChannelHHNernst
	: public IChannel<TDomain>
{
	public:
		using IChannel<TDomain>::m_pVMDisc;

		/// @copydoc IChannel<TDomain>::IChannel(const char*)
		ChannelHHNernst(const char* functions, const char* subsets)
		try : IChannel<TDomain>(functions, subsets),
		m_g_K(3.6e-4), m_g_Na(1.2e-3), m_g_I(3.0e-6),
		m_R(8.314), m_T(310.0), m_F(96485.0),
		m_accuracy(1e-12) {}
		UG_CATCH_THROW("Error in ChannelHHNernst initializer list.");


		/// @copydoc IChannel<TDomain>::IChannel(const std::vector<std::string>&)
		ChannelHHNernst(const std::vector<std::string>& functions, const std::vector<std::string>& subsets)
		try : IChannel<TDomain>(functions, subsets),
		m_g_K(3.6e-4), m_g_Na(1.2e-3), m_g_I(3.0e-6),
		m_R(8.314), m_T(310.0), m_F(96485.0),
		m_accuracy(1e-12) {}
		UG_CATCH_THROW("Error in ChannelHHNernst initializer list.");

		/// destructor
		virtual ~ChannelHHNernst() {};

		/// create attachments and accessors
		void init_attachments();

		// functions for setting some HH params
		void set_accuracy(double ac);

		void set_conductivities(number Na, number K, number L);

		// inherited from IChannel
		virtual void init(const LocalVector& u, Edge* e);
		virtual void update_gating(number newTime, const LocalVector& u, Edge* e);
		virtual void ionic_current(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues);
		virtual void vm_disc_available();
		//virtual void Jacobi_sets(Vertex* vrt, const std::vector<number>& vrt_values, std::vector<number>& outJFlux);


	private:
		// membrane conductivities
		number m_g_K;
		number m_g_Na;
		number m_g_I;

		// constants in Nernst equilibrium
		number m_R;
		number m_T;
		number m_F;

		// params gating
		number m_accuracy;

		// one attachment per state variable
		ADouble m_MGate;
		ADouble m_HGate;
		ADouble m_NGate;

		Grid::AttachmentAccessor<Vertex, ADouble> m_aaMGate;
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaHGate;
		Grid::AttachmentAccessor<Vertex, ADouble> m_aaNGate;
};



} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__
