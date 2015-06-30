/*
 * channel_interface.h
 *
 *  Created on: 29.10.2014
 *      Author: ppgottmann
 */

#ifndef __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__
#define __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__


#include "VM_Disc.h"


namespace ug {
namespace cable {


// forward declaration
template <typename TDomain>
class VMDisc;


template <typename TDomain>
class IChannel
{
	public:
		///	constructor with comma-separated c-string
		IChannel(const char* functions, const char* subsets);

		/// constructor with vector of string
		IChannel(const std::vector<std::string>& functions, const std::vector<std::string>& subsets);

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
		virtual void init(Vertex* vrt, const std::vector<number>& vrt_values) = 0;

		/// updates the gating parameters
		virtual void update_gating(number newTime, Vertex* vrt, const std::vector<number>& vrt_values) = 0;

		/// provides the ionic current (mol*s^-1) at a given vertex
		virtual void ionic_current(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outCurrentValues) = 0;

		/// called when access to the root VM disc is possible (i.e. after call to set_vm_disc)
		virtual void vm_disc_available() {};

		/// getting values of internal channel states
		virtual std::vector<number> state_values(number x, number y, number z) = 0;

		/// adding some Jacobian infos at given vertex
		//virtual void Jacobi_sets(Vertex* v, const std::vector<number>& vrt_values, std::vector<number>& outJFlux) = 0;

		const std::vector<std::string>& write_fcts() {return m_vWFct;}
		const std::vector<std::string>& write_subsets() {return m_vSubset;}
		void set_vm_disc(VMDisc<TDomain>*  vmdisc) {m_pVMDisc = vmdisc;}

	protected:
		/// functions whose defect will be written to by this channel
		std::vector<std::string> m_vWFct;

		/// vector of subsets this channel is declared on
		std::vector<std::string> m_vSubset;

		/// joint VMDisc
		VMDisc<TDomain>* m_pVMDisc;

};


} // namespace cable
} // namespace ug

#endif // __UG__PLUGINS__EXPERIMENTAL__CABLE__CHANNEL_INTERFACE_H__
