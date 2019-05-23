/*!
 * \file geometry_information.h
 * \brief struct to store geometry information used during compare in tests
 *
 *  Created on: Jul 20, 2015
 *      Author: stephan
 */

#ifndef __H__UG__NEURONAL_TOPOLOGY_IMPORTER__GEOMETRY_INFORMATION__
#define __H__UG__NEURONAL_TOPOLOGY_IMPORTER__GEOMETRY_INFORMATION__

namespace ug {
namespace cable_neuron {
namespace neuronal_topology_importer {

/*!
 * \brief holds information for geometry we import
 */
struct GeometryInformation {

	private:
		size_t m_vertices;
		size_t m_edges;
		size_t m_exp2Postsyn;
		size_t m_alphaPostsyn;
		size_t m_onsetPresyn;
		size_t m_thresholdPostsyn;

	public:
		GeometryInformation(size_t vertices, size_t edges, size_t nExp2, size_t nAlpha,
			size_t nOnset, size_t nthresh)
		: m_vertices(vertices), m_edges(edges), m_exp2Postsyn(nExp2),
		  m_alphaPostsyn(nAlpha), m_onsetPresyn(nOnset), m_thresholdPostsyn(nthresh)
		{}

		const size_t num_vertices() const {
			return m_vertices;
		}

		const size_t num_edges() const {
			return m_edges;
		}

		const size_t num_exp2syn() const {
			return m_exp2Postsyn;
		}

		const size_t num_alphasyn() const {
			return m_alphaPostsyn;
		}

		const size_t num_onset_presyn() const {
			return m_onsetPresyn;
		}

		const size_t num_threshold_presyn() const {
			return m_thresholdPostsyn;
		}
};

} // namespace neuronal_topology_importer
} // namespace cable_neuron
} // namespace ug

#endif // __H__UG__NEURONAL_TOPOLOGY_IMPORTER__GEOMETRY_INFORMATION__
