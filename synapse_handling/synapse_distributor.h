/*
 * synapse_distributor.h
 *
 *  Created on: Mar 24, 2016
 *      Author: lreinhardt
 */

#ifndef UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_DISTRIBUTOR_H_
#define UG__PLUGINS__CABLE_NEURON__SYNAPSE_HANDLING__SYNAPSE_DISTRIBUTOR_H_

// system includes
#include <vector>
#include <string>
#include <utility> // std::pair

// ugcore includes
#include "common/util/smart_pointer.h"
#include "lib_grid/lib_grid.h"
#include "common/types.h"

// synapses
#include "../synapse_handling/synapses/pre_synapse.h"
#include "../synapse_handling/synapses/post_synapse.h"


namespace ug {
namespace cable_neuron {
namespace synapse_handler {


/**
 * This class can be used to distribute primary post-synapses of arbitrary type on a grid,
 * i.e., synapses that are activated by a simple OnsetPreSynapse which is located on the
 * same edge as the post-synapse (not on the axon of another neuron).
**/
class SynapseDistributor
{

private:
	typedef Attachment<std::vector<IBaseSynapse*> > AVSSynapse;

	bool m_bDomBased;

	size_t m_nextID;

	std::string m_InputFile;
	std::string m_LastMessage;

	/// pointers to heap memory, have to be freed in destructor
	SmartPtr<MultiGrid> m_spGrid;
	SmartPtr<MGSubsetHandler> m_spSubsetHandler;
	AVSSynapse m_aSSyn;
	Grid::EdgeAttachmentAccessor<AVSSynapse> m_aaSSyn;
	Grid::VertexAttachmentAccessor<APosition> m_aaPosition;

public:
	///	constructor with file name
	SynapseDistributor(const std::string& infile);

	/// constructor with domain
	SynapseDistributor(SmartPtr<Domain3d> dom);

	/// destructor
	virtual ~SynapseDistributor();

	///	multigrid synapse transfer methods
	void CopySynapsesToAllLevels();
	void CopySynapsesFromParentToChild(Edge* parent);

	/// Removes all synapses from the grid.
	void clear();

	/// Removes all synapses from a given subset.
	void clear(int subsetIndex);

// synapse placement //
	/**
	 * @brief Add a primary post-synapse at a given location.
	 * The synapse is put to the (base-level) edge nearest to the given location.
	 * Local coordinates are calculated accordingly.
	 *
	 * @note As this method loops the complete set of base level edges,
	 * it should only be used for the addition of very few synapses.
	 *
	 * @param vCoords location of synapse
	 * @param pre pre-synapse
	 * @param post post-synapse
	**/
	void place_synapse_at_coords(const std::vector<number>& vCoords, IPreSynapse* pre, IPostSynapse* post);

    /// Place given number of synapses, uniformly distributed.
	void place_synapses_uniform(size_t numSyn, const std::string& type);

    /// Place given number of synapses on edges of given subset, uniformly distributed.
    void place_synapses_uniform(int si, size_t numSyn, const std::string& type);

    /// Place given number of synapses on edges of given subset, uniformly distributed.
    void place_synapses_uniform(const char* subset, size_t numSyn, const std::string& type);

    /// Place synapses only on edges of given subset, uniformly distributed according to set density.
    void place_synapses_uniform(int si, number density, const std::string& type);

    /// Place synapses only on edges of given subset, uniformly distributed according to set density.
    void place_synapses_uniform(const char* subset, number density, const std::string& type);

    /// Place synapses according to given distribution on the subsets.
    void place_synapses(const std::vector<number>& distr, size_t numSynapses, const std::string& type);

    /**
     * @brief Place synapses uniformly on all edges contained within a ball.
     * @param[in] density desired synapse density
     * @param[in] x x coordinate of center
     * @param[in] y y coordinate of center
     * @param[in] z z coordinate of center
     * @param[in] radius radius of the ball
     * @param[in] type synapse type
     */
	void place_synapses_uniform(number density, number x, number y, number z, number radius, const std::string& type);

// synapse degeneration //
	/// Randomly remove a given percentage of all synapses from the grid.
	void degenerate_uniform(number p);

	/// Randomly remove a given percentage of all synapses from a subset (given by index).
	void degenerate_uniform(number p, int si);

    /// Randomly remove a given percentage of all synapses from a subset (given by name).
	void degenerate_uniform(number p, const char* subset);


	/// number of synapses in grid
    size_t num_synapses() const;

    /// number of synapses on subset (given by index)
	size_t num_synapses(int si) const;

    /// number of synapses on subset (given by name)
    size_t num_synapses(const char* subset) const;


    /// Prints the status of the given grid, i.e., subsets and number of vertices in those subsets.
    void print_status();

    /// Save grid with synapses to given file.
    bool export_grid(const std::string& outfile);

protected:
    /// Removes the last synapse (pre and post part) from edge e.
    void remove_synapse(Edge* e);

    /// Removes all synapses (pre and post part) from edge e.
    void remove_all_synapses(Edge* e);

    /// create a pair of OnsetPreSynapse and synapse of type
    void create_synapse_pair(const std::string& type, IPreSynapse** preOut, IPostSynapse** postOut);

    /// Place a single synapse of given type on a specified edge.
    void place_synapse(Edge* e, const std::string& type);

    /// Place given number of synapses on given set of edges, uniformly distributed.
    void place_synapses_uniform(const std::vector<Edge*>& vEdges, size_t numSyn, const std::string& type);

    /// Randomly remove a given number of synapses from given edges.
    void degenerate_uniform(const std::vector<Edge*>& vEdges, size_t numSynapses);

    /// number of synapses in given edge set
    size_t num_synapses(const std::vector<Edge*>& vEdges) const;

    /// number of synapses in given edge set
    size_t max_id() const;
};

} // namespace synapse_handler
} // namespace cable_neuron
} // namespace ug

#endif // UG__PLUGINS__CABLE_NEURON__SYNAPSE_DISTRIBUTOR__SYNAPSE_DISTRIBUTOR_H_
