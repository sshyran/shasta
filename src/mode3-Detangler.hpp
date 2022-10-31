#ifndef SHASTA_MODE3_DETANGLER_HPP
#define SHASTA_MODE3_DETANGLER_HPP

// Shasta.
#include "ReadId.hpp"

// Standard library.
#include <list>
#include <map>
#include "utility.hpp"
#include "vector.hpp"

/*******************************************************************************

Class mode3::Detangler contains data structures and code used to detangle the
mode3::AssemblyGraph.

In the Detangler, each oriented read is represented by the sequence
of AssemblyGraph segments it visits. This sequence is not necessarily a path
in the AssemblyGraph, unless the assembly graph was created with
minCoverage for links <=1.

This sequence is called a Journey. In the AssemblyGraph,
it is represented as a sequence of AssemblyGraphJourneyEntry objects
and is stored in AssemblyGraph::assemblyGraphJourneys.

In Detangler code, the journey is represented as a sequence of Step objects.
Step(s) are grouped into Cluster(s). All Step(s) in a Cluster refer to
the same segmentId, but there can be more than one Cluster for each segmentId.
At the beginning, there is exactly one Cluster for each segmentId,
but during the detangling process Cluster(s) can be split.

Each Step stores the segmentId it refers to, and an iterator pointing to
the Cluster the Step currently belongs to. The segmentId for a Step never
changes, but the Cluster it points to can change during the detangling process.

*******************************************************************************/

namespace shasta {
    namespace mode3 {
        class Detangler;

        class AssemblyGraph;
    }
}



class shasta::mode3::Detangler {
public:

    // See the comments at the top of this file for the meanings
    // of Step, Journey, Cluster.

    class Cluster;

    class Step {
    public:
        const uint64_t segmentId;
        const Cluster* cluster = 0;

        Step(uint64_t segmentId) :
            segmentId(segmentId) {}
    };

    using Journey = vector<Step>;

    // The journey of each oriented read.
    // Obtained from the AssemblyGraph::assemblyGraphJourneys.
    // Indexed by OrientedReadId::getValue().
     vector<Journey> journeys;

    // Type used to identify a step in a journey.
    class StepInfo {
    public:
        OrientedReadId orientedReadId;

        // The position of this entry in the journey of this oriented read.
        uint64_t position;

        StepInfo() {}
        StepInfo(OrientedReadId orientedReadId, uint64_t position) :
            orientedReadId(orientedReadId),
            position(position) {}
    };

    // A cluster is a set of Step(s) all corresponding to the same
    // segment id.
    class Cluster {
    public:
        uint64_t segmentId;
        uint64_t id = 0;        // Within that segmentId.
        vector<StepInfo> steps; // Sorted by orientedReadId.
        Cluster(uint64_t segmentId, uint64_t id) :
            segmentId(segmentId), id(id) {}
        string stringId() const
        {
            return to_string(segmentId) + "." + to_string(id);
        }
    };

    // Store the clusters keyed by segmentId.
    // Clusters are never removed.
    // However, during detangling, the steps of a cluster
    // can be moved to other clusters for the same segmentId.
    // We use a list so pointers to Cluster(s) are not invalidated
    // when elements are added.
    using ClusterContainer = std::map<uint64_t, std::list<Cluster> >;
    ClusterContainer clusters;


    Detangler(const AssemblyGraph&);
private:
    void createJourneys(const AssemblyGraph&);
    void createInitialClusters();

    // Find the next/previous cluster for each of the steps in a given cluster.
    // The output vector has size equal to the number of steps in this cluster,
    // and the corresponding OrientedReadId(s) are the same
    // as the ones in the steps vector for the given cluster.
    // Some of the pointers returned can be zero. This can happen if this
    // cluster is the first or last cluster in the journey of an oriented read.
    void findNextClusters(
        const Cluster*,
        vector<const Cluster*>&
        ) const;
    void findPreviousClusters(
        const Cluster*,
        vector<const Cluster*>&
        ) const;

    // Simple, classical detangling of a single cluster.
    bool simpleDetangle(Cluster*, uint64_t minLinkCoverage);

    void writeGfa(const string& fileName, uint64_t minLinkCoverage) const;
    void writeGfa(ostream&, uint64_t minLinkCoverage) const;
};



#endif

