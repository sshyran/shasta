#ifndef SHASTA_MODE3_DETANGLER_HPP
#define SHASTA_MODE3_DETANGLER_HPP

// Shasta.
#include "ReadId.hpp"

// Standard library.
#include <list>
#include <map>
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
    };

    // A cluster is a set of Step(s) all corresponding to the same
    // segment id.
    class Cluster {
    public:
        uint64_t segmentId;
        vector<StepInfo> steps;
        Cluster(uint64_t segmentId) :
            segmentId(segmentId) {}
    };

    // Store the clusters keyed by segmentId.
    // Use a list so we can easily remove clusters that are being split.
    using ClusterContainer = std::map<uint64_t, std::list<Cluster> >;
    ClusterContainer clusters;


    Detangler(const AssemblyGraph&);
private:
    void createJourneys(const AssemblyGraph&);
    void createInitialClusters();
};



#endif

